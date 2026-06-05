/**
 * CP1(b) — invariance test (dynamical, printed verdict + data dump).
 *
 * The Poisson/Ott–Antonsen submanifold is INVARIANT (not attracting) under the
 * sinusoidal mean-field coupling of the two-population Sakaguchi–Kuramoto model
 * (Ott & Antonsen 2008; Pikovsky & Rosenblum 2008). This script verifies that
 * invariance *in this implementation*: start the full integrator from an
 * on-manifold chimera-like IC — both populations are Möbius pushes of a uniform
 * grid, so D = 0 to machine precision at t=0 in the continuum and O(1/N) at
 * finite N — evolve over several breath periods, and confirm D(t) stays at its
 * finite-N floor with no secular growth.
 *
 * This is the physics precondition for CP2–CP4: if the coupling here did NOT
 * preserve the manifold, "distance from the manifold" would not be a meaningful
 * dynamical coordinate. A FAIL here is a STOP condition (do not proceed to CP2).
 *
 * Writes manifold_results/cp1_invariance.jsonl (t, D_incoh, D_sync, R1, R2) and
 * prints the verdict. Deterministic from the committed config.
 *
 * Run:  node tools/manifold-probe/cp1_invariance.mjs
 */
import {
  appendFileSync,
  openSync,
  closeSync,
  mkdirSync,
  readFileSync,
} from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  makeScratch,
  rk4StepInPlace,
  mulberry32,
  TAU,
} from '../chimera-campaign/integrator.mjs';
import {
  mobiusBlock,
  bothPopulations,
  poissonDistance,
  DEFAULT_M,
} from './moments.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'manifold.config.json'), 'utf8'),
);

const M = cfg.observable.M ?? DEFAULT_M;
const { beta, dt, sampleStride } = cfg.model;
const { Np, A, a_sync_mod, a_incoh_mod, horizon_s } = cfg.cp1_invariance;

const mu = (1 + A) / 2;
const nu = (1 - A) / 2;
const alpha = Math.PI / 2 - beta;
const n = 2 * Np;

// On-manifold IC: each population is a Möbius push of N RANDOM-UNIFORM pre-image
// phases — i.e. a representative finite-N *sample* from the Poisson distribution,
// not the measure-zero exact grid (whose D is anomalously ≈0 and whose |a|→1 grid
// correction ∝|a|^N is a finite-grid artifact, not physics). pop1 (sync) |a|→1,
// pop2 (incoherent) moderate |a|. Real-axis a (arg 0) WLOG. Seeded for determinism.
const IC_SEED = 424242;
const rng = mulberry32(IC_SEED);
const psiSync = Array.from({ length: Np }, () => rng() * TAU);
const psiIncoh = Array.from({ length: Np }, () => rng() * TAU);
const phases = new Float64Array(n);
mobiusBlock(phases, 0, Np, { re: a_sync_mod, im: 0 }, psiSync); // pop1
mobiusBlock(phases, Np, Np, { re: a_incoh_mod, im: 0 }, psiIncoh); // pop2

// Empirical finite-N on-manifold floor: mean D over many independent random-uniform
// Möbius pushes at this (N, |a|). This is the realistic noise scale a single
// finite-N sample of the Poisson distribution carries — the trajectory's D should
// fluctuate AROUND this, not climb away from it. Deterministic (seeded).
function floorMC(N, amod, reps, seed) {
  const r = mulberry32(seed);
  let acc = 0;
  for (let k = 0; k < reps; k++) {
    const psi = Array.from({ length: N }, () => r() * TAU);
    const ph = new Float64Array(N);
    mobiusBlock(ph, 0, N, { re: amod, im: 0 }, psi);
    acc += poissonDistance(ph, 0, N, M).D;
  }
  return acc / reps;
}
const floorIncoh = floorMC(Np, a_incoh_mod, 2000, 1234);
const floorSync = floorMC(Np, a_sync_mod, 2000, 5678);

const scratch = makeScratch(n);
const nSteps = Math.round(horizon_s / dt);
const sampleEvery = Math.max(1, Math.round(sampleStride / dt));

const outDir = resolve(process.cwd(), cfg.output_dir);
mkdirSync(outDir, { recursive: true });
const outPath = resolve(outDir, 'cp1_invariance.jsonl');
const fd = openSync(outPath, 'w');

const dIncoh = [];
const dSync = [];
const record = (step) => {
  const b = bothPopulations(phases, Np, M);
  dIncoh.push(b.D_incoh);
  dSync.push(b.D_sync);
  appendFileSync(
    fd,
    JSON.stringify({
      t: step * dt,
      D_incoh: b.D_incoh,
      D_sync: b.D_sync,
      R1: b.R1,
      R2: b.R2,
    }) + '\n',
  );
};

record(0);
for (let step = 1; step <= nSteps; step++) {
  rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch);
  if (step % sampleEvery === 0) record(step);
}
closeSync(fd);

// --- verdict ---------------------------------------------------------------- //
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;
const max = (a) => a.reduce((x, y) => Math.max(x, y), -Infinity);

// Split into halves to test for secular drift (skip the first 10% as a settling
// transient). Compare each population's time-mean to its empirical finite-N floor.
const half = Math.floor(dIncoh.length / 2);
const settle = Math.floor(dIncoh.length * 0.1);
const meanIncoh = mean(dIncoh.slice(settle));
const meanSync = mean(dSync.slice(settle));
const earlyIncoh = mean(dIncoh.slice(settle, half));
const lateIncoh = mean(dIncoh.slice(half));
const earlySync = mean(dSync.slice(settle, half));
const lateSync = mean(dSync.slice(half));
const maxIncoh = max(dIncoh);
const maxSync = max(dSync);

// PASS criteria (per population), all relative to the EMPIRICAL finite-N floor:
//   (i)  time-mean D stays at/below the floor (≤ FLOOR_MULT × floor): the
//        trajectory is no further off-manifold than a random finite-N sample is.
//   (ii) no secular drift: late-half mean within DRIFT× of early-half mean — D
//        fluctuates around a stationary level, it does not climb toward O(1).
//   (iii) D never approaches the O(1) collapse / off-manifold scale.
const FLOOR_MULT = 2.0;
const DRIFT = 3.0;
const COLLAPSE_SCALE = 1.0;
const verdict = (m, e, l, mx, floor) =>
  m <= FLOOR_MULT * floor &&
  l <= DRIFT * Math.max(e, 0.2 * floor) &&
  mx < COLLAPSE_SCALE;
const incohPass = verdict(
  meanIncoh,
  earlyIncoh,
  lateIncoh,
  maxIncoh,
  floorIncoh,
);
const syncPass = verdict(meanSync, earlySync, lateSync, maxSync, floorSync);
const pass = incohPass && syncPass;

console.log('CP1(b) — Poisson-manifold INVARIANCE test');
console.log(
  `On-manifold IC (random-uniform pre-images): Np=${Np}, A=${A}, ` +
    `|a_sync|=${a_sync_mod}, |a_incoh|=${a_incoh_mod}, β=${beta}, dt=${dt}, ` +
    `horizon=${horizon_s}s, M=${M}`,
);
console.log(
  `Empirical finite-N floor (mean D, 2000 random draws): ` +
    `incoh=${floorIncoh.toExponential(2)}, sync=${floorSync.toExponential(2)}`,
);
console.log('');
console.log(
  '  population  D(0)       mean(D)    early½     late½      max(D)     floor      mean/floor  verdict',
);
const fmt = (x) => x.toExponential(2).padStart(9);
console.log(
  `  incoh       ${fmt(dIncoh[0])} ${fmt(meanIncoh)} ${fmt(earlyIncoh)} ${fmt(lateIncoh)} ${fmt(maxIncoh)} ${fmt(floorIncoh)} ${(meanIncoh / floorIncoh).toFixed(2).padStart(10)}  ${incohPass ? 'PASS' : 'FAIL'}`,
);
console.log(
  `  sync        ${fmt(dSync[0])} ${fmt(meanSync)} ${fmt(earlySync)} ${fmt(lateSync)} ${fmt(maxSync)} ${fmt(floorSync)} ${(meanSync / floorSync).toFixed(2).padStart(10)}  ${syncPass ? 'PASS' : 'FAIL'}`,
);
console.log('');
console.log(
  pass
    ? 'CP1(b) GATE: PASS ✅ — D stays at the finite-N floor with no secular growth; ' +
        'the sinusoidal mean-field coupling preserves the Poisson manifold here.'
    : 'CP1(b) GATE: FAIL ❌ — secular growth or floor breach; STOP, diagnose before CP2.',
);
console.log(`Wrote ${outPath}`);
process.exit(pass ? 0 : 1);
