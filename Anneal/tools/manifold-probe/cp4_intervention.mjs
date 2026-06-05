/**
 * CP4 (data) — interventional seed families: correlation → causation.
 *
 * The CP2 correlation (D₀ vs lifetime) cannot, on its own, distinguish "manifold
 * distance causes early collapse" from "some hidden seed property drives both."
 * Here we INTERVENE: hold the synchronized population at the canonical chimera
 * envelope and inject a CONTROLLED initial manifold distance into the incoherent
 * population only, then measure lifetime. If lifetime falls monotonically with
 * injected D₀, the manifold-escape mechanism is causally supported; if lifetime is
 * insensitive to injected D₀, it is refuted/incomplete. We report it straight.
 *
 * Construction (per family f, seeded mulberry32(seed0+f), RNG consumed in fixed
 * order so every run is reproducible):
 *   - sync pop (pop1): canonical tight cluster about a random anchor (±jitter),
 *     i.e. the shipped seedChimera envelope — on-manifold, D_sync ≈ 0.
 *   - incoh pop (pop2): base = random-uniform pre-images ψ_j (canonical uniform,
 *     on-manifold), then a controlled m=2 & m=3 harmonic distortion of the
 *     pre-images:  ψ_j → ψ_j + ε·[sin(2ψ_j+φ₂) + sin(3ψ_j+φ₃)]  (φ₂,φ₃ random).
 *     m=2,3 distortions inject exactly the Z₂,Z₃ closure defects that D measures
 *     (MMS 2009: harmonics of the pre-image map to higher OA moments).
 *
 * Levels are ε in the MONOTONE band (config `eps_levels`, default [0,0.5]); each
 * family is run at EVERY level (a PAIRED design: same sync pop + same pre-images +
 * same φ across levels), so a baseline-vs-max paired test removes the dominant
 * family-to-family finite-N floor noise. We drive ε rather than bisecting to an
 * absolute D₀ target because (i) the finite-N sampling floor sets a hard lower
 * bound an absolute target cannot cross, and (ii) D₀(ε) folds over (decreases)
 * above ε≈0.5, so absolute-D₀ bisection is ill-posed there. Realized D₀ rises
 * monotonically over [0,0.5] (~floor → ~2×floor), spanning the upper CP2 range.
 *
 * Output: manifold_results/cp4_intervention.jsonl — campaign schema + `eps`,
 * `d0_target` (per-level nominal D₀; null at ε=0), realized `d0_incoh`, `family`,
 * and the harmonic phases φ₂,φ₃ for full reproducibility.
 *
 * Run:  node tools/manifold-probe/cp4_intervention.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { execSync } from 'node:child_process';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  mulberry32,
  TAU,
  RUNNER_VERSION,
  SEED_JITTER,
} from '../chimera-campaign/integrator.mjs';
import { runFromPhases } from './trace.mjs';
import { poissonDistance, DEFAULT_M } from './moments.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'manifold.config.json'), 'utf8'),
);
const M = cfg.observable.M ?? DEFAULT_M;
const { beta, dt, sampleStride, t_max, theta, W } = cfg.model;
const { Np, A, nFamilies, eps_levels, seed0 } = cfg.cp4;

const wrap = (v) => ((v % TAU) + TAU) % TAU;

const gitHash = (() => {
  try {
    return execSync('git rev-parse --short HEAD', { encoding: 'utf8' }).trim();
  } catch {
    return 'unknown';
  }
})();

// CP2 incoherent-D₀ range (this N,A) — for context only (we parametrize by ε, not
// by an absolute D₀ target, because D₀(ε) is non-monotonic above the band and the
// finite-N floor sets a hard lower bound that absolute targeting cannot cross).
const d0Path = resolve(process.cwd(), cfg.output_dir, 'cp2_d0.jsonl');
const d0rows = readFileSync(d0Path, 'utf8')
  .split('\n')
  .filter((l) => l.trim())
  .map((l) => JSON.parse(l))
  .filter((r) => r.N === Np && r.A === A);
const d0vals = d0rows.map((r) => r.d0_incoh).sort((a, b) => a - b);
const pct = (p) =>
  d0vals.length
    ? d0vals[Math.min(d0vals.length - 1, Math.floor(p * d0vals.length))]
    : NaN;
console.log(
  `CP4 intervention: N=${Np}, A=${A}, ${nFamilies} PAIRED families × ${eps_levels.length} ε-levels`,
);
console.log(
  `CP2 incoherent-D₀ range at this point (context): p10=${pct(0.1).toFixed(3)}, ` +
    `median=${pct(0.5).toFixed(3)}, p95=${pct(0.95).toFixed(3)}`,
);
console.log(`ε-levels (monotone band): ${eps_levels.join(', ')}`);

// --- per-family IC builder -------------------------------------------------- //
// Realized incoherent population for a given ε (m=2 & m=3 pre-image distortion).
function buildIncoh(psi, eps, phi2, phi3) {
  const inc = new Float64Array(psi.length);
  for (let j = 0; j < psi.length; j++) {
    const d = eps * (Math.sin(2 * psi[j] + phi2) + Math.sin(3 * psi[j] + phi3));
    inc[j] = wrap(psi[j] + d);
  }
  return inc;
}

// Per-level nominal D₀ label (cross-family median realized D₀ at that ε): a
// meaningful "d0_target" ladder even though we drive ε, not D₀ directly.
function levelMedianD0(eps) {
  const vals = [];
  for (let f = 0; f < nFamilies; f++) {
    const rng = mulberry32(seed0 + f);
    rng(); // anchor
    for (let j = 0; j < Np; j++) rng(); // sync jitter
    const psi = Array.from({ length: Np }, () => rng() * TAU);
    const phi2 = rng() * TAU;
    const phi3 = rng() * TAU;
    vals.push(poissonDistance(buildIncoh(psi, eps, phi2, phi3), 0, Np, M).D);
  }
  vals.sort((a, b) => a - b);
  return vals[Math.floor(vals.length / 2)];
}
const levelD0 = eps_levels.map((e) => levelMedianD0(e));
console.log(
  `Per-level median realized D₀: ${levelD0.map((x) => x.toFixed(3)).join(', ')}`,
);

const outDir = resolve(process.cwd(), cfg.output_dir);
mkdirSync(outDir, { recursive: true });
const outPath = resolve(
  outDir,
  cfg.cp4.output.replace(/^.*\//, '') || 'cp4_intervention.jsonl',
);
const out = [];

let nRuns = 0;
for (let f = 0; f < nFamilies; f++) {
  const seed = seed0 + f;
  const rng = mulberry32(seed);
  // RNG order (fixed): anchor → Np sync jitter → Np incoh pre-images → φ₂ → φ₃.
  const anchor = rng() * TAU;
  const sync = new Float64Array(Np);
  for (let j = 0; j < Np; j++)
    sync[j] = wrap(anchor + SEED_JITTER * (rng() - 0.5));
  const psi = new Float64Array(Np);
  for (let j = 0; j < Np; j++) psi[j] = rng() * TAU;
  const phi2 = rng() * TAU;
  const phi3 = rng() * TAU;
  const dSync = poissonDistance(sync, 0, Np, M).D;

  // Run the SAME family at every ε-level (paired design).
  for (let li = 0; li < eps_levels.length; li++) {
    const eps = eps_levels[li];
    const inc = buildIncoh(psi, eps, phi2, phi3);
    const phases0 = new Float64Array(2 * Np);
    phases0.set(sync, 0);
    phases0.set(inc, Np);

    const r = runFromPhases({
      phases0,
      Np,
      A,
      beta,
      t_max,
      dt,
      theta,
      W,
      M,
      sampleStride,
    });
    out.push(
      JSON.stringify({
        N: Np,
        A,
        beta,
        dt,
        seed,
        family: f,
        eps,
        d0_target: eps === 0 ? null : levelD0[li], // per-level nominal D₀ (null = ε=0 baseline)
        d0_incoh: r.d0_incoh,
        d0_sync: dSync,
        phi2,
        phi3,
        lifetime: r.lifetime,
        censored: r.censored,
        theta,
        W,
        t_max,
        git_hash: gitHash,
        runner_version: RUNNER_VERSION,
      }),
    );
    nRuns++;
  }
}
writeFileSync(outPath, out.join('\n') + '\n');
console.log(`Wrote ${outPath} (${nRuns} runs)`);
