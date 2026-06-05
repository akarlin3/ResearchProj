/**
 * Time-resolved integration with the Poisson-manifold observable.
 *
 * Thin layer over the SHIPPED integrator primitives (`rk4StepInPlace`,
 * `seedChimera`, `mulberry32`, `orderParam`, `makeScratch`) from
 * `../chimera-campaign/integrator.mjs` — it adds nothing to the dynamics, it only
 * also samples D_incoh(t)/D_sync(t)/R₁(t)/R₂(t) on the campaign's model-time
 * stride, and (for CP4) lets the trajectory start from a caller-supplied phase
 * vector instead of the canonical seed. Numerically identical RK4 to the campaign.
 *
 * Two entry points:
 *   traceFromSeed(...)  — seed via seedChimera (matches the campaign exactly),
 *                          integrate to t_max OR until collapse, return full series.
 *   runFromPhases(...)  — start from a supplied phase vector, early-exit on the
 *                          same collapse criterion, return lifetime/censored (the
 *                          CP4 mini-campaign runner for intervention seeds).
 */
import {
  mulberry32,
  seedChimera,
  makeScratch,
  rk4StepInPlace,
  DEFAULT_DT,
} from '../chimera-campaign/integrator.mjs';
import { DEFAULT_THETA, DEFAULT_W } from '../chimera-campaign/collapse.mjs';
import { bothPopulations, DEFAULT_M } from './moments.mjs';

/**
 * Integrate one trajectory (from the canonical seed) and record the manifold
 * observable + order parameters at a model-time stride. Optionally stop a few
 * samples after the collapse criterion fires (so CP3 can event-align without
 * paying full t_max on long survivors).
 *
 * @returns {{ t, D_incoh, D_sync, R1, R2, R_incoh, dt, sampleEvery,
 *             lifetime, censored, collapseIndex }}
 */
export function traceFromSeed({
  Np,
  A,
  beta,
  seed,
  t_max,
  dt = DEFAULT_DT,
  theta = DEFAULT_THETA,
  W = DEFAULT_W,
  M = DEFAULT_M,
  sampleStride = 0.1,
  tailAfterCollapse = 5.0, // s of extra samples to keep past the collapse mark
}) {
  const mu = (1 + A) / 2;
  const nu = (1 - A) / 2;
  const alpha = Math.PI / 2 - beta;
  const n = 2 * Np;

  const phases = seedChimera(Np, mulberry32(seed));
  const scratch = makeScratch(n);

  const nSteps = Math.round(t_max / dt);
  const sampleEvery = Math.max(1, Math.round(sampleStride / dt));
  const sampleDt = sampleEvery * dt;
  const need = Math.max(1, Math.round(W / sampleDt));
  const tailSamples = Math.max(1, Math.round(tailAfterCollapse / sampleDt));

  const t = [];
  const D_incoh = [];
  const D_sync = [];
  const R1 = [];
  const R2 = [];
  const R_incoh = [];

  let consecutive = 0;
  let runStartStep = -1;
  let lifetime = t_max;
  let censored = true;
  let collapseIndex = -1;
  let stopAtIndex = -1;

  const record = (step) => {
    const b = bothPopulations(phases, Np, M);
    t.push(step * dt);
    D_incoh.push(b.D_incoh);
    D_sync.push(b.D_sync);
    R1.push(b.R1);
    R2.push(b.R2);
    R_incoh.push(b.R_incoh);
    const idx = t.length - 1;
    // Collapse criterion on R_incoh (= min(R1,R2)), identical to runner.mjs.
    if (collapseIndex < 0) {
      if (b.R_incoh > theta) {
        if (consecutive === 0) runStartStep = step;
        consecutive++;
        if (consecutive >= need) {
          lifetime = runStartStep * dt;
          censored = false;
          // collapse is dated to the run start; mark the sample index nearest it
          collapseIndex = Math.max(0, idx - (consecutive - 1));
          stopAtIndex = idx + tailSamples;
        }
      } else {
        consecutive = 0;
        runStartStep = -1;
      }
    }
  };

  record(0);
  for (let step = 1; step <= nSteps; step++) {
    rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch);
    if (step % sampleEvery === 0) {
      record(step);
      if (stopAtIndex >= 0 && t.length - 1 >= stopAtIndex) break;
    }
  }

  return {
    t: Float64Array.from(t),
    D_incoh: Float64Array.from(D_incoh),
    D_sync: Float64Array.from(D_sync),
    R1: Float64Array.from(R1),
    R2: Float64Array.from(R2),
    R_incoh: Float64Array.from(R_incoh),
    dt,
    sampleEvery,
    lifetime,
    censored,
    collapseIndex,
  };
}

/**
 * CP4 runner: integrate from a CALLER-SUPPLIED initial phase vector and early-exit
 * on the collapse criterion. Returns the same shape the campaign driver logs, so
 * the intervention mini-campaign reuses the campaign schema. Also returns the
 * initial manifold distances so the actual injected D₀ can be logged alongside the
 * target.
 *
 * @param {Float64Array} phases0  initial phase vector (length 2·Np); copied.
 */
export function runFromPhases({
  phases0,
  Np,
  A,
  beta,
  t_max,
  dt = DEFAULT_DT,
  theta = DEFAULT_THETA,
  W = DEFAULT_W,
  M = DEFAULT_M,
  sampleStride = 0.1,
}) {
  const mu = (1 + A) / 2;
  const nu = (1 - A) / 2;
  const alpha = Math.PI / 2 - beta;
  const n = 2 * Np;
  if (phases0.length !== n) {
    throw new Error(`phases0 length ${phases0.length} != 2*Np=${n}`);
  }
  const phases = Float64Array.from(phases0);
  const scratch = makeScratch(n);

  const b0 = bothPopulations(phases, Np, M);
  const d0_incoh = b0.D_incoh;
  const d0_sync = b0.D_sync;

  const nSteps = Math.round(t_max / dt);
  const sampleEvery = Math.max(1, Math.round(sampleStride / dt));
  const sampleDt = sampleEvery * dt;
  const need = Math.max(1, Math.round(W / sampleDt));

  let consecutive = 0;
  let runStartStep = -1;

  const evalSample = (step) => {
    const bb = bothPopulations(phases, Np, M);
    if (bb.R_incoh > theta) {
      if (consecutive === 0) runStartStep = step;
      consecutive++;
      if (consecutive >= need) return runStartStep * dt;
    } else {
      consecutive = 0;
      runStartStep = -1;
    }
    return null;
  };

  let life = evalSample(0);
  if (life !== null) {
    return { lifetime: life, censored: false, d0_incoh, d0_sync, dt, theta, W };
  }
  for (let step = 1; step <= nSteps; step++) {
    rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch);
    if (step % sampleEvery === 0) {
      life = evalSample(step);
      if (life !== null) {
        return {
          lifetime: life,
          censored: false,
          d0_incoh,
          d0_sync,
          dt,
          theta,
          W,
        };
      }
    }
  }
  return { lifetime: t_max, censored: true, d0_incoh, d0_sync, dt, theta, W };
}
