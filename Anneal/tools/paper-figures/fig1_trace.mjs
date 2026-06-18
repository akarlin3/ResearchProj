/**
 * Fig. 1 trace generator — regenerates the chosen example run DETERMINISTICALLY
 * from its logged seed using the SHIPPED-identical RK4 core
 * (chimera-campaign/integrator.mjs), records the per-population order parameters
 * R1, R2 (so the figure can show them faintly behind min(R1,R2)), drives the
 * single-source-of-truth absorption Labeler online so t_graze / t_abs reproduce
 * the campaign bit-for-bit, and replays the shipped supervisor's collapse
 * detector to mark where its 2-s detector would have fired.
 *
 * Adds nothing to the dynamics: same mulberry32 seed, same seedChimera, same
 * rk4StepInPlace, same min(R1,R2) sampling cadence as the campaign. Reads the
 * chosen run + criterion from figures.config.json; writes the trace JSON that
 * fig1.py plots. No prior result file is modified.
 *
 * Usage: node tools/paper-figures/fig1_trace.mjs > paper_figures/fig1_trace.json
 */
import { readFileSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import {
  mulberry32,
  seedChimera,
  orderParam,
  makeScratch,
  rk4StepInPlace,
} from '../chimera-campaign/integrator.mjs';
import { makeLabeler } from '../absorption-recampaign/labeling.mjs';
import { breathPeriod } from '../absorption-recampaign/breath.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'figures.config.json'), 'utf8'),
);
const F = cfg.fig1;

const Np = F.N;
const A = F.A;
const beta = F.beta;
const dt = F.dt;
const seed = F.seed;
const sampleStride = F.sampleStride;
const tMax = F.trace_horizon_s;

const mu = (1 + A) / 2;
const nu = (1 - A) / 2;
const alpha = Math.PI / 2 - beta;
const n = 2 * Np;

const label = {
  theta: F.criterion.theta,
  W: F.criterion.W,
  T_v: F.criterion.T_v,
  recThresh: F.criterion.recoveryThreshold,
  recWin: F.criterion.recoveryWindowSec,
};

const phases = seedChimera(Np, mulberry32(seed));
const scratch = makeScratch(n);

const nSteps = Math.round(tMax / dt);
const sampleEvery = Math.max(1, Math.round(sampleStride / dt));
const sampleDt = sampleEvery * dt;

const labeler = makeLabeler(sampleDt, label);
const R1 = [];
const R2 = [];
const Rincoh = [];
const Rsync = [];

function sampleAt() {
  const r1 = orderParam(phases, 0, Np).R;
  const r2 = orderParam(phases, Np, Np).R;
  R1.push(r1);
  R2.push(r2);
  const lo = r1 < r2 ? r1 : r2;
  const hi = r1 < r2 ? r2 : r1;
  Rincoh.push(lo);
  Rsync.push(hi);
  labeler.push(lo);
}

sampleAt(); // initial condition (step 0)
for (let step = 1; step <= nSteps; step++) {
  rk4StepInPlace(phases, Np, mu, nu, alpha, dt, scratch);
  if (step % sampleEvery === 0) sampleAt();
}

const lab = labeler.result(tMax);

// Breath period over the pre-absorption window (excludes the absorbing cycle),
// identical recipe to the campaign tracer.
const absIdx = lab.abs_censored
  ? Rincoh.length
  : Math.round(lab.t_abs / sampleDt);
const pre = Rincoh.slice(0, Math.min(absIdx, Rincoh.length));
const bp = breathPeriod(pre, sampleDt);

// ---------------------------------------------------------------------------
// Shipped supervisor replay: alive <=> Rsync > SYNC_HI && Rincoh < INCOH_LO;
// a "firing" = not-alive sustained >= COLLAPSE_HOLD_S. Emit firing intervals
// (model-time seconds) so the figure can tick where the 2-s detector triggers.
// ---------------------------------------------------------------------------
const { SYNC_HI, INCOH_LO, COLLAPSE_HOLD_S } = F.supervisor;
const holdNeed = Math.max(1, Math.round(COLLAPSE_HOLD_S / sampleDt));
const firings = [];
{
  // maximal not-alive episodes
  let start = -1;
  const flushEpisode = (s, e) => {
    // e is exclusive end index of the not-alive run
    if (e - s >= holdNeed) {
      // fires at confirmation (s + holdNeed - 1), persists to episode end
      firings.push({
        fire_idx: s + holdNeed - 1,
        fire_t: (s + holdNeed - 1) * sampleDt,
        start_t: s * sampleDt,
        end_t: (e - 1) * sampleDt,
      });
    }
  };
  for (let i = 0; i < Rincoh.length; i++) {
    const alive = Rsync[i] > SYNC_HI && Rincoh[i] < INCOH_LO;
    const notAlive = !alive;
    if (notAlive && start < 0) start = i;
    if (!notAlive && start >= 0) {
      flushEpisode(start, i);
      start = -1;
    }
  }
  if (start >= 0) flushEpisode(start, Rincoh.length);
}

// Breath-peak indices (for the T_b bracket) over the pre-absorption window.
const breathPeaks = (bp.peaks || []).map((p) => ({ idx: p, t: p * sampleDt }));

const out = {
  meta: {
    N: Np,
    A,
    beta,
    dt,
    seed,
    sampleDt,
    tMax,
    theta: label.theta,
    recThresh: label.recThresh,
    W: label.W,
    T_v: label.T_v,
    recWin: label.recWin,
    SYNC_HI,
    INCOH_LO,
    COLLAPSE_HOLD_S,
  },
  labels: {
    t_graze: lab.t_graze,
    graze_censored: lab.graze_censored,
    t_abs: lab.t_abs,
    abs_censored: lab.abs_censored,
    n_grazes_before_abs: lab.n_grazes_before_abs,
    grazeIndex: lab.graze_censored ? -1 : Math.round(lab.t_graze / sampleDt),
    absIndex: lab.abs_censored ? -1 : Math.round(lab.t_abs / sampleDt),
  },
  breath: {
    T_b: bp.Tb,
    n_peaks: bp.nPeaks,
    cycles: bp.cyclesCompleted,
    pAuto: bp.pAuto,
    peaks: breathPeaks,
  },
  supervisor_firings: firings,
  n: Rincoh.length,
  R1,
  R2,
  R_incoh: Rincoh,
  R_sync: Rsync,
};

process.stdout.write(JSON.stringify(out));
