/**
 * CP3 helper — measure the finite-N A=0.2 never-absorbers' mean R_incoh level.
 *
 * The A=0.2 persistent (abs_censored) runs were never trace-dumped, so we
 * re-trace a deterministic sample with the SHIPPED-identical tracer (the exact
 * same traceRun the transient-tests CP1 contrast used) and report the mean and
 * median of min(R1,R2) over the steady window. This is the finite-N counterpart
 * of the reduced stable-chimera fixed point r* — read-only inputs, deterministic.
 *
 * Usage: node tools/reduced-ode/trace_a02_level.mjs > reduced_results/cp3_a02_level.json
 */
import { readFileSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { traceRun } from '../absorption-recampaign/tracer.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = resolve(__dirname, '../..');
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'reduced.config.json'), 'utf8'),
);

const campaignPath = resolve(ROOT, 'absorption_results/absorption_campaign.jsonl');
const rows = readFileSync(campaignPath, 'utf8')
  .split('\n')
  .filter((s) => s.trim())
  .map((s) => JSON.parse(s));

const beta = cfg.beta_ours; // 0.05
const Ns = [8, 16, 32, 64];
const perN = 20; // deterministic sample (lowest-id persistent seeds)
const tmax = 2000;
const dt = cfg.timescale.dt;
const sampleStride = cfg.timescale.sampleStride;

// Settle window: drop the first 200s, average the steady remainder.
const settleSec = 200;

const label = {
  theta: cfg.boundary.theta,
  W: 5.0,
  T_v: 120.0,
  recThresh: cfg.boundary.recoveryThreshold,
  recWin: cfg.boundary.recoveryWindowSec,
};

const perRun = [];
for (const N of Ns) {
  const persistent = rows
    .filter((r) => r.A === 0.2 && r.N === N && r.abs_censored)
    .sort((a, b) => a.seed - b.seed)
    .slice(0, perN);
  for (const r of persistent) {
    const tr = traceRun({
      Np: N,
      A: 0.2,
      beta,
      seed: r.seed,
      t_max: tmax,
      dt,
      label,
      sampleStride,
      earlyExit: false,
      keepTrace: true,
    });
    const series = tr.R_incoh;
    const i0 = Math.round(settleSec / tr.sampleDt);
    const steady = series.slice(i0);
    if (steady.length < 10) continue;
    const sorted = [...steady].sort((a, b) => a - b);
    const mean = steady.reduce((s, x) => s + x, 0) / steady.length;
    const median = sorted[Math.floor(sorted.length / 2)];
    perRun.push({ N, seed: r.seed, mean, median, max: sorted[sorted.length - 1] });
  }
}

const agg = {};
for (const N of Ns) {
  const sub = perRun.filter((p) => p.N === N);
  const means = sub.map((p) => p.mean);
  const meds = sub.map((p) => p.median);
  agg[N] = {
    n: sub.length,
    mean_Rincoh: means.reduce((s, x) => s + x, 0) / means.length,
    median_Rincoh: meds.reduce((s, x) => s + x, 0) / meds.length,
  };
}
const allMeans = perRun.map((p) => p.mean);
const pooled = allMeans.reduce((s, x) => s + x, 0) / allMeans.length;

console.log(
  JSON.stringify(
    {
      note: 'Finite-N A=0.2 persistent (never-absorber) mean R_incoh over the steady window (>200s), re-traced with the shipped-identical tracer. Compare to the reduced stable-chimera fixed point r*.',
      beta,
      Ns,
      perN,
      settleSec,
      per_N: agg,
      pooled_mean_Rincoh: pooled,
      n_runs: perRun.length,
    },
    null,
    2,
  ),
);
