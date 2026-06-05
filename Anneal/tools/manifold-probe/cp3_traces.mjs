/**
 * CP3 (data) — time-resolved escape traces.
 *
 * For each configured (N, A) point, pick `nSeeds` campaign seeds that SPAN the
 * lifetime distribution (quantile-stratified, so we get short, medium and long
 * survivors), then re-integrate each with the manifold observable, recording
 * D_incoh(t), D_sync(t), R₁(t), R₂(t), R_incoh(t) on the campaign stride until the
 * collapse criterion fires (+ a short tail). The Python side event-aligns these at
 * t_collapse to test whether manifold escape PRECEDES the order-parameter collapse.
 *
 * Output: manifold_results/cp3_traces.jsonl — one row per run with the full series
 * plus lifetime, collapseIndex (sample index nearest t_collapse) and the seed.
 *
 * Determinism: seeds are read from the committed campaign + selected by a fixed
 * quantile rule; the integration is the shipped RK4.
 *
 * Run:  node tools/manifold-probe/cp3_traces.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { traceFromSeed } from './trace.mjs';
import { DEFAULT_M } from './moments.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'manifold.config.json'), 'utf8'),
);
const M = cfg.observable.M ?? DEFAULT_M;
const { beta, dt, sampleStride, t_max, theta, W } = cfg.model;
const { points, nSeeds } = cfg.cp3;

const campaignPath = resolve(process.cwd(), cfg.inputs.campaign_jsonl);
const rows = readFileSync(campaignPath, 'utf8')
  .split('\n')
  .filter((l) => l.trim())
  .map((l) => JSON.parse(l));

const outDir = resolve(process.cwd(), cfg.output_dir);
mkdirSync(outDir, { recursive: true });
const outPath = resolve(outDir, 'cp3_traces.jsonl');

/** Quantile-stratified seed pick: sort by lifetime, take nSeeds evenly spaced. */
function pickSeeds(N, A, k) {
  const pool = rows
    .filter((r) => r.N === N && r.A === A)
    .sort((a, b) => a.lifetime - b.lifetime);
  if (pool.length <= k) return pool.map((r) => r.seed);
  const out = [];
  for (let i = 0; i < k; i++) {
    const idx = Math.round((i * (pool.length - 1)) / (k - 1));
    out.push(pool[idx].seed);
  }
  return out;
}

const lines = [];
const round = (x) => Number(x.toFixed(6));
for (const { Np, A } of points) {
  const seeds = pickSeeds(Np, A, nSeeds);
  console.log(
    `CP3: N=${Np}, A=${A} — ${seeds.length} quantile-stratified seeds`,
  );
  for (const seed of seeds) {
    const tr = traceFromSeed({
      Np,
      A,
      beta,
      seed,
      t_max,
      dt,
      theta,
      W,
      M,
      sampleStride,
      tailAfterCollapse: 5.0,
    });
    lines.push(
      JSON.stringify({
        N: Np,
        A,
        seed,
        lifetime: tr.lifetime,
        censored: tr.censored,
        collapseIndex: tr.collapseIndex,
        sampleDt: tr.sampleEvery * tr.dt,
        t: Array.from(tr.t, round),
        D_incoh: Array.from(tr.D_incoh, round),
        D_sync: Array.from(tr.D_sync, round),
        R1: Array.from(tr.R1, round),
        R2: Array.from(tr.R2, round),
        R_incoh: Array.from(tr.R_incoh, round),
      }),
    );
  }
}
writeFileSync(outPath, lines.join('\n') + '\n');
console.log(`Wrote ${outPath} (${lines.length} traces)`);
