/**
 * beta=0.10 A_stable control trace.
 *
 * Re-traces the A_stable (A=0.2) seeds at beta=0.10 to t_max WITHOUT early exit
 * (keepTrace), and records the mean R_incoh over the steady (second-half) window
 * for the PERSISTENT (abs_censored) runs — so the finite-N incoherent level can
 * be compared against the reduced stable fixed point r* (the never-absorber
 * control). Absorbed runs are tagged but their steady level is not meaningful.
 *
 * Output: absorption_results/stable_level_beta010.jsonl
 * Run   : node tools/beta010-slice/trace_stable_beta010.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { traceRun } from '../absorption-recampaign/tracer.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = process.cwd();
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, '../absorption-recampaign/absorption.beta010.config.json'), 'utf8'),
);
const { beta, dt, sampleStride, t_max } = cfg.model;
const label = {
  theta: cfg.graze_criterion.theta,
  W: cfg.graze_criterion.W,
  T_v: cfg.absorption_criterion.T_v,
  recThresh: cfg.absorption_criterion.recoveryThreshold,
  recWin: cfg.absorption_criterion.recoveryWindowSec,
};
const breath = cfg.breath;
const sweep = cfg.sweeps.find((s) => s.label === 'A_stable');
const OUT = resolve(ROOT, cfg.output_dir);
mkdirSync(OUT, { recursive: true });

const mean = (a) => a.reduce((s, x) => s + x, 0) / a.length;
const lines = [];
for (const Np of sweep.Ns) {
  for (let s = 0; s < sweep.seeds; s++) {
    const seed = sweep.seed0 + s;
    const r = traceRun({
      Np, A: sweep.A, beta, seed, t_max, dt, sampleStride,
      label, breath, earlyExit: false, keepTrace: true,
    });
    const half = Math.floor(r.R_incoh.length / 2);
    const steady = Array.from(r.R_incoh.slice(half));
    lines.push(JSON.stringify({
      N: Np, A: sweep.A, seed,
      abs_censored: r.abs_censored,
      mean_Rincoh_steady: mean(steady),
      mean_Rincoh_full: mean(Array.from(r.R_incoh)),
    }));
  }
  const persist = lines.filter((l) => JSON.parse(l).N === Np && JSON.parse(l).abs_censored);
  console.log(`A=0.2 N=${Np}: ${persist.length} persistent of ${sweep.seeds}`);
}
const outPath = resolve(OUT, 'stable_level_beta010.jsonl');
writeFileSync(outPath, lines.join('\n') + '\n');
console.log(`Wrote ${outPath} (${lines.length} rows)`);
