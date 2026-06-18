/**
 * beta=0.10 ratchet M_k extraction.
 *
 * Mirrors tools/transient-tests/extract.mjs's CP1 A=0.5 path EXACTLY (same
 * canonical breath detector breath.mjs detectPeaks + movingAverage on the
 * smoothed pre-absorption min(R1,R2) prefix), but applied to the beta=0.10
 * post-homoclinic phase traces. Emits one M_k jsonl per post-homoclinic A.
 *
 * Input : absorption_results/phase_traces_beta010.jsonl
 * Output: transient_results/cp1_mk_beta010_A<A>.jsonl
 * Run   : node tools/beta010-slice/extract_mk_beta010.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { detectPeaks, movingAverage } from '../absorption-recampaign/breath.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const ROOT = process.cwd();
const cfg = JSON.parse(
  readFileSync(
    resolve(__dirname, '../absorption-recampaign/absorption.beta010.config.json'),
    'utf8',
  ),
);
const BREATH = cfg.breath;
const THETA = cfg.graze_criterion.theta;
const TRACES = resolve(ROOT, 'absorption_results/phase_traces_beta010.jsonl');
const OUT = resolve(ROOT, 'transient_results');
mkdirSync(OUT, { recursive: true });

// Post-homoclinic A values (the ratchet/transient analysis points).
const POST_A = cfg.sweeps
  .filter((s) => s.label === 'A_post' || s.label === 'A_depth')
  .map((s) => s.A);

function mkSequence(rPre, sampleDt) {
  const { peaks } = detectPeaks(rPre, sampleDt, BREATH);
  const w = Math.max(1, Math.round(BREATH.smoothWindowSec / sampleDt));
  const sm = movingAverage(rPre, w);
  const Mk = peaks.map((p) => sm[p]);
  const MkRaw = peaks.map((p) => rPre[p]);
  const ptimes = peaks.map((p) => p * sampleDt);
  return { nPeaks: peaks.length, cycles: Math.max(0, peaks.length - 1), peakTimes: ptimes, Mk, MkRaw };
}

const rows = [];
for (const line of readFileSync(TRACES, 'utf8').split('\n')) {
  const s = line.trim();
  if (s) rows.push(JSON.parse(s));
}

for (const A of POST_A) {
  const lines = [];
  let n = 0;
  for (const t of rows) {
    if (Math.abs(t.A - A) > 1e-9) continue;
    if (!(Array.isArray(t.R_incoh) && t.R_incoh.length)) continue;
    if (t.abs_censored) continue;
    const absIdx = t.absIndex > 0 ? t.absIndex : t.R_incoh.length;
    const pre = t.R_incoh.slice(0, Math.min(absIdx, t.R_incoh.length));
    const seq = mkSequence(pre, t.sampleDt);
    lines.push(
      JSON.stringify({
        N: t.N,
        A: t.A,
        seed: t.seed,
        t_abs: t.t_abs,
        abs_censored: t.abs_censored,
        sampleDt: t.sampleDt,
        theta: THETA,
        n_peaks: seq.nPeaks,
        cycles: seq.cycles,
        peak_times: seq.peakTimes,
        Mk: seq.Mk,
        Mk_raw: seq.MkRaw,
      }),
    );
    n++;
  }
  const tag = String(A).replace('.', 'p');
  const outPath = resolve(OUT, `cp1_mk_beta010_A${tag}.jsonl`);
  writeFileSync(outPath, lines.join('\n') + '\n');
  console.log(`A=${A}: wrote ${n} M_k sequences -> ${outPath}`);
}
