/**
 * CP2 (data) — regenerate the initial manifold distance D₀ for every logged
 * campaign run, directly from its seed. No integration: D₀ depends only on the
 * initial phase vector, which is `seedChimera(N, mulberry32(seed))` — the exact
 * shipped seed construction (verified bit-identical by the campaign crosscheck).
 *
 * For each campaign row we emit the same identifying fields plus:
 *   d0_incoh — D of the seeded INCOHERENT population (pop2, the uniform draw): the
 *              physically relevant distance (the campaign's collapse is the loss of
 *              the incoherent population).
 *   d0_sync  — D of the seeded SYNCHRONIZED population (pop1, tight cluster).
 *   R1_0, R2_0 — initial order parameters (sanity / stratification).
 *
 * Output: manifold_results/cp2_d0.jsonl (one row per campaign row, same order).
 * Determinism: identical seeds + the committed M ⇒ byte-stable output.
 *
 * Run:  node tools/manifold-probe/cp2_d0.mjs
 */
import { readFileSync, writeFileSync, mkdirSync } from 'node:fs';
import { resolve, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { mulberry32, seedChimera } from '../chimera-campaign/integrator.mjs';
import { poissonDistance, DEFAULT_M } from './moments.mjs';

const __dirname = dirname(fileURLToPath(import.meta.url));
const cfg = JSON.parse(
  readFileSync(resolve(__dirname, 'manifold.config.json'), 'utf8'),
);
const M = cfg.observable.M ?? DEFAULT_M;
const campaignPath = resolve(process.cwd(), cfg.inputs.campaign_jsonl);
const outDir = resolve(process.cwd(), cfg.output_dir);
mkdirSync(outDir, { recursive: true });
const outPath = resolve(outDir, 'cp2_d0.jsonl');

const lines = readFileSync(campaignPath, 'utf8')
  .split('\n')
  .filter((l) => l.trim());
console.log(`CP2 D₀: ${lines.length} campaign rows, M=${M}`);

const out = [];
for (const line of lines) {
  const r = JSON.parse(line);
  const Np = r.N;
  const phases = seedChimera(Np, mulberry32(r.seed));
  // pop1 = [0,Np) synchronized; pop2 = [Np,2Np) incoherent (the uniform draw).
  const p1 = poissonDistance(phases, 0, Np, M);
  const p2 = poissonDistance(phases, Np, Np, M);
  out.push(
    JSON.stringify({
      N: r.N,
      A: r.A,
      seed: r.seed,
      lifetime: r.lifetime,
      censored: r.censored,
      d0_incoh: p2.D,
      d0_sync: p1.D,
      R1_0: p1.R,
      R2_0: p2.R,
    }),
  );
}
writeFileSync(outPath, out.join('\n') + '\n');
console.log(`Wrote ${outPath} (${out.length} rows)`);
