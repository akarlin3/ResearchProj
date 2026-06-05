import { describe, expect, it } from 'vitest';
// @ts-expect-error — plain ESM module, no type declarations needed for the probe.
import {
  circularMoments,
  poissonDistance,
  mobiusBlock,
  cpow,
  DEFAULT_M,
} from './moments.mjs';

const TAU = 2 * Math.PI;

/**
 * CP1 validation gate — unit-testable parts.
 *
 *   (a) Identity test: Z_m = (Z_1)^m to machine precision for a Möbius push of a
 *       uniform N-grid (the analytic Poisson submanifold), and D = O(1/N) for a
 *       Möbius push of random-uniform pre-image draws.
 *   (c) Contrast test: a deliberately clustered (bimodal) population must show a
 *       large D (far off the submanifold).
 *
 * The invariance test (b) is a dynamical script (cp1_invariance.mjs), not a unit
 * test. If (a) fails as written, that is a STOP condition for the human to check
 * the Möbius convention against MMS 2009 — the test reports the discrepancy
 * rather than masking it.
 *
 * Identity provenance: ⟨e^{imθ}⟩ = a^m for the Möbius push of the uniform
 * distribution (residue calc); MMS 2009. D-observable: standard circular moments.
 */

function cabs(z: { re: number; im: number }) {
  return Math.hypot(z.re, z.im);
}

describe('CP1(a) — Poisson identity Z_m = (Z_1)^m on the Möbius-pushed grid', () => {
  // Large grid so the discrete grid-correction (∝ |a|^N) is below machine eps for
  // the |a| values tested. N=512, |a|≤0.8 ⇒ 0.8^512 ≈ 1e-50, negligible.
  const N = 512;
  const M = 4;

  for (const aMod of [0.0, 0.2, 0.5, 0.8]) {
    for (const phase of [0, 0.7, 2.1]) {
      it(`|a|=${aMod}, arg(a)=${phase}: moments close to within machine precision`, () => {
        const a = { re: aMod * Math.cos(phase), im: aMod * Math.sin(phase) };
        const phases = new Float64Array(N);
        mobiusBlock(phases, 0, N, a);
        const Z = circularMoments(phases, 0, N, M);
        const z1 = Z[0];
        // Z_1 must equal a (the defining property of the push).
        expect(cabs({ re: z1.re - a.re, im: z1.im - a.im })).toBeLessThan(
          1e-10,
        );
        // Closure: Z_m = z1^m for every m.
        for (let m = 2; m <= M; m++) {
          const pred = cpow(z1, m);
          const dev = cabs({
            re: Z[m - 1].re - pred.re,
            im: Z[m - 1].im - pred.im,
          });
          expect(dev).toBeLessThan(1e-9);
        }
        const { D } = poissonDistance(phases, 0, N, M);
        expect(D).toBeLessThan(1e-16);
      });
    }
  }
});

describe('CP1(a) — D = O(1/N) for a Möbius push of random-uniform pre-images', () => {
  // Deterministic LCG so the test is reproducible.
  function lcg(seed: number) {
    let s = seed >>> 0;
    return () => {
      s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
      return s / 4294967296;
    };
  }

  it('D shrinks ~1/N as N grows (off-grid sampling fluctuation only)', () => {
    const a = { re: 0.4, im: 0.2 };
    const M = DEFAULT_M;
    const rng = lcg(12345);
    const meanDover = (N: number, reps: number) => {
      let acc = 0;
      for (let r = 0; r < reps; r++) {
        const psi = Array.from({ length: N }, () => rng() * TAU);
        const phases = new Float64Array(N);
        mobiusBlock(phases, 0, N, a, psi);
        acc += poissonDistance(phases, 0, N, M).D;
      }
      return acc / reps;
    };
    const dSmall = meanDover(64, 40);
    const dLarge = meanDover(1024, 40);
    // O(1/N): a 16× increase in N should cut D by roughly ~16× (allow slack).
    expect(dLarge).toBeLessThan(dSmall);
    expect(dSmall / dLarge).toBeGreaterThan(4); // clearly decreasing toward 1/N
    // Absolute floor sanity: at N=1024 the closure defect is tiny.
    expect(dLarge).toBeLessThan(0.05);
  });
});

describe('CP1(c) — contrast: a clustered (bimodal) population has large D', () => {
  it('two tight clusters are far off the Poisson submanifold', () => {
    const N = 256;
    const M = DEFAULT_M;
    const phases = new Float64Array(N);
    // Bimodal: half near 0, half near π, tight jitter.
    for (let j = 0; j < N; j++) {
      const center = j < N / 2 ? 0 : Math.PI;
      phases[j] = (((center + 0.05 * (j % 7) - 0.15) % TAU) + TAU) % TAU;
    }
    const { D } = poissonDistance(phases, 0, N, M);
    expect(D).toBeGreaterThan(0.1); // far from the O(1/N) on-manifold floor

    // And it dwarfs an on-manifold population of the same N.
    const a = { re: 0.5, im: 0.0 };
    const onMan = new Float64Array(N);
    mobiusBlock(onMan, 0, N, a);
    const dOn = poissonDistance(onMan, 0, N, M).D;
    expect(D / Math.max(dOn, 1e-18)).toBeGreaterThan(1e6);
  });
});
