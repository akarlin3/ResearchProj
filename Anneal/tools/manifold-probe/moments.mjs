/**
 * CP1 — circular moments + Poisson-manifold distance observable.
 *
 * For a population of N phases {θ_j}, the m-th circular (Kuramoto–Daido) moment is
 *
 *     Z_m = (1/N) Σ_j exp(i·m·θ_j),        m = 1 … M.
 *
 * (Z_1 is the usual complex order parameter R·e^{iΦ}.) The Ott–Antonsen / Poisson
 * submanifold is the set of phase distributions that are a Möbius image of the
 * uniform distribution; on it the moments satisfy the closure identity
 *
 *     Z_m = (Z_1)^m   for all m                         (continuum, exact).
 *
 * Provenance of the identity: for ψ uniform on the circle and the Möbius push
 *     e^{iθ} = (e^{iψ} + a) / (1 + ā·e^{iψ}),   |a| < 1,
 * a contour-integral / residue calculation gives ⟨e^{imθ}⟩ = a^m, hence
 * Z_m = a^m = (Z_1)^m. This is the Watanabe–Strogatz / Ott–Antonsen Poisson
 * kernel structure — see Marvel, Mirollo & Strogatz, Chaos 19, 043104 (2009)
 * (Möbius-group action on the WS constants) and Pikovsky & Rosenblum, PRL 101,
 * 264103 (2008) (partial integrability / OA manifold). The closure is what we
 * VERIFY NUMERICALLY in the CP1 gate; the WS-constant interpretation is left
 * PENDING-READING in the final report.
 *
 * The Poisson-manifold distance of a population is the squared closure defect
 *
 *     D = Σ_{m=2}^{M} | Z_m − (Z_1)^m |².
 *
 * On the submanifold D = 0 in the continuum; at finite N every squared term is a
 * sampling fluctuation of size O(1/N), so D itself is O(1/N) on-manifold.
 *
 * This is a standard circular-moment construction; pure dynamics, no audio path.
 */

import { TAU, orderParam } from '../chimera-campaign/integrator.mjs';

/** Default highest moment order tracked. */
export const DEFAULT_M = 4;

// --- minimal complex helpers (represented as {re, im}) --------------------- //

function cmul(a, b) {
  return { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re };
}

/** Integer power z^p, p ≥ 0, by repeated multiplication (exact for our M≤~8). */
export function cpow(z, p) {
  let r = { re: 1, im: 0 };
  for (let k = 0; k < p; k++) r = cmul(r, z);
  return r;
}

function cabs2(z) {
  return z.re * z.re + z.im * z.im;
}

/**
 * Circular moments Z_1…Z_M of a contiguous block of `count` phases starting at
 * `start`. Returns an array of length M of {re, im} (index k holds Z_{k+1}).
 */
export function circularMoments(phases, start, count, M = DEFAULT_M) {
  const acc = new Array(M);
  for (let m = 0; m < M; m++) acc[m] = { re: 0, im: 0 };
  for (let j = 0; j < count; j++) {
    const th = phases[start + j];
    // accumulate cos/sin of m·θ for m=1..M
    for (let m = 1; m <= M; m++) {
      acc[m - 1].re += Math.cos(m * th);
      acc[m - 1].im += Math.sin(m * th);
    }
  }
  for (let m = 0; m < M; m++) {
    acc[m].re /= count;
    acc[m].im /= count;
  }
  return acc;
}

/**
 * Poisson-manifold distance D = Σ_{m=2}^{M} |Z_m − (Z_1)^m|² for a contiguous
 * block of phases. Also returns the moments and the per-m defect terms so callers
 * can report structure (which harmonic carries the defect).
 *
 * @returns {{ D:number, Z:Array<{re,im}>, terms:number[], R:number }}
 *          terms[k] = |Z_{k+2} − Z_1^{k+2}|²  (k = 0 … M-2); R = |Z_1|.
 */
export function poissonDistance(phases, start, count, M = DEFAULT_M) {
  const Z = circularMoments(phases, start, count, M);
  const z1 = Z[0];
  const terms = [];
  let D = 0;
  for (let m = 2; m <= M; m++) {
    const pred = cpow(z1, m);
    const d = { re: Z[m - 1].re - pred.re, im: Z[m - 1].im - pred.im };
    const t = cabs2(d);
    terms.push(t);
    D += t;
  }
  return { D, Z, terms, R: Math.sqrt(cabs2(z1)) };
}

/**
 * Convenience: D for both populations of a two-population state vector, plus the
 * order parameters, and the role-swap-robust "incoherent" pick (the population
 * with the lower R, matching the campaign collapse criterion). pop1 = [0,Np),
 * pop2 = [Np,2Np).
 *
 * @returns {{ D1, D2, R1, R2, D_incoh, D_sync, R_incoh }}
 */
export function bothPopulations(phases, Np, M = DEFAULT_M) {
  const p1 = poissonDistance(phases, 0, Np, M);
  const p2 = poissonDistance(phases, Np, Np, M);
  const R1 = p1.R;
  const R2 = p2.R;
  // Incoherent = weaker order parameter (role-swap robust). On the seed this is
  // pop2 (the uniform population); after a role swap it can become pop1.
  const incohIsP2 = R2 <= R1;
  return {
    D1: p1.D,
    D2: p2.D,
    R1,
    R2,
    D_incoh: incohIsP2 ? p2.D : p1.D,
    D_sync: incohIsP2 ? p1.D : p2.D,
    R_incoh: Math.min(R1, R2),
  };
}

// --- Möbius push (used by the CP1 identity/invariance constructions) ------- //

/**
 * Möbius push of a phase ψ by complex parameter a (|a| < 1):
 *     e^{iθ} = (e^{iψ} + a) / (1 + ā·e^{iψ}).
 * Returns θ wrapped to [0, 2π). MMS 2009 §II (Möbius-group action).
 *
 * @param {number} psi    pre-image phase
 * @param {{re:number,im:number}} a  Möbius parameter, |a| < 1
 */
export function mobiusPush(psi, a) {
  const w = { re: Math.cos(psi), im: Math.sin(psi) }; // e^{iψ}
  const num = { re: w.re + a.re, im: w.im + a.im }; // e^{iψ} + a
  // 1 + ā·e^{iψ}, ā = conj(a)
  const abar = { re: a.re, im: -a.im };
  const aw = cmul(abar, w);
  const den = { re: 1 + aw.re, im: aw.im };
  const d2 = cabs2(den);
  const q = {
    re: (num.re * den.re + num.im * den.im) / d2,
    im: (num.im * den.re - num.re * den.im) / d2,
  };
  const th = Math.atan2(q.im, q.re);
  return ((th % TAU) + TAU) % TAU;
}

/**
 * Build a Möbius-pushed phase block of length N into `out[start..start+N)`.
 * `grid` true ⇒ uniform grid ψ_j = 2πj/N (on-manifold to machine precision);
 * false ⇒ caller supplies pre-image phases in `psi`.
 */
export function mobiusBlock(out, start, N, a, psi = null) {
  for (let j = 0; j < N; j++) {
    const ps = psi ? psi[j] : (TAU * j) / N;
    out[start + j] = mobiusPush(ps, a);
  }
}

export { orderParam };
