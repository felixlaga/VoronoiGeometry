# DEBT — open items, deferred work, known limitations

Rules of this file (CLAUDE.md): deferred work lives here, never in TODO
comments. Every item is either open, or closed with the commit that closed it.
No silent disappearances.

Seeded at T0 from IMPLEMENTATION_SPEC Sec. 6; updated as tasks complete.

- Tangential enumeration is complete for **Bravais lattices only**. Non-lattice
  structures are covered only by the FCC sublattice scan to `k <= 20`;
  non-sublattice perfect coverings, HCP-derived depletions and general
  multi-orbit structures are not enumerated.
- `(k-1) | 12` observed, not derived — to be CLOSED at T2 by the counting
  identity + spectral bound (`geometry/coloring.py`), with independence forced
  by tangentiality.
- Whether uniform-K degeneracy is EXTENSIVE (entropy density > 0) is open;
  attack by `enumerate_uniform_K` on growing tori (T4). Highest-value theory
  item.
- Existence of uniform-K states on HCP is open (spectral condition passes for
  all K<=4); if they exist, density degenerates across stackings and topology
  classification becomes mandatory, not optional.
- No closed-form 3D analogue of `q_r(n)` exists; any generalised face number is
  a documented choice. Default: do not use one; work with `Q_iso` directly.
- `eps*` is composition-independent but has NOT been shown to resolve the
  transition; the 2D control was consistent with a featureless power law.
- `f_c` at fixed tolerance is unmeasurable (exact contact has measure zero) and
  its `f_c(delta)` repair is degenerate with the median gap.
- MC dynamics is a Brownian proxy; `eta_a` from it is not the EDMD `eta_a`.
- `g(sigma+)` linear extrapolation degrades above `eta ~ 0.4` (4% at 0.45).
- `eps* = C (eta_m - eta)^p` is the wrong functional form near contact
  percolation; window spread 0.15. A free-volume-motivated form is needed.
- Whether any of these structures is thermodynamically relevant at `T > 0` is
  entirely open, in 3D as in 2D.
