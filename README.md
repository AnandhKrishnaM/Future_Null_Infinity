# Future Null Infinity — Higher Dimensional Analysis

In this repository, you will find the codes written for analysis of higher
dimensional asymptotically flat spacetime in the **xAct** framework for
Mathematica.  This work was carried out during a Masters thesis and is
provided as **version 1** of the analysis; alternate approaches exist, but the
current implementation is sufficient for the cases studied.

---

## Contents

| File | Description |
|------|-------------|
| `HigherDim_NullInfinity_v1.nb` | Main Mathematica notebook — version 1 of the analysis |

---

## Overview

The notebook `HigherDim_NullInfinity_v1.nb` implements the analysis of
**D-dimensional (D ≥ 4) asymptotically flat spacetimes** using the Bondi–Sachs
formalism.  The key topics covered are:

1. **Manifold and tensor setup** — D-dimensional spacetime manifold, Lorentzian
   metric, and the round metric on the transverse (D−2)-sphere, all defined
   with xTensor.
2. **Bondi–Sachs metric** — the general D-dimensional Bondi–Sachs line element
   in retarded coordinates (u, r, x^A) together with the determinant condition
   on the angular metric h_{AB}.
3. **Asymptotic falloff conditions** — systematic 1/r expansion of the metric
   functions β, V, U^A, h_{AB} appropriate for D ≥ 4.
4. **Curvature computation** — Riemann, Ricci and Einstein tensors computed
   symbolically via xTensor, with Bianchi identity verified automatically.
5. **Einstein equations in Bondi gauge** — decomposition into hypersurface,
   evolution and supplementary equations; Bondi mass-loss and angular
   momentum-loss formulae extracted at leading order.
6. **Generalised BMS symmetry** — supertranslations and (for D = 4)
   superrotations; for D > 4 only the global Lorentz group SO(D−1,1) acts as
   the superrotation part.  Killing equation on S^{D−2} coded explicitly.
7. **Conserved charges** — Bondi mass Q[f] and angular momentum J[Y] expressed
   as integrals over future null infinity (ℐ⁺) via the Wald–Zoupas formalism.
8. **Gravitational memory** — displacement memory tensor ΔC_{AB} = ∫ N_{AB} du
   and its scalar measure.
9. **xPert linearised perturbation theory** — first-order perturbations of the
   Ricci and Einstein tensors around flat Minkowski space; transverse-traceless
   gauge conditions set up.
10. **Concrete 5-dimensional example** — explicit background Bondi metric,
    3-sphere angular metric, leading-order Christoffel symbols.

---

## Prerequisites

- **Mathematica** 12.0 or later (tested on 13.x).
- **xAct** suite (free, open-source):
  - `xTensor` — abstract tensor algebra
  - `xCoba`  — component-based (coordinate) calculations
  - `xPert`  — perturbation theory

  Install from [http://www.xact.es/](http://www.xact.es/) or via the
  Mathematica Package Manager.

---

## Usage

1. Install the xAct packages so that Mathematica can find them on its
   `$Path`.
2. Open `HigherDim_NullInfinity_v1.nb` in Mathematica.
3. Adjust the spacetime dimension at the top of **Section 1**:
   ```mathematica
   D = 5;  (* change to any integer >= 4 *)
   ```
4. Evaluate all cells (Evaluation → Evaluate Notebook).

---

## Physical Background

The future null infinity **ℐ⁺** of an asymptotically flat spacetime is the
boundary to which outgoing null geodesics converge.  In D = 4 it carries the
celebrated **BMS₄ symmetry** — an infinite-dimensional group comprising
supertranslations and superrotations — which underlies Weinberg's soft graviton
theorem and the gravitational-wave memory effect.

For **D > 4** the situation changes qualitatively:

- The asymptotic falloff is faster (metric perturbations decay as r^{−(D−3)/2}
  rather than r^{-1} in D = 4).
- Only a **finite-dimensional** asymptotic symmetry group survives:
  supertranslations ⋊ SO(D−1,1).
- Logarithmic soft factors and the infinite tower of super-Lorentz charges
  present in D = 4 are absent.
- The memory effect is still present but its amplitude is suppressed.

---

## Key Equations

**Bondi–Sachs metric (D dimensions)**

```
ds² = −(e^{2β} V/r − r² h_{AB} U^A U^B) du²
      − 2 e^{2β} du dr
      − 2 r² h_{AB} U^B du dx^A
      + r² h_{AB} dx^A dx^B
```

**Bondi mass-loss formula**

```
∂_u M = − 1 / (2(D−2))  N_{AB} N^{AB}  +  (angular divergence terms)
```

**Displacement memory**

```
ΔC_{AB} = ∫_{−∞}^{+∞} N_{AB} du
```

---

## References

- Bondi, van der Burg, Metzner (1962) *Proc. R. Soc. London* A269, 21
- Sachs (1962) *Phys. Rev.* 128, 2851
- Hollands & Ishibashi (2005) *Commun. Math. Phys.* 261, 1–7
- Kapec, Raclariu & Strominger (2015) [arXiv:1507.04352](https://arxiv.org/abs/1507.04352)
- Campiglia & Laddha (2020) [arXiv:2011.01908](https://arxiv.org/abs/2011.01908)
- xAct framework: [http://www.xact.es/](http://www.xact.es/)

---

## Status

This is **version 1** — it covers the foundational setup and the leading-order
analysis.  Planned extensions for v2 include:

- Full subleading 1/r expansion of all Bondi functions.
- Explicit Bondi mass-loss formula with angular-momentum flux.
- Relation between D > 4 soft theorems and absence of soft divergences.
- Treatment of odd-dimensional spacetimes (half-integer falloff rates).
- Numerical cross-checks of the analytic expressions.
