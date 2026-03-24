# CGNC_v1 — Conformal Gaussian Null Analysis in Higher Dimensions
## Notebook Documentation

**File:** `CGNC_v1.nb`  
**Created with:** Wolfram Mathematica 14.3  
**Backend:** xAct framework (xTras, TexAct, ShowTime1)  
**Date of last computation:** February 2026  

---

## Overview

This notebook demonstrates a fully covariant, index-abstract computation similar to that of 
Bondi–Metzner–Sachs (BMS) analysis in *D* spacetime dimensions, particularly even higher dimensions. It uses the
**xAct** suite for abstract tensor calculus and automated canonicalization, and
carries out the calculation in a *D*-dimensional Conformal Gaussian Null metric ansatz without
fixing the dimension. The main deliverables are:

- Explicit metric components and Christoffel symbols in these coordinates,
- All independent components of the Einstein equations expanded in powers of
  `1/r` to the order needed to extract the BMS asymptotic symmetries,
- The Killing equations for a generic vector field, reduced to their angular and
  radial components,
- The system of equations that defines the higher-dimensional asymptotic symmetry (BMS) group.

---

## 1. Initialization

### xAct Packages Loaded

```mathematica
Block[{Print}, << xAct`xTras`]
Block[{Print}, << xAct`TexAct`]
<< xAct`ShowTime1`
```

| Package | Role |
|---|---|
| `xTras` | Extended tensor algebra on top of xTensor (contraction, symmetrization, canonicalization) |
| `TexAct` | LaTeX output of xAct expressions via `TexPrint` / `TexPrintAlignedEquations` |
| `ShowTime1` | Performance profiling; `$ShowTimeThreshold = 0.1` prints timing for cells taking > 0.1 s |

### Global Settings

```mathematica
$CovDFormat        = "Prefix"       (* ∇_a T instead of T_{;a} *)
$Pre               = ScreenDollarIndices
$CommuteCovDsOnScalars = True
$DefInfoQ          = False           (* suppress def-time messages *)
$ShowTimeThreshold = 0.1
```

`SetOptions[ToCanonical, UseMetricOnVBundle -> None]` prevents the
canonicalizer from automatically raising/lowering indices with the metric,
keeping expressions in the chosen valence throughout.

---

## 2. Manifold and Index Structure

### 2.1 Spacetime Dimension

```mathematica
DefConstantSymbol[Dim, PrintAs -> "D"]
```

`Dim` (printed as *D*) is treated as a symbolic constant throughout, so every
result is valid for arbitrary *D ≥ 4*. No numerical value is ever substituted.

### 2.2 Manifold Decomposition

The full *D*-dimensional spacetime **M** is taken to be a direct product of two
factor manifolds:

| Symbol | Dimension | Index range | Physical role |
|---|---|---|---|
| `M1` | 2 | `a … m` (Latin lower) | Null sector: coordinates (*u*, *r*) |
| `M2` | *D* − 2 | `A, B, C, D, F, G, H` (Latin upper) | Compact (angular) sector: sphere S^{D-2} |
| `M` | *D* = `{M1, M2}` | `α … ν` (Greek) | Full spacetime |

```mathematica
DefManifold[M1, 2,        IndexRange[a, m]]
DefManifold[M2, Dim - 2,  {A, B, C, D, F, G, H}]
DefManifold[M, {M1, M2},  IndexRange[α, ν]]
DefParameter[u];  DefParameter[r]
```

- `u` is the retarded (Bondi) time, `r` is the luminosity distance / radial
  coordinate.
- Indices on `M1` and `M2` are kept *strictly separate* throughout: the
  `indiceOf` / `indicesOf` / `indices3Of` helper functions enforce this when
  building substitution rules.

### 2.3 Covariant Derivatives

| Symbol | Base manifold | Associated metric |
|---|---|---|
| `CD` | M (full spacetime) | `metD` $(g_{αβ})$ |
| `sphCD` | M2 | `metq` $(γ_{AB})$ |
| `PD` | M | Partial derivative ∂ |

Rules encoded in the `"Rules on PD and sphCD"` subsubsection translate
`PD[LI[u]]` and `PD[LI[r]]` into `ParamD[u]` and `ParamD[r]` (parametric
derivatives) so that they can be handled correctly in the `1/r` series.

---

## 3. Tensor Fields (Bondi–Sachs Parametrization)

### 3.1 Metrics Defined

```mathematica
DefMetric[-1, metD[-α, -β], CD,  PrintAs -> "g"]    (* full spacetime metric *)
DefMetric[+1, metq[-A, -B],  sphCD, PrintAs -> "γ"] (* angular metric on M2 *)
DefMetric[+1, met[-A, -B],   …,    PrintAs -> "q"]   (* reference / sphere metric *)
```

`Invmetq[A, B]` $(γ^{AB})$ is the inverse of `metq`.
`Detsphmetq` stores $det(γ_{AB})$ used for the Bondi gauge condition.

### 3.2 Bondi–Sachs Dynamical Fields

All tensors depend on *(M, u, r)* or a subset thereof.

| Name | Symbol | Indices | Printed | Physical meaning |
|---|---|---|---|---|
| `alp` | α | scalar | α | Bondi lapse / mass-aspect function |
| `be` | $β^A$ | `[A]` on M2 | β | Angular shift vector |
| `metq[-A,-B]` | $γ_{AB}$ | symmetric on M2 | γ | Angular (sphere) metric |
| `sphmetq[-A,-B]` | — | symmetric on M2 | — | Reference sphere metric (Bondi gauge) |

### 3.3 Basis and Projection Tensors

These are auxiliary tensors that appear in the metric rules. They are defined
*for completeness* and are not required for the main Einstein-equation
computations.

| Name | Indices | Role |
|---|---|---|
| `uv[α]`, `rv[β]` | `[α]` on M | "Up-index" basis co-vectors for the *u* and *r* directions; vanish on M2 |
| `ud[α]`, `rd[β]` | `[α]` on M | "Down-index" basis vectors for *u* and *r*; vanish on M2 |
| `basv[α, A]` ($= e^α_A$) | `[α]` on M, `[A]` on M2 | Projection from M2 into M (vielbein-like) |
| `basd[α, A]` | `[α]` on M, `[A]` on M2 | Dual projection |

---

## 4. The Bondi–Sachs Metric

### 4.1 Metric Components

The full *D*-dimensional metric is encoded through the replacement rule
`Allmetricrule` (combining `nullmetricrule` and `metricrules`). Its non-zero
components in the coordinate basis *(u, r, $x^A$)* are:

| Component | Value |
|---|---|
| $g_{uu}$ | $2α + γ_{AB} β^A β^B$ |
|$g_{ur} = g_{ru} $| −1 |
| $g_{rr}$ | 0 |
| $g_{uA} = g_{Au}$ | $−r γ_{AB} β^B$ |
| $g_{rA} = g_{Ar}$ | 0 |
| $g_{AB}$ |$ r² γ_{AB}$ |

This is the standard *D*-dimensional Bondi–Sachs metric in retarded
coordinates, with the null condition g_{rr} = 0 and the cross term g_{ur} = −1
normalizing the radial null geodesics.

### 4.2 Inverse Metric (Contravariant Components)

From the above, the non-zero contravariant components are:

| Component | Value |
|---|---|
| $g^{uu}$ | 0 |
| $g^{ur} = g^{ru}$ | −1 |
| $g^{rr}$ | $−2α − γ_{AB} β^A β^B$ |
| $g^{rA} = g^{Ar}$ | $−β^A / r$ |
| $g^{uA} = g^{Au}$| 0 |
| $g^{AB}$ | $γ^{AB} / r²$ |

### 4.3 Bondi Gauge Condition

```mathematica
Detmetq[] := Detsphmetq
```

This enforces $det(γ_{AB})$ = det(sphmetq), the standard Bondi gauge fixing that
eliminates the trace of the angular metric as a dynamical degree of freedom.
The consequence is verified through the Christoffel symbol:

```
Γ^B_{AB} = (1/2) γ^{BC} ∂_A γ_{BC}  =  0
```

under the gauge condition, which is checked explicitly in the **Bondi Gauge
Conditions** subsubsection.

---

## 5. Auxiliary Rules and Helper Functions

### 5.1 Trace Contraction (`TraceNull`, `TraceProductDummy`)

Because indices on `M1` (the null sector) are discrete—they only take values
`u` and `r`—summing over a `TangentM1` index must be done by hand. The
`TraceNull` function replaces a dummy index `a ∈ TangentM1` by the explicit
sum over `{LI[u], LI[r]}`:

```mathematica
TraceNull[expr_] := TraceDummy[expr,
    _?TangentM1`Q :> IndexList[LI[u], LI[r]]]
```

`TraceProductDummy` applies a similar expansion for products. These are called
after projecting Einstein equations onto specific components.

### 5.2 Kronecker Delta Rules (`drule`)

```mathematica
drule = {PD[a_] @ delta[-b_, c_] /; LIndexQ[b] || LIndexQ[c] :> 0}
```

Derivatives of the Kronecker delta that carry a label index (an `LI[u]` or
`LI[r]` index) vanish, since those directions are parametric rather than
coordinate directions in the abstract manifold sense.

### 5.3 `ToCanonicalSym`

A convenience wrapper calling `ToCanonical` after imposing the metric
symmetries. Used after each major simplification step to keep expressions in a
normal form that allows term cancellation.

---

## 6. Series Expansion in 1/r

All dynamical fields are expanded in a formal `1/r` series. The rule `series`
performs this substitution simultaneously on all fields:

```mathematica
series = {
    alp[]       :> Evaluate[expra],
    be[A_]      :> Evaluate[exprb],
    metq[-A,-B] :> Evaluate[exprc],   (* γ_{AB} *)
    ...
}
```

where `expra`, `exprb`, etc., are the truncated `Series[…, {r, ∞, n}]`
expansions. The **truncation order** `n` is set in a dedicated cell labelled
*"The order of truncation : n"* and controls how many subleading terms are
retained in all subsequent computations.

The expansion of α, for example, takes the form:

$$
α(u, r, x^A) = α^(0)(u, x^A) + α^(1)(u, x^A)/r + α^(2)(u, x^A)/r² + O(1/r³)
$$

with the coefficients `α[0]`, `Derivative[1][α][0]`, `Derivative[2][α][0]/2`
appearing as independent data in the asymptotic phase space.

---

## 7. Computing Christoffel Symbols

**Section:** *Computing Christoffel Symbols*

This is the heaviest symbolic section. Christoffel symbols Γ^α_{βγ} of the
full metric `metD` are obtained by the xAct pipeline:

```mathematica
RiemannToChristoffel[RicciCD[-α, -β]] // ChristoffelToGradMetric
```

followed by `Allmetricrule` to substitute metric components, `TraceNull` to
sum over null-sector indices, `ToCanonicalSym` for normal-form reduction, and
`series` for the `1/r` expansion.

Key intermediate objects computed in this section:

- `Christoffelrule` — replacement rules expressing all Christoffel components
  in terms of the Bondi–Sachs fields and their derivatives.
- **Determinant perturbation** (`Perturbation of Determinant`): the variation
  of det(g) appearing in the Einstein–Hilbert action density, used in some
  cross-checks.
- **Deprecated expansions** (marked as such): older versions of some Christoffel
  components kept for reference only.

The Christoffel computation is timed via `ShowTime1` and typically constitutes
the majority of the notebook's runtime.

---

## 8. Einstein Equations

**Section:** *Einstein Equations*

### 8.1 Setup — TEX Configuration

Before computing Ricci components, `TexAct` labels are assigned to all tensors
so that the LaTeX output aligns with standard GR notation:

```mathematica
Tex[Invmetq] = "\\gamma"
```

and analogously for `metD`, `be`, `alp`, etc.

### 8.2 Workflow for Each Component

Every Ricci component R_{μν} is extracted by the same pipeline:

```mathematica
(* 1. Compute Ricci in abstract form *)
(RicciCD[-α, -β] // RiemannToChristoffel // ChristoffelToGradMetric) // Expand;
(* 2. Project onto desired component: e.g. α→LI[r], β→LI[u] *)
% /. {α -> LI[r], β -> LI[u]};
(* 3. Trace over null-sector contractions *)
% // TraceProductDummy // TraceNull;
(* 4. Substitute metric values *)
% //. Allmetricrule /. metricrules;
(* 5. Canonicalize *)
% // ToCanonicalSym;
(* 6. Convert parametric derivatives and store *)
Einsteinru = % /. (PD[-LI[a_]] @ A_) :> (ParamD[a] @ A) // ToCanonical;
(* 7. Apply 1/r series *)
result = Einsteinru /. series // Expand;
```

### 8.3 Components Computed

| Subsection | Component | Stored as |
|---|---|---|
| $R_{uu}$ | Ricci scalar in the uu direction | `Einsteinuu` |
| $R_{ru}$ (TEX config) | The *ru* vacuum equation | `Einsteinru` |
| $R_{AB}$ | Angular–angular equation | `EinsteinAB` |
| $R_{uA}$ | Mixed retarded-time / angular | `EinsteinuA` |

All components are reduced to expressions involving `α`, `β^A`, `γ_{AB}`,
and their u-, r-, and angular derivatives, expanded to the chosen truncation
order in `1/r`.

---

## 9. Killing Vectors and the BMS Group

**Section:** *Killing vectors and BMS group*

### 9.1 Killing Vector Field

```mathematica
DefTensor[ξ[α], {M, r, u}]
```

$ξ^α$ is a generic vector field on the full spacetime. The Killing equation is:

$$
ℒ_ξ g_{βγ} = ∇_β ξ_γ + ∇_γ ξ_β = 0
$$

This is computed in xAct as:

```mathematica
LieD[ξ[α]][metD[-β, -γ]]
// LieDToCovD[#, CD] & // ContractMetric
// ChangeCovD[#, CD, PD] & // ToCanonical
```

yielding the explicit expression in terms of partial derivatives and Christoffel
symbols:


$$∂_β ξ_γ + ∂_γ ξ_β − 2 Γ^α_{βγ} ξ_α = 0$$

### 9.2 Killing Equations (Component Projections)

**Subsection:** *Killing equations*

The abstract Killing tensor equation is projected onto each pair of coordinate
directions by substituting `α → LI[?], β → LI[?]`. The pipeline then applies
`Allmetricrule` and `metricrules` to substitute all metric components, giving
explicit PDEs in *u*, *r*, *x^A*.

The components analyzed are:

| Projection | Physical equation |
|---|---|
| `(α, β) = (LI[r], A)` | Angular shift equation for $ξ^A$ |
| `(α, β) = (LI[r], LI[r])` | Radial equation for $ξ^r$ |
| `(α, β) = (LI[u], A)` | Time-angular equation |
| `(A, B)` | Pure angular Killing equation on M2 |

### 9.3 Lower-Index Killing Equations

**Subsubsection:** *Lower indices*

The covariant (lowered-index) version of the Killing equations is also
extracted, which is often more convenient for matching against the standard BMS
analysis:

$$
ξ_α = g_{αβ} ξ^β
$$

substituted via `Allmetricrule`, giving $ξ_u, ξ_r, ξ_A$ in terms of the
component fields. The resulting system determines the BMS supertranslation
and superrotation generators in *D* dimensions.

---

## 10. Computational Notes

### 10.1 Index Conventions

| Range | Sector | Type |
|---|---|---|
| α, β, γ, δ, ε, ζ, η, θ, ι, κ, λ, μ, ν | Full spacetime M | Greek |
| a, b, c, … m | Null sector M1 | Latin lower |
| A, B, C, D, F, G, H | Angular sector M2 | Latin upper |
| `LI[u]`, `LI[r]` | Label indices for the two null coordinates | — |

Negative index notation (`-α`) denotes a covariant (down) index in xAct
convention.

### 10.2 Known Symbol Clashes

At initialization, xAct raises `ValidateSymbol::capital` warnings for `C` and
`D`, since these overlap with built-in Mathematica symbols. The notebook
resolves this by calling:

```mathematica
Unprotect[D]; Unprotect[C]
PrintAs[C] ^= "C";  PrintAs[D] ^= "D"
```

and `BaseOfVBundle[Labels] ^= {}` to suppress spurious print statements during
canonicalization.

### 10.3 Performance

The `ShowTime1` package monitors each cell. Cells computing Christoffel symbols
and projecting the Einstein tensor are the bottlenecks. The `$ShowTimeThreshold
= 0.1` setting means any cell exceeding 100 ms is timed in the output.

### 10.4 Metric Rule Application Order

The canonical workflow for applying metric values is:

```mathematica
expr //. Allmetricrule /. metricrules
```

`Allmetricrule` (using `//. ` = `ReplaceRepeated`) handles the mixed and
inverse components, which may require multiple passes. `metricrules` (using
`/. ` = `ReplaceAll`) then substitutes the pure angular metric components.

---

## 11. File Structure Summary

```
CGNC_v1.nb
├── Initialization
│   └── xAct load, global options
├── Defining Manifold and Metrics
│   ├── Checks
│   ├── Defining Metrics
│   │   ├── Examples
│   │   └── Rules on PD and sphCD
│   ├── Bondi-Sachs like quantities
│   │   ├── Index matching rules
│   │   └── Examples
│   ├── Defining coordinates in the Null sector
│   │   ├── Metric Values in Null Sector and Compact Sector
│   │   └── Bondi Gauge Conditions
│   ├── Some rules on Determinant
│   │   └── Setting up the Bondi Metric parameters
│   │       └── Examples
│   └── Metric values and Calculations
├── Computing Christoffel Symbols
│   ├── Perturbation of Determinant (used for some derivation)
│   ├── Series Expansion
│   └── Some deprecated expansions
├── Einstein Equations
│   ├── TEX Configuration  [subsection]
│   ├── [uu component]     [subsection — titled via FormBox R_{uu}]
│   ├── ru Component
│   ├── AB component
│   └── uA component
└── Killing vectors and BMS group
    ├── Killing equations
    └── Lower indices
```

---

## 12. Dependencies

| Dependency | Version / Notes |
|---|---|
| Mathematica | 14.3 (Wolfram Language) |
| xAct | xTras, xTensor (included in xTras), xCore |
| TexAct | For LaTeX rendering via `TexPrint` |
| ShowTime1 | Performance timing (custom xAct utility) |

All xAct packages are loaded silently via `Block[{Print}, <<...]` to suppress
the package banner. No external numerical libraries are used.

## 13. TO DO
We also note that this is a quite naive way to implement the computation by involving LabelIndices(which is not usually meant for this use). A more concrete way is to define an implementation
similar to that of Basis, like, 
```mathematica
g[{-a_, -CG}, {-b_, -CG}] instead of g[LI[u],LI[u]] ; for example
```
