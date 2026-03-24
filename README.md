# General EARSM — GCEARSM

A general implementation of an **Explicit Algebraic Reynolds Stress Model (EARSM)** in OpenFOAM 2.3.1.  
The model is called **GCEARSM** (*Generalized Coefficients EARSM*) and extends the standard k-ω SST turbulence model with a data-driven, nonlinear Reynolds-stress correction driven by a user-supplied 84-element coefficient vector **Θ**.

---

## Table of Contents

1. [Installation](#installation)
2. [Repository Structure](#repository-structure)
3. [Physics and Governing Equations](#physics-and-governing-equations)
   - [RANS Equations](#rans-equations)
   - [Constitutive Relation — Reynolds Stress](#constitutive-relation--reynolds-stress)
   - [Nonlinear Correction bΔ — EARSM Ansatz](#nonlinear-correction-b%CE%94--earsm-ansatz)
   - [Turbulent Kinetic Energy Equation (k)](#turbulent-kinetic-energy-equation-k)
   - [Specific Dissipation Rate Equation (ω)](#specific-dissipation-rate-equation-%CF%89)
   - [Eddy Viscosity](#eddy-viscosity)
   - [Blending Functions](#blending-functions)
4. [How to Provide Your Expression (Theta Vector)](#how-to-provide-your-expression-theta-vector)
   - [Theta Index Table](#theta-index-table)
   - [Using the Python Script](#using-the-python-script)
   - [Editing RASProperties Directly](#editing-rasproperties-directly)
5. [Output Fields](#output-fields)
6. [R Correction to the k Equation](#r-correction-to-the-k-equation)
7. [Model Coefficients](#model-coefficients)

---

## Installation

OpenFOAM 2.3.1 can be downloaded here:  
<https://sourceforge.net/projects/openfoam/files/2.3.1/OpenFOAM-2.3.1.tgz/download>

Installation guide for Ubuntu 18.04:  
<https://openfoamwiki.net/index.php/Installation/Linux/OpenFOAM-2.3.1/Ubuntu#Ubuntu_18.04>

Build the custom turbulence model library:

```bash
cd src_turbulenceModels
wmake
```

This compiles `libGCEARSM.so` into `$FOAM_USER_LIBBIN`.

Run the case:

```bash
cd case
simpleFoam
```

Or use the provided Python workflow (see [Using the Python Script](#using-the-python-script)).

---

## Repository Structure

```
.
├── src_turbulenceModels/
│   ├── GCEARSM/
│   │   ├── GCEARSM.C          # Main turbulence model source
│   │   ├── GCEARSM.H          # Class declaration
│   │   └── nonLinearModel.H   # Polynomial expression for bΔ (included at compile time)
│   └── Make/
│       ├── files              # Compilation target list
│       └── options            # Include/library paths
├── case/
│   ├── constant/
│   │   └── RASProperties      # ← Model selection + Theta coefficients go here
│   ├── system/
│   │   └── controlDict, fvSchemes, fvSolution, ...
│   └── 0/                     # Initial conditions
├── documentation/
│   └── info.tex               # LaTeX derivation of the model equations
├── script.py                  # Python helper to write Theta into RASProperties
└── pyfoam/                    # Python utilities to run/post-process OpenFOAM cases
```

---

## Physics and Governing Equations

### RANS Equations

The incompressible, constant-density Reynolds-Averaged Navier-Stokes (RANS) equations read:

```
∂ᵢUᵢ = 0

Uⱼ ∂ⱼUᵢ = ∂ⱼ [ −P/ρ + ν ∂ⱼUᵢ − τᵢⱼ ]
```

where `Uᵢ` is the mean velocity, `P` is the mean pressure, `ρ` is the constant density, `ν` is the kinematic viscosity, and `τᵢⱼ` is the **Reynolds stress** that requires modelling.

The Reynolds stress is decomposed into anisotropic and isotropic parts:

```
τᵢⱼ = 2k ( bᵢⱼ + ¹⁄₃ δᵢⱼ )
```

---

### Constitutive Relation — Reynolds Stress

The standard **Boussinesq** (linear) baseline is:

```
bᵢⱼ⁰ = −(νₜ/k) Sᵢⱼ
```

GCEARSM **augments** this with an additive nonlinear correction `bᵢⱼᐩ` (called `bDelta` in the code):

```
bᵢⱼ = −(νₜ/k) Sᵢⱼ  +  bᵢⱼᐩ
       └─────────────┘   └─────┘
        linear (Boussinesq)   nonlinear correction
```

The full Reynolds stress field computed by the solver is therefore:

```
τᵢⱼ = (2/3) k δᵢⱼ  −  2 νₜ Sᵢⱼ  +  2k bᵢⱼᐩ
```

This additive approach leaves the linear (Boussinesq) part of the solver matrix **unchanged**, giving enhanced numerical stability.

---

### Nonlinear Correction bΔ — EARSM Ansatz

Following Pope (1975) and the Cayley-Hamilton theorem, the most general symmetric, traceless tensor formed from `Sᵢⱼ` and `Ωᵢⱼ` can be written as a polynomial in three **basis tensors** with **scalar invariant** coefficients.

A local timescale is defined as:

```
τ = 1 / max( S/a₁ + ωₘᵢₙ,  ω + ωₘᵢₙ )
```

#### Basis Tensors

| Symbol | Formula | Code variable |
|--------|---------|---------------|
| T¹ᵢⱼ | τ Sᵢⱼ | `T1 = tau * Sij` |
| T²ᵢⱼ | τ² symm(Sᵢₖ Ωₖⱼ − Ωᵢₖ Sₖⱼ) | `T2 = tau² * symm(S·Ω − Ω·S)` |
| T³ᵢⱼ | τ² symm(Sᵢₖ Sₖⱼ) − ¹⁄₃ δᵢⱼ I₁ | `T3 = tau² * symm(S·S) − ¹⁄₃ I * l1` |

where `Sᵢⱼ = dev(symm(∇U))` and `Ωᵢⱼ = −½(∇U − ∇Uᵀ)`.

#### Scalar Invariants

| Symbol | Formula | Code variable |
|--------|---------|---------------|
| I₁ (l1) | τ² tr(SᵢⱼSⱼᵢ) | `l1 = tau² * tr(S·S)` |
| I₂ (l2) | τ² tr(ΩᵢⱼΩⱼᵢ) | `l2 = tau² * tr(Ω·Ω)` |

#### The General Polynomial Expression

The nonlinear correction `bᵢⱼᐩ` is expressed as a polynomial of degree up to D=6 in I₁ and I₂, multiplied by each basis tensor.  
The 28 polynomial terms per basis tensor are:

```
D=0:                     1
D=1:              I₁,          I₂
D=2:         I₁²,     I₁I₂,       I₂²
D=3:    I₁³,    I₁²I₂,   I₁I₂²,      I₂³
D=4:  I₁⁴,  I₁³I₂, I₁²I₂², I₁I₂³,     I₂⁴
D=5: I₁⁵, I₁⁴I₂, I₁³I₂², I₁²I₂³, I₁I₂⁴,  I₂⁵
D=6: I₁⁶, I₁⁵I₂, I₁⁴I₂², I₁³I₂³, I₁²I₂⁴, I₁I₂⁵, I₂⁶
```

28 terms × 3 basis tensors = **84 coefficients** — the vector **Θ**.

The complete expression implemented in `nonLinearModel.H`:

```
bᵢⱼᐩ = T¹ᵢⱼ · (Θ[0]·1 + Θ[1]·I₁ + Θ[2]·I₂ + Θ[3]·I₁² + ... + Θ[27]·I₂⁶)
      + T²ᵢⱼ · (Θ[28]·1 + Θ[29]·I₁ + Θ[30]·I₂ + ... + Θ[55]·I₂⁶)
      + T³ᵢⱼ · (Θ[56]·1 + Θ[57]·I₁ + Θ[58]·I₂ + ... + Θ[83]·I₂⁶)
```

---

### Turbulent Kinetic Energy Equation (k)

```
∂k/∂t + Uⱼ ∂k/∂xⱼ = min(Pₖ, c₁ β* ω k)  −  β* ω k
                     + ∂/∂xⱼ [ (ν + σₖ νₜ) ∂k/∂xⱼ ]
```

The **augmented production** term includes the nonlinear correction:

```
Pₖ = νₜ S²  −  2k (bᵢⱼᐩ : Sᵢⱼ)
     └──────┘   └─────────────────┘
      linear      nonlinear correction (PkDelta)
```

In code:
```cpp
PkDelta_ = 2.0*k_*(bDelta_) && symm(tgradU());   // 2k bΔ : S
Pk_      = nut_*S2 - PkDelta_;                    // Pₖ = νₜ S² − PkDelta
```

---

### Specific Dissipation Rate Equation (ω)

```
∂ω/∂t + Uⱼ ∂ω/∂xⱼ = (γ/νₜ) · min(Pₖ/(νₜ+ε), c₁/a₁ · β* ω · max(a₁ω, b₁ F₂₃ S))
                     − β ω²
                     + ∂/∂xⱼ [ (ν + σω νₜ) ∂ω/∂xⱼ ]
                     + CD_{kω}
```

where the cross-diffusion term is:

```
CD_{kω} = max( 2 σω₂ (1/ω) ∂k/∂xᵢ ∂ω/∂xᵢ,  10⁻¹⁰ )
```

---

### Eddy Viscosity

```
νₜ = a₁ k / max( a₁ ω,  b₁ F₂₃ S )

S = sqrt(2 Sᵢⱼ Sᵢⱼ)
```

---

### Blending Functions

The standard k-ω SST blending functions are used unchanged:

```
arg₁ = min( min( max( √k/(β* ω y),  500ν/(y² ω) ),  4σω₂ k/(CD_{kω} y²) ),  10 )
F₁   = tanh(arg₁⁴)

arg₂ = min( max( 2√k/(β* ω y),  500ν/(y² ω) ),  100 )
F₂   = tanh(arg₂²)
```

Coefficients are blended: `Φ = F₁·Φ₁ + (1−F₁)·Φ₂`

---

## How to Provide Your Expression (Theta Vector)

You supply `bᵢⱼᐩ` by filling in the 84-element `Theta` list in `case/constant/RASProperties`.

### Theta Index Table

Each group of 28 entries in `Theta` corresponds to one basis tensor.  
Within each group the ordering of polynomial terms is:

| Local index | Term | Global index T¹ | Global index T² | Global index T³ |
|:-----------:|------|:---------------:|:---------------:|:---------------:|
| 0  | 1              |  0 | 28 | 56 |
| 1  | I₁             |  1 | 29 | 57 |
| 2  | I₂             |  2 | 30 | 58 |
| 3  | I₁²            |  3 | 31 | 59 |
| 4  | I₁I₂           |  4 | 32 | 60 |
| 5  | I₂²            |  5 | 33 | 61 |
| 6  | I₁³            |  6 | 34 | 62 |
| 7  | I₁²I₂          |  7 | 35 | 63 |
| 8  | I₁I₂²          |  8 | 36 | 64 |
| 9  | I₂³            |  9 | 37 | 65 |
| 10 | I₁⁴            | 10 | 38 | 66 |
| 11 | I₁³I₂          | 11 | 39 | 67 |
| 12 | I₁²I₂²         | 12 | 40 | 68 |
| 13 | I₁I₂³          | 13 | 41 | 69 |
| 14 | I₂⁴            | 14 | 42 | 70 |
| 15 | I₁⁵            | 15 | 43 | 71 |
| 16 | I₁⁴I₂          | 16 | 44 | 72 |
| 17 | I₁³I₂²         | 17 | 45 | 73 |
| 18 | I₁²I₂³         | 18 | 46 | 74 |
| 19 | I₁I₂⁴          | 19 | 47 | 75 |
| 20 | I₂⁵            | 20 | 48 | 76 |
| 21 | I₁⁶            | 21 | 49 | 77 |
| 22 | I₁⁵I₂          | 22 | 50 | 78 |
| 23 | I₁⁴I₂²         | 23 | 51 | 79 |
| 24 | I₁³I₂³         | 24 | 52 | 80 |
| 25 | I₁²I₂⁴         | 25 | 53 | 81 |
| 26 | I₁I₂⁵          | 26 | 54 | 82 |
| 27 | I₂⁶            | 27 | 55 | 83 |

> **Note:** The solver enforces that `Theta` contains **exactly 84 entries**.  
> A `FatalError` is raised otherwise.

---

### Using the Python Script

`script.py` provides a helper to build and inject the coefficient vector:

```python
from script import create_theta, change_coef

# Specify which Theta indices are nonzero and their values
theta_ind = [0, 10, 22, 25]         # active term indices
theta_val = [0.1, 0.2, 0.3, 0.9]   # their coefficient values

# Build the dense 84-element vector (all other entries = 0)
theta_dense = create_theta(theta_ind, theta_val)
# Theta[0]=0.1  →  bΔ += 0.1 · T¹
# Theta[10]=0.2 →  bΔ += 0.2 · I₁⁴ · T¹
# Theta[22]=0.3 →  bΔ += 0.3 · I₁⁵I₂ · T¹
# Theta[25]=0.9 →  bΔ += 0.9 · I₁²I₂⁴ · T¹

# Write into case/constant/RASProperties
change_coef('case', theta_dense)

# Run the solver
import pyfoam
pyfoam.run_case('case', 'simpleFoam')
```

---

### Editing RASProperties Directly

Open `case/constant/RASProperties` and fill in the 84 values under `GCEARSMCoeffs`:

```c++
RASModel   GCEARSM;
turbulence on;
printCoeffs on;

GCEARSMCoeffs
{
    // Standard k-omega SST coefficients (defaults shown)
    alphaK1     0.85;
    alphaK2     1.0;
    alphaOmega1 0.5;
    alphaOmega2 0.856;
    beta1       0.075;
    beta2       0.0828;
    betaStar    0.09;
    gamma1      0.5556;   // 5/9
    gamma2      0.44;
    a1          0.31;
    b1          1.0;
    c1          10.0;
    F3          no;

    // 84-element EARSM coefficient vector
    // Theta[0..27]  → T¹ terms
    // Theta[28..55] → T² terms
    // Theta[56..83] → T³ terms
    Theta (
        0.1    // [0]  T¹ · 1
        0.0    // [1]  T¹ · I₁
        0.0    // [2]  T¹ · I₂
        0.0    // [3]  T¹ · I₁²
        // ... (must be exactly 84 entries)
        0.0    // [83] T³ · I₂⁶
    );
}
```

---

## Output Fields

The model writes the following fields to each time directory:

| Field     | Type        | Description |
|-----------|-------------|-------------|
| `k`       | Scalar      | Turbulent kinetic energy |
| `omega`   | Scalar      | Specific dissipation rate |
| `nut`     | Scalar      | Eddy viscosity: `a₁k / max(a₁ω, b₁F₂₃S)` |
| `bDelta`  | SymmTensor  | Nonlinear correction `bᵢⱼᐩ` |
| `aDelta`  | SymmTensor  | `aᵢⱼᐩ = 2k bᵢⱼᐩ` |
| `bij`     | SymmTensor  | Full normalised anisotropy: `bᵢⱼ = −(νₜ/k)Sᵢⱼ + bᵢⱼᐩ` |
| `aij`     | SymmTensor  | Full Reynolds stress anisotropy: `aᵢⱼ = 2k bᵢⱼ` |
| `Rall`    | SymmTensor  | Full Reynolds stress: `τᵢⱼ = (2/3)k δᵢⱼ + aᵢⱼ` |
| `Pk`      | Scalar      | Augmented production: `νₜ S² − 2k (bᵢⱼᐩ : Sᵢⱼ)` |
| `PkDelta` | Scalar      | Nonlinear production correction: `2k (bᵢⱼᐩ : Sᵢⱼ)` |

---

## R Correction to the k Equation

The predecessor model (`dataDrivenKOmegaSST_modKEq`) used an **explicit** separate field `R` added as a source term to the k equation.  

In GCEARSM this concept is **absorbed into the augmented production term** `Pₖ`:

```
Pₖ = νₜ S²  −  2k (bᵢⱼᐩ : Sᵢⱼ)
               └──────────────────┘
                   ≡ PkDelta (the R-correction)
```

The term `−PkDelta` is precisely the R-correction to the k equation.  
- If `bᵢⱼᐩ` is **co-aligned** with `Sᵢⱼ` (positive double contraction) → production is **reduced**.  
- If `bᵢⱼᐩ` is **anti-aligned** with `Sᵢⱼ` (negative double contraction) → production is **increased**.

This allows the nonlinear model to locally correct the turbulent kinetic energy by changing its net production rate, without requiring a separate transport equation for `R`.

---

## Model Coefficients

Default k-ω SST coefficients (can be overridden in `GCEARSMCoeffs`):

| Coefficient | Default | Description |
|-------------|---------|-------------|
| `alphaK1`   | 0.85    | σₖ inner zone |
| `alphaK2`   | 1.0     | σₖ outer zone |
| `alphaOmega1` | 0.5   | σω inner zone |
| `alphaOmega2` | 0.856 | σω outer zone |
| `beta1`     | 0.075   | β inner zone |
| `beta2`     | 0.0828  | β outer zone |
| `betaStar`  | 0.09    | β* |
| `gamma1`    | 5/9     | γ inner zone |
| `gamma2`    | 0.44    | γ outer zone |
| `a1`        | 0.31    | Bradshaw constant |
| `b1`        | 1.0     | |
| `c1`        | 10.0    | Production limiter |
| `F3`        | no      | Rough-wall correction (Hellsten 1998) |

---

## References

- Menter, F., Esch, T., *"Elements of Industrial Heat Transfer Prediction"*, COBEM 2001.
- Menter, F. R., Kuntz, M., Langtry, R., *"Ten Years of Industrial Experience with the SST Turbulence Model"*, THMT-4, 2003.
- Pope, S. B., *"A more general effective-viscosity hypothesis"*, J. Fluid Mech., 1975.
- Gatski, T. B., Jongen, T., *"Nonlinear eddy viscosity and algebraic stress models"*, Progress in Aerospace Sciences, 2000.
- Weatheritt, J., Sandberg, R., *"A novel evolutionary algorithm applied to algebraic modifications of the RANS stress–strain relationship"*, JCP, 2016.
- Hellsten, A., *"Some Improvements in Menter's k-omega-SST turbulence model"*, AIAA-98-2554, 1998.
