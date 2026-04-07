# Arruda-Boyce Hyperelastic Material Model in Scilab

> **Advanced Material Characterization & FEM Implementation**  
>  
> A comprehensive Scilab implementation of the Arruda-Boyce hyperelastic constitutive model, featuring experimental data fitting, analytical tangent stiffness computation, and 3D anisotropy visualization for large-strain elastomers.

---

## Table of Contents

- [Overview](#overview)
- [Project Structure](#project-structure)
- [Key Features](#key-features)
- [Installation & Setup](#installation--setup)
- [Core Modules](#core-modules)
- [Material Characterization Workflow](#material-characterization-workflow)
- [Documentation](#documentation)
- [License](#license)

---

## Overview

This project implements a **complete hyperelastic material characterization pipeline** based on the Arruda-Boyce model, a state-of-the-art constitutive framework for elastomers (rubber-like materials). The model captures the characteristic S-shaped stress-strain response and strain-induced stiffening (locking effect) observed in real rubber under large deformations.

### What is Arruda-Boyce?

The Arruda-Boyce model provides an energy-based description of hyperelastic materials through:

$$W(\mathbf{C}) = W_{\text{dev}}(\bar{I}_1) + W_{\text{vol}}(J)$$

Where:
- **$W_{\text{dev}}$**: Deviatoric (shape-changing) energy based on an inverse Langevin function
- **$W_{\text{vol}}$**: Volumetric (compression-resisting) energy

This dual-potential structure accurately reproduces both the nonlinear material stiffening at high strains and the incompressibility typical of elastomers.

---

## Project Structure

```
arruda_boyce/
├── README.md                             # This file
├── arruda_boyce_erklärung.md             # Complete theoretical documentation (German)
├── Herleitung_Tangensteifigkeit[...].md  # Detailed derivations (German)
│
├── 00_General                      # Aufgabenstellung, Organisation in der Gruppe (gitignored)
│
├── 01_Doku/                        # Dokumentation, Bericht (gitignored)
│    └── Bilder/...
│
├── 02_Code/                        # (gitignored) 
│   ├── Main.sci                    # Primary workflow orchestrator, Plotting
│   ├── Material_Model.sci          # Core material routine (F → σ)
│   ├── Curve_Fit.sci               # Parameter optimization (lsqrsolve)
│   ├── Tangentialsteifigkeit.sci   # Tangent stiffness (E_T) & glyphs
│   ├── [...]_Test.txt              # Extracted Experimental stress-strain data from Hyperelastic.txt
│   └── ...
│
├── 03_References/                     # (gitignored)
│   ├── 01_Unterlagen von Prof/        # Lecture notes & reference implementations
│   │   ├── ET-NeoHook.sci             # 3D glyph example
│   │   └── ...
│   └── Hyperelastic.txt               # Experimental stress-strain data
└── ...
```

---

## Key Features

### 1. **Material Routine (`F_to_sig`)**
   - Converts deformation gradient **F** to Cauchy stress **σ**
   - Implements Taylor-expanded inverse Langevin function for efficiency
   - Separates deviatoric and volumetric contributions

### 2. **Multi-Mode Loading Simulations**
   - **Uniaxial tension**: Standard tensile test (free lateral contraction)
   - **Biaxial stress**: Membrane inflation scenario
   - **Pure shear**: Optimal for parameter identification
   - Automatic correction for pressure and reference-to-nominal stress transformation

### 3. **Advanced Parameter Fitting**
   - Levenberg-Marquardt optimization (`lsqrsolve`) for μ, β, K identification
   - **Combined fitting** across all three load cases for robust, physically meaningful parameters
   - Automatic initial guess generation
   - Detailed residual analysis

### 4. **Analytical Tangent Stiffness (E_T)**
   - Full 4th-order tensor computation via exact partial derivatives
   - Both material (Lagrangian) and spatial (Eulerian) formulations
   - Push-forward transformation from reference to current configuration
   - Zero approximation error (unlike finite differences)

### 5. **3D Anisotropy Visualization**
   - Directional stiffness glyphs: radius = E_nnnn(n) in direction n
   - Side-by-side comparison of undeformed (isotropic sphere) vs. deformed (anisotropic ellipsoid)
   - Automatic colorbar scaling with adaptive formatting
   - Isoview perspective for physical proportionality preservation

### 6. **Comprehensive Documentation**
   - Theoretical foundations (German, detailed)
   - Derivation of key equations with intermediate steps
   - Code-to-math correspondence table
   - Complete nomenclature and symbol reference

---

## Installation & Setup

### Requirements
- **Scilab** ≥ 2025.1.0 (or compatible version)
- **Text editor** with UTF-8 support (for .md files)
- **Data file**: `Hyperelastic.txt` with experimental measurements

### Setup Steps

1. **Clone or download this repository**
   ```bash
   git clone https://github.com/bochsimon-sketch/arruda_boyce.git
   cd arruda_boyce
   ```

---

## 🔧 Core Modules

### `02_Material_Model.sci`
- **`F_to_sig(F, mu, beta, K)`**: Main routine (F → σ)
- **`calculate_taylor_expansion_coefficients(beta)`**: Langevin Taylor series
- Computes first invariant: $\bar{I}_1 = J^{-2/3} \text{tr}(\mathbf{C})$

### `03_Simulation_Routines.sci`
| Function | Load Case | Deformation Gradient | 
|----------|-----------|----------------------|
| `simulate_uniaxial_tension()` | Uniaxial | $\text{diag}(\lambda, \lambda^{-1/2}, \lambda^{-1/2})$ |
| `simulate_biaxial_tension()` | Biaxial | $\text{diag}(\lambda, \lambda, \lambda^{-2})$ |
| `simulate_pure_shear()` | Shear | $\text{diag}(\lambda, 1, \lambda^{-1})$ |

### `04_Curve_Fit.sci`
- **`fit_combined_data(exp_data, initial_guess)`**
  - Levenberg-Marquardt solver via `lsqrsolve`
  - All three load cases simultaneously
  - Returns: [μ*, β*, K*], residuals, convergence info

### `05_Tangentialsteifigkeit.sci`
- **`compute_material_stiffness_tensor(F, mu, beta, K)`**: $\mathbb{C}_T = 2\frac{\partial \mathbf{P}}{\partial \mathbf{C}}$
- **`push_forward_to_spatial(C_T, F, J)`**: Transform to spatial stiffness
- **`visualize_stiffness_glyph_3d(...)`**: 3D rendering with isotropic view

---

## Material Characterization Workflow

```
1. LOAD EXPERIMENTAL DATA
        ↓
2. INITIAL GUESS (μ₀, β₀, K₀)
        ↓
3. PLOT INITIAL GUESS
        ↓
4. PARAMETER FITTING (lsqrsolve)
        ↓
5. PLOT FITTED CURVES
        ↓
6. COMPUTE TANGENT STIFFNESS (E_T)
        ↓
7. VISUALIZE STIFFNESS GLYPHS
        ↓
✓ CHARACTERIZATION COMPLETE
```

---

## Documentation

### Within This Repository
1. **`arruda_boyce_erklärung.md`** (German)
   - Complete theoretical background
   - Section-by-section code explanation
   - Derivations with intermediate steps
   - Nomenclature table (50+ symbols)

2. **`Herleitung_Tangensteifigkeit_Zwischenschritte.md`** (German)
   - Deep-dive tangent stiffness derivation
   - Step-by-step partial derivatives
   - Volumetric and deviatoric contributions

### External Resources
- **Scilab Documentation**: [https://www.scilab.org/](https://help.scilab.org/docs/2025.1.0/en_US/index.html)
- **Reiter, T. J. (2025). Lecture Notes:** Solid mechanics of continua
- **Holzapfel, G.A. (2000). Nonlinear Solid Mechanics:** A Continuum Approach for Engineering. New York: John Wiley & Sons. 
---

## Learning Outcomes

**Continuum Mechanics**
- Deformation gradients and strain tensors
- Stress measures (Cauchy, PK1, PK2)
- Hyperelastic constitutive models

**Computational Methods**
- Nonlinear optimization (Levenberg-Marquardt)
- Tensor calculus and 4th-order tensors
- Numerical differentiation

**Scientific Computing in Scilab**
- Advanced 3D visualization
- Matrix operations and linear algebra
- Modular code design

**Material Characterization**
- Experimental data fitting
- Parameter identification
- Model validation

**FEM Preparation**
- Constitutive model implementation
- Stiffness matrix computation
- Integration with implicit time-stepping

---

## Tips & Tricks

### Debugging Material Routine
```scilab
F = diag([2, 1/sqrt(2), 1/sqrt(2)]);
J_manual = det(F);
C = F' * F;
I1 = trace(C);
I1_bar = J_manual^(-2/3) * I1;
```

### Improving Fit Convergence
- Start with well-estimated initial values
- Normalize stresses to similar magnitudes
- Use combined fitting (all 3 modes)
- Inspect residual distribution (should be random)

### 3D Glyph Interpretation
- **Sphere at λ=1**: Isotropic undeformed material
- **Ellipsoid at λ>1**: Strain-induced anisotropy
- **Elongation factor**: Locking effect magnitude

---

## License

Educational resource for *"Solid Mechanics of Continua"*

**Usage:**
- ✅ Personal study and learning
- ✅ Academic coursework
- ✅ Non-commercial research

---

## Author

**Simon Boch**  
[GitHub](https://github.com/bochsimon-sketch)

**Last Updated:** April 7, 2026
