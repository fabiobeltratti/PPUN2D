# **PPUN2D** – Post-Processing Tool for Unstructured 2D CFD Solutions

**PPUN2D** (*Post-Processing Unstructured 2D Data*) is a C++ software developed for the aerodynamic force computation and drag decomposition from CFD solutions on **unstructured 2D grids**. It supports both **near-field** and **far-field** approaches and includes Green-Gauss gradient reconstruction and domain decomposition techniques. The software is intended for academic and research purposes.

> Developed at the *Italian Aerospace Research Centre (CIRA)*, in collaboration with the *University of Naples Federico II* as part of a master's thesis project.

---

## 🚀 Overview

PPUN2D processes binary **Tecplot (`.szplt`) files** obtained from CFD solvers and provides:

- Total aerodynamic force evaluation using the near-field approach.
- Support for **Lamb vector-based** and **thermodynamic** far-field methods.
- Drag breakdown into viscous, wave, and spurious components.
- Gradient and divergence computation via **Green-Gauss** techniques.
- Sensor-based **region tagging** (boundary layer, shock, wake).

Results are exported in both **Tecplot** and **plain-text** formats for visualization and further analysis.

---

## 📦 Input Requirements

- **Tecplot `.szplt` file** containing:
  - Primitive variables: `ρ`, `u`, `v`, `p`, `T`.
  - Optional turbulence fields: `μ_t`, `ω`.
  - Mesh connectivity and coordinates.

- **Configuration file** (`ppun2d_config.cas`) with:
  - Freestream conditions (Mach, Temperature, Density, etc.).
  - Flags for method activation.
  - Control volume specifications.
  - Sensor thresholds and surface definition parameters.

---

## 🧮 Implemented Methods

### Near-field method
Computes the aerodynamic force through surface integration of pressure and viscous stresses.
- It does not allow decomposition into physical and spurious drag contributions.
- Although often provided by the CFD solver, it is recalculated here to ensure consistent comparison with the implemented far-field methods.

### Far-field methods
Compute the aerodynamic force through integration of the net momentum flux on a control volume surrounding the body. These methods allow decomposition into physical and spurious contributions. The following formulations are implemented:
- **Paparone & Tognaccini method** → computes the irreversible drag as a function of the non-dimensional entropy.
- **Destarac & van der Vooren method** → computes the irreversible drag as a function of both entropy and total enthalpy.
- **Vorticity-based Lamb vector method**, where `l = ω × V` → computes both lift and drag.
- **Thermodynamic-based Lamb vector method**, where `l = T∇s − ∇h + ∇⋅τ/ρ` → computes lift and drag with reduced numerical sensitivity.

### Irreversible and parasite drag breakdown
PPUN2D supports detailed domain decomposition to isolate the sources of drag:
- **Viscous**, **shock-associated**, and **spurious** contributions can be separated.
- Decomposition is performed via:
  - *Internal sensors* based on eddy-viscosity ratio, specific dissipation rate, and pressure gradient.
  - *User-defined margins* extending each tagged region, set in the config file.
- Viscous regions are detected via high values of $(\mu + \mu_t)/\mu$ or $\omega$.
- Shock regions are flagged using pressure gradients.
- Spurious drag is computed as the residual from the total.
- All tagged regions are exported for direct visualization in Tecplot.

---

## 📐 Gradient and Divergence Schemes

- **Green-Gauss (GG) method**:
  - Both **cell-based** and **node-based** implementations are available.
  - Gradients and divergences are allocated at cell centers.
  - Used to compute $\nabla \phi$, $\nabla \cdot \boldsymbol{\tau}$, and other key quantities.

- Volume integrals are computed via cell-wise summation.


- Surface integrals use interpolated face values and outward normal vectors.

---

## 🗂 Output Files

- `output.szplt`: Tecplot file containing:
- Input and computed variables (Cl, Cd, Lamb vector, region tags, etc.).
- Compatible with Tecplot 360 for visualization.

- `forces_summary.txt`: structured table including:
- Lift and drag coefficients.
- Near-field and far-field results.
- Drag-breakdown by region (viscous, shock, spurious).
- Results for all evaluated wake-plane locations.

---

## 🔍 Applications

PPUN2D is suitable for:
- Post-processing of compressible 2D steady-state RANS solutions.
- Investigating drag contributions from different flow features.
- Benchmarking Lamb vector-based and thermodynamic drag breakdown methods.
- Supporting research in unstructured mesh CFD and aerodynamic optimization.

---

## 📖 References

Theoretical background and numerical formulations are based on:

Lamb vector-based methods:
- Minervino, M. & Tognaccini, R., A unified thermodynamic/Lamb-vector-based analysis of the aerodynamic force, Physics of Fluids, 2023. DOI: 10.1063/5.0164384
- Minervino, M. & Tognaccini, R., On the spurious effects in Lamb-vector-based force decomposition methods, Aerospace Science and Technology, 2023. DOI: 10.1016/j.ast.2023.108674
- Wu, J.-Z., Lu, X.-Y., & Zhuang, L.-X., Integral force acting on a body due to local flow structures, Journal of Fluid Mechanics, Vol. 576, 2007, pp. 265–286. DOI: 10.1017/S0022112006004551
- Mele, B. & Tognaccini, R., Aerodynamic force by Lamb vector integrals in compressible flow, Physics of Fluids, Vol. 26, No. 5, 2014, p. 056104. DOI: 10.1063/1.4875015

Thermodynamic formulations:
- Paparone, L. & Tognaccini, R., Computational Fluid Dynamics-Based Drag Prediction and Decomposition, AIAA Journal, Vol. 41, No. 9, 2003, pp. 1647–1657. DOI: 10.2514/2.7300
- Destarac, D. & van der Vooren, J., Drag/thrust analysis of jet-propelled transonic transport aircraft, Aerospace Science and Technology, 2004. DOI: 10.1016/j.ast.2004.03.004

Ph.D. Thesis:
- Minervino, M., A unified thermodynamic/vortical theory for the aerodynamic force analysis, Ph.D. thesis, University of Naples Federico II, 2025. hdl: 11588/989058

---

## 👨‍💻 Development and License

This code was developed by **Fabio Beltratti** as part of a master’s thesis in aerospace engineering.  
Supervisors: Ing. **Mauro Minervino** (CIRA), Prof. **Renato Tognaccini** (University of Naples "Federico II")  
License: Research use only – contact the author for usage permissions.

---

## 🙏 Acknowledgments

Special thanks to Ing. Mauro Minervino (CIRA) for his guidance and technical insights, and to Prof. Renato Tognaccini for scientific supervision and support throughout the development.

---

## 🛠 Build & Run

# Clone the repository
git clone https://github.com/your-username/PPUN2D.git
cd PPUN2D

# Compile with Intel compiler (icc)
icc -std=c++11 -o PPUN2D PPUN2D.cpp /path/to/tecio/libtecio.so

# Run the executable with the following syntax:
# ./ppun2d <input_solution.szplt> <config_file.cas> <output_file.szplt>
./ppun2d input_CFDsolution.szplt ppun2d_config.cas output_results.szplt
