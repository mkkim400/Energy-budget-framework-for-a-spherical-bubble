# Bubble Dynamics Simulation and Energy Budget Analysis

This repository contains MATLAB scripts to simulate the radial dynamics of a gas bubble in a liquid medium using the **Keller-Miksis (KM)** equation. The model considers effects such as compressibility, viscosity, and surface tension. These simulations are useful for research in cavitation, sonoluminescence, and acoustic bubble dynamics.

---

## Contents

- `run_RPsolver.m`: Main driver script to run multiple KM simulations over a range of pressure ratios.
- `KM_solver.m`: Function that implements the KM model and numerically integrates it over time using `ode15s`.

---

## Overview

### 1. `run_RPsolver.m`

This script performs a parametric study by varying the **pressure ratio (PR)** and solving the Keller-Miksis equation for each case.

#### Key Features:
- Supports multiple physical models:
  - With/without viscosity
  - With/without surface tension
- Supports bubble **growth** or **collapse**
- Automatically saves:
  - Raw simulation data (in physical units and non-dimensional form)
  - Bubble radius evolution over time
- Computes and records:
  - Maximum bubble pressure
  - Minimum radius
  - Corresponding gas temperature inside the bubble

---

### 2. `KM_solver.m`

This function solves the **non-dimensional** form of the Keller-Miksis equation. It returns physical output variables by scaling the results after solving the ODE.

#### Output:
- `t`: Time vector [s]
- `y1`: Bubble radius [m]
- `y2`: Wall velocity [m/s]
- `y3`: Bubble pressure [Pa]
- `y4`: Wall acceleration [m/s²]
- `y5`: Pressure derivative [Pa/s]

#### Input Parameters:
- `PR`: Pressure ratio
- `equilb`: Flag for equilibrium (initial wall velocity)
- `sig`: Surface tension [N/m]
- `mu_l`: Liquid viscosity [Pa·s]
- `c_l`: Speed of sound in liquid [m/s]
- `gam_v`: Specific heat ratio of gas
- `p_init`: Initial pressure inside bubble [Pa]
- `p_v`: Vapor pressure [Pa]
- `R0`: Initial bubble radius [m]
- `grow`: Growth (0) or collapse (1)
- `rho_l`: Liquid density [kg/m³]
- `t0`: Initial time [s]

---

## How to Use

1. Open `run_RPsolver.m` in MATLAB.
2. Set desired parameter ranges (e.g., `PRi`, `R0`, model flags).
3. Run the script to generate results.
4. Output files will be saved in automatically generated directories based on simulation settings.

> **Note:** Ensure `KM_solver.m` is in the MATLAB path or in the same directory as `run_RPsolver.m`.

---

## Output Files

Each simulation case generates:
- `R_Q.txt`: Bubble radius, velocity, acceleration, pressure, and pressure rate (in physical units)
- `R_Q_non.txt`: Same variables, but non-dimensionalized
- Directories are auto-labeled with physical parameters for traceability

---

## References

- Keller, J.B., and Miksis, M. (1980). *Bubble oscillations of large amplitude*. Journal of the Acoustical Society of America.
- Prosperetti, A., and Lezzi, A. (1986). *Bubble dynamics in a compressible liquid. Part 1. First-order theory.* Journal of Fluid Mechanics.

---

## License

This project is for academic and research purposes. Contact the author for redistribution or publication rights.

