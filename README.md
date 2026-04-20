# Molecular Dynamics Simulation of Lennard-Jones Particles

A C++ header-only molecular dynamics (MD) engine for simulating a 3D Lennard-Jones
fluid. Developed as a university project for the *Computational Statistical Mechanics*
(CSM) course.

Author: **Luciano Lanzotti**

> **Note:** The base implementations of `sim.hpp`, `randnumgen.hpp`, `particle.hpp`, `params.hpp`,  and `pvector.hpp`
> were provided by **Prof. Cristiano De Michele** — University of Rome "La Sapienza" — as course material
> for the *Computational Statistical Mechanics* course.

---

## Features

| Feature | Description |
|---|---|
| **Lennard-Jones potential** | Shifted and unshifted, with configurable cutoff radius |
| **Nosé-Hoover thermostat** | NVT ensemble via extended-system equations of motion |
| **Brownian dynamics** | Langevin integrator (stochastic, NVT) |
| **NVE integrator** | Plain velocity-Verlet for microcanonical runs |
| **Periodic boundary conditions** | Minimum-image convention in all three directions |
| **Radial distribution function** | g(r) accumulated over production steps |
| **Static structure factor** | S(k) via random wave-vector sampling |

---

## File Overview

```
.
├── pvector.hpp      # Generic N-dimensional fixed-size vector class
├── randnumgen.hpp   # Mersenne Twister (mt19937_64) random number generator
├── params.hpp       # Simulation parameters container (simpars)
├── particle.hpp     # Base particle class + LJ particle (forces, potential)
├── nh.hpp           # Nosé-Hoover thermostat propagators
├── brownian.hpp     # Brownian dynamics (Langevin) propagator
├── radial.hpp       # Radial distribution function g(r)
├── static.hpp       # Static structure factor S(k)
├── sim.hpp          # Core simulation engine (sim + mdsim)
├── test_dt.cpp      # Time-step convergence test (NH integrator)
├── test_U.cpp       # Energy conservation test (NH and Brownian)
├── test_g.cpp       # g(r) at 4 densities, NH and Brownian
└── test_ssf.cpp     # S(k) at 4 higher densities, NH and Brownian
```

---

## Requirements

- C++17-compatible compiler (GCC ≥ 7, Clang ≥ 5, MSVC ≥ 19.14)
- Standard library only — no external dependencies

---

## Building

Each test program is a single translation unit that includes all headers. Compile with:

```bash
# Time-step convergence test
g++ -O2 -std=c++17 test_dt.cpp -o test_dt

# Energy conservation test
g++ -O2 -std=c++17 test_U.cpp -o test_U

# Radial distribution function
g++ -O2 -std=c++17 test_g.cpp -o test_g

# Static structure factor
g++ -O2 -std=c++17 test_ssf.cpp -o test_ssf
```

---

## Running

Before running, create the required output directories:

```bash
mkdir -p dati/dt dati/u dati/rdf dati/ssf
```

Then execute any compiled binary:

```bash
./test_dt
./test_U
./test_g
./test_ssf
```

Each program prints progress to stdout and writes results to tab-separated `.dat` files
inside `dati/`. These files can be read directly with NumPy (`np.loadtxt`) or similar.

---

## Default Simulation Parameters

| Parameter | Value | Description |
|---|---|---|
| `nx × ny × nz` | 5 × 5 × 5 = 125 | Number of particles |
| `sigma` | 1 | LJ length scale |
| `epsilon` | 1 | LJ energy scale |
| `rc` | 2.5 σ | Cutoff radius |
| `rho_star` | 0.4 | Reduced density ρ\* = ρ σ³ |
| `T_eq` | 1.4 | Target temperature (LJ units) |
| `dt` | 0.01 τ | Integration time step |
| `t_tot` | 3 τ | Total simulation time |

All parameters can be changed programmatically via the `simpars` setters before calling
`prepare_initial_conf()` and `run()`.

---

## Output Format

The main output file (default `dati/results.dat`) has a comment header line followed by
one row per saved step:

```
#dt:<dt> totstep:<N> tau:<tau> t_tot:<t_tot> rho:<rho> NH:<0/1> B:<0/1>
<step>   <time>   <K+U_shifted>   <E_NH_conserved>   <T_instantaneous>
```

The `g(r)` and `S(k)` files contain two columns: the observable argument (r or k) and
the corresponding function value.

---

## Physical Model

The particles interact via the truncated (and optionally shifted) Lennard-Jones potential:

```
V(r) = 4 ε [ (σ/r)^12 − (σ/r)^6 ]   for r < rc
V(r) = 0                               for r ≥ rc
```

The shifted variant subtracts V(rc) to make the potential continuous at the cutoff.

---

## License

This code was written for educational purposes as part of a university course assignment.
