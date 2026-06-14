# Solar Sprint Vehicle Simulations

Fortran simulations written for the Pitsco Solar Sprint competition. The vehicle is a ~300 g motor-driven car running a 20 m straight course, powered by a Pitsco 280 DC motor on two AA batteries (3 V nominal). Each program targets a different analysis question — from single-run force accounting to multi-thousand-iteration Monte Carlo sweeps.

## Requirements

- `gfortran` (GCC 10 or later recommended)
- OpenMP (`-fopenmp`) for `gear_optimizer.f90`, `speed_project.f90`, and `speed_project_2.f90`

On macOS with Homebrew:
```bash
brew install gcc
```

---

## Files

### `elastic_inelastic.f90`
**Two-cart momentum calculator**

Interactive program for conservation-of-momentum labs. Prompts for the mass and initial/final velocity of two carts, then prints individual momenta before and after collision and the ratio p_final / p_initial for the system.

```bash
gfortran -o momentum elastic_inelastic.f90
./momentum
```

---

### `friction_analysis.f90`
**Deterministic force simulation — 20 m run**

Steps through the 20 m course at 1 ms intervals using mean parameter values (no randomization). At each step it computes the motor torque from a linear torque–speed curve, applies a power cap, then evaluates bearing friction, rolling resistance, aerodynamic drag (rectangular box body, Cd = 1.05), and expected guide wire friction. Prints a force table every 100 steps and a full breakdown of average forces and drag share at the end.

Useful for a quick sanity check on the physics or for tracing exactly how each force evolves over the run.

```bash
gfortran -O2 -o friction friction_analysis.f90
./friction
```

---

### `gear_optimizer.f90`
**Monte Carlo gear ratio sweep — OpenMP parallelized**

Sweeps gear ratios from 1.0 to 8.0 in steps of 0.125 (57 ratios total). For each ratio it runs 10,000 simulations with randomized bearing friction, rolling resistance, motor no-load speed, stall torque, and power output. Guide wire contact is modeled as a 30% per-step probability. Reports average finish time, average exit speed, and valid-run percentage for every ratio, then identifies the best.

Runs across 100 OpenMP threads by default; adjust `num_threads` at the top of the file to match your hardware.

```bash
gfortran -O2 -fopenmp -o gear_optimizer gear_optimizer.f90
./gear_optimizer
```

---

### `sensitivity.f90`
**One-at-a-time local sensitivity analysis**

Perturbs each of 9 physical parameters by ±10% from the baseline, runs the deterministic simulation for each perturbation, and computes normalized sensitivity (% change in finish time per 10% parameter change). Results are sorted by absolute impact and annotated with the direction of improvement (increase or reduce).

Parameters analyzed: vehicle mass, drag coefficient, frontal area, motor stall torque, motor no-load RPM, gear efficiency, rolling friction coefficient, bearing friction coefficient, and guide wire friction.

```bash
gfortran -O2 -o sensitivity sensitivity.f90
./sensitivity
```

---

### `speed_project.f90`
**Early Monte Carlo prototype — 4000 m course**

Original simulation written before the Solar Sprint configuration was finalized. Models a heavier vehicle (1.2 kg) over a 4000 m course with a quadratic torque curve. Does not include aerodynamic drag or guide wire friction. Kept for reference; the physics assumptions here diverge from the final vehicle.

```bash
gfortran -O2 -fopenmp -o speed_project speed_project.f90
./speed_project
```

---

### `speed_project_2.f90`
**Production Monte Carlo — 20 m sprint, 100k iterations**

Full-fidelity simulation of the competition run. Randomizes the same six parameters as `gear_optimizer.f90` and additionally tracks per-step guide wire contact stochastically. Runs 100,000 iterations across 100 OpenMP threads and reports average, minimum, and maximum finish time and exit speed.

This is the primary tool for estimating expected race performance and quantifying run-to-run variance.

```bash
gfortran -O2 -fopenmp -o speed_project_2 speed_project_2.f90
./speed_project_2
```

---

### `terminal_finder.f90`
**Terminal velocity detector — extended distance run**

Runs the deterministic simulation over up to 500 m (finer 0.1 ms time step) and watches for acceleration dropping below 0.005 m/s² for 100 consecutive steps. When that condition is met it records and reports the position, elapsed time, and speed at which the vehicle reached terminal velocity. If the flatline is never detected within the distance limit the program reports the final speed instead.

```bash
gfortran -O2 -o terminal_finder terminal_finder.f90
./terminal_finder
```

---

## Physical model summary

All programs share the same core physics unless noted:

| Parameter | Value |
|---|---|
| Motor | Pitsco 280, linear torque–speed curve |
| Stall torque | 12.7–13.0 mN·m |
| No-load speed | 5000–5400 RPM at 3 V |
| Supply | 3 V (2× AA), 1.55 A at max efficiency |
| Motor efficiency | 60% |
| Gear efficiency | 98% |
| Wheel diameter | 39.2 mm |
| Bearings | 4× ZrO2 full ceramic, μ = 0.0005–0.002 |
| Rolling resistance | μ = 0.004–0.008 (latex on paper) |
| Aerodynamic drag | Cd = 1.05, A = 25 mm × 100 mm |
| Guide wire friction | μ = 0.05 on nylon, 30% contact probability |
| Vehicle mass | 297 g (Monte Carlo baseline) |
| Course distance | 20 m |
| Time step | 1 ms (0.1 ms in terminal_finder) |
