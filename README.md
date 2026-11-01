# Solar Sprint Vehicle Simulations

Monte Carlo simulation in Fortran, parallelized with OpenMP, for a Pitsco Solar Sprint car. It models motor torque, drag and friction across thousands of car configurations. The simulation set the car's gearing and predicted the real car's performance within 10%. The car was 4th fastest of 16 in Senior Solar Sprint at the 2025 TSA National Leadership Conference.

Written by Vivaan Rochlani.

## Programs

- `gear_optimizer.f90`: sweeps gear ratios from 1.0 to 8.0, running 10,000 randomized simulations per ratio (OpenMP).
- `speed_project_2.f90`: 100,000-iteration Monte Carlo of the 20 m run (OpenMP).
- `friction_analysis.f90`: deterministic force breakdown over the 20 m course.
- `sensitivity.f90`: one-at-a-time sensitivity of finish time to nine parameters.
- `terminal_finder.f90`: terminal velocity detector over an extended distance.
- `speed_project.f90`: early prototype, kept for reference.
- `elastic_inelastic.f90`: two-cart momentum calculator.

## Build

```bash
gfortran -O2 -fopenmp -o gear_optimizer gear_optimizer.f90
./gear_optimizer
```

Requires `gfortran`; OpenMP (`-fopenmp`) for the Monte Carlo programs.

## Physical model

All programs share the same core physics unless noted.

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
| Time step | 1 ms (0.1 ms in `terminal_finder.f90`) |
