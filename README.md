# Underexpanded Supersonic Jet: Godunov Marching Method

Computation of a planar underexpanded supersonic jet (p_j/p_a > 1) exhausting from a nozzle into a quiescent medium. The Godunov marching method is used: the coordinate x along the jet plays the role of time, and the fluxes through cell faces are obtained from the exact solution of the steady discontinuity breakup problem (oblique shocks and Prandtl–Meyer waves). The boundary conditions are also imposed through breakups: on the symmetry axis and on the free jet boundary, where p = p_a.

The code is written in GNU Octave. It will not run in MATLAB without changes: it uses Octave-specific syntax, for example `figure(1, 'Position', ...)`.

## Parameters

The parameters are set at the beginning of [src/main.m](src/main.m):

| Parameter | Default | Meaning                                           |
| --------- | ------- | ------------------------------------------------- |
| `gamma`   | 1.4     | ratio of specific heats                           |
| `M`       | 2.0     | jet Mach number at the nozzle exit                |
| `npr`     | 2       | pressure ratio p_j/p_a                            |
| `N`       | 100     | number of cells across the half-jet               |
| `kappa`   | 0.8     | safety factor in the step condition for dx        |
| `every`   | 1       | print every `every`-th step (0 — do not print)    |
| `w_noz`   | 0.01 m  | nozzle half-width at the exit                     |
| `x_len`   | 0.2 m   | length of the computational domain along the axis |

Internally the computation uses dimensionless quantities: lengths are scaled by `w_noz`, and p_a = 1, ρ_j = 1.

## Output

During the run, the step number, x, dx and the jet boundary position y_e are printed to the console. At the end, the Mach number field in the upper half-jet is plotted together with its boundary.

The computation stops with an error in three cases:

- u ≤ c somewhere in the flow, and the marching method becomes inapplicable;
- a shock detaches or a vacuum forms in a breakup;
- the discriminant becomes negative when recovering the parameters from the fluxes (choking). This usually means the step is too large; reducing `kappa` helps.

## Structure

```
docs/physics.md     description of the method: equations, grid, scheme, discontinuity
                    breakup, boundary conditions, step choice, algorithm, references
docs/physics-en.md  English version of the same description
src/main.m          entry point: parameters and the marching loop
src/riemann_solve.m, wave_f.m, shock_rmax.m, ...   discontinuity breakup
src/riemann_flux.m, face_flux.m                    fluxes through faces
src/axis_boundary.m, free_boundary.m               boundary conditions
src/recover.m, flux_x.m                            conversion between F and (p, ρ, u, v)
src/plot_results.m                                 plotting
```

The method is described in detail in [docs/physics-en.md](docs/physics-en.md) (Russian original: [docs/physics.md](docs/physics.md)); the comments in the code refer to its sections.

## Limitations

- The flow must remain supersonic everywhere (u > c). Mach reflection, subsonic regions and reverse flows cannot be computed with the marching method.
- The model is inviscid: mixing of the jet with the ambient medium and boundary layers are not described.
- The gas is perfect, γ = const.
