# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Computation of a planar underexpanded supersonic jet (p_j/p_a > 1) by the Godunov marching method with the exact solution of the steady discontinuity breakup problem. The code is GNU Octave (Octave-specific syntax is used, e.g. `figure(1, 'Position', ...)`); there is no build, no tests and no linter.

## Running

The only entry point is the script `src/main.m` (it starts with `clear`; the problem parameters are at its beginning: `gamma`, `M`, `npr`, `N`, `kappa`, `every`, geometry `w_noz`, `x_len`):

```sh
cd src && octave --persist main.m      # --persist keeps the plot window open
```

or from an Octave session: `cd src; main`. All functions live in `src/` and are on the path only when run from that folder (or after `addpath('src')`).

## Documentation

`docs/physics.md` is a detailed description of the method (equations, grid, scheme, discontinuity breakup, boundary conditions, step choice, algorithm). The English version is `docs/physics-en.md`; section numbering is the same. Code comments refer to its sections ("section 5.6", "section 6.3", etc.) — when changing the algorithm, check against it and keep the references consistent.

## Architecture

The marching coordinate x plays the role of time. The state on a layer is N cells across the upper half-jet between the axis (y = 0) and the free boundary y_e(x); nodes are distributed uniformly, `y = eta*ye`. The conservative variables are the flux vector through the cross-section x = const, `F = [ρu, p+ρu², ρuv, ρuH]`, not (ρ, ρu, …).

The loop in `main.m` at each step (corresponds to section 8 of physics.md):
1. Marching check `u > c` in every cell (error otherwise).
2. Exact breakups on interior faces: `riemann_solve(WL, WR)` → states `WLS`, `WRS` on both sides of the slip line.
3. Boundary conditions also via breakups: `axis_boundary` (axis as a mirror image of cell 1, θ* = 0 → pressure p_w) and `free_boundary` (C− wave to p = p_a → boundary angle θ_b).
4. Step dx — a priori estimate from Mach line slopes with the factor `kappa` (section 7).
5. New layer: boundary shift `ye + dx*tan(θ_b)`, face slopes `s`.
6. Fluxes through faces: `riemann_flux` selects the zone of the self-similar solution (L | C− | L* | slip line | R* | C+ | R, including the interior of a fan) by the ray ξ = s and computes the flux via `face_flux`.
7. Scheme `Fn = (h·F − dx·ΔΦ)/h_new`, then `recover` inverts F → (p, ρ, u, v) using the larger (supersonic) root of a quadratic equation; `Dmin < 0` is an error (choking).

Conventions shared by the functions:
- Inside breakups the state is `W = [p, rho, q, theta]` (speed magnitude and angle), in cells — `p, rho, u, v`; conversion via `flow_speed`/`flow_angle`.
- Functions are vectorized: `[n x 4]` arrays, one row per face/cell.
- `wave_f` is the turning angle f_k(p) and its derivative: oblique shock for p ≥ p_k, Prandtl–Meyer wave for p < p_k. `riemann_solve` (Newton safeguarded by bisection), `axis_boundary` and `free_boundary` all rely on it.
- Shock detachment is checked via `shock_rmax` (limiting p2/p1 at δ_max); shock detachment and vacuum errors are raised with `error`.
- Comments, error messages, console output and plot labels are in English (terminology as in `docs/physics-en.md`); each function header describes inputs/outputs in a uniform format (`Input:` / `Output:`) — keep it in new functions.
