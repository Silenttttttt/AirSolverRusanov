# AirSolverRusanov — gas (air) + heat, split for TPT

**Two files:**
- **Gas (air)** — `air_solver_rusanov.cpp`: convection (Rusanov) + viscosity. Used by TPT’s air.
- **Heat** — `heat.hpp` / `heat.cpp`: thermal diffusion k∇²T, heat capacity c_v/c_p, temperature. Used by particles (later) and optionally by air.

**Stable by design:** first-order Rusanov + viscosity; heat diffusion in separate module. Conservative: mass and energy conserved (periodic BCs).

**Physics:** density ρ, pressure p, velocity (u,v), temperature T, viscosity μ, thermal conductivity k, heat capacities c_v and c_p (ideal gas, γ=1.4). Gas: convection + viscous stress τ. Heat: k∇²T.

## Why this exists

FluidSim (MacCormack) is hard to keep stable and performant. This solver trades **accuracy** for **stability and simplicity**:

- **First-order** → no reconstruction, so no negative density/energy from overshoots.
- **Rusanov flux** → very dissipative but stable; CFL ≤ 0.5 in 1D, we use 0.35 in 2D.
- **Viscosity & heat** → full Navier–Stokes: μ=1.8e-5 Pa·s, k from Pr=0.71; dt limited by CFL and diffusion.
- **Conservative** → same physics (compressible ideal gas); mass and energy conserved.

## Build and run

**Headless (conservation test):**
```bash
g++ -O2 -std=c++17 -o air_rusanov air_solver_rusanov.cpp heat.cpp
./air_rusanov [ny] [nx]   # default 40 60
```
Output: step count, sim time, relative mass and energy error (should be ~1e-14).

**With SDL2 visualization (pressure colormap, click to add bumps, Space = pause):**
```bash
g++ -O2 -std=c++17 -DUSE_SDL -o air_rusanov air_solver_rusanov.cpp heat.cpp $(pkg-config --cflags --libs sdl2)
./air_rusanov
```
Requires SDL2 dev package (e.g. `sudo pacman -S sdl2` on Manjaro). **Left-click** = add pressure bump; **right-click or drag** = add wall cells; **middle-click or drag** = remove wall cells; **Space** = pause; window title and terminal show mass/energy err.

## Use in your project (e.g. Powder Toy)

- **No external deps** — only C++17 and STL.
- **Gas (air):** `air_solver_rusanov.cpp` — `State(ny, nx, dx)`; `step()` (convection + viscosity); `apply_heat_diffusion(dt)` (calls heat module). Constants: `R_gas`, `c_v`, `c_p`, `mu`.
- **Heat (particles later):** `heat.hpp` / `heat.cpp` — `heat::heat_diffusion_step(ny, nx, dx, dt, rho, rhou, rhov, E)`; `heat::temperature_from_e(e)`; `heat::max_heat_dt(...)`; constants: `heat::R_gas`, `heat::c_v`, `heat::c_p`, `heat::k_thermal`, `heat::Pr`.
- **API (gas):** `set_uniform(rho, ux, uy, p)`; `add_bump(cy, cx, dP, rad)`; `double dt = step()` (returns dt used, or 0 if rejected); `total_mass()`, `total_energy()`; `primitives(iy, ix, rho, ux, uy, e, p)`; `temperature(iy, ix)` (Kelvin); read/write `U0[l][iy][ix]` for conserved [ρ, ρu, ρv, E]. **Walls:** `set_wall(iy, ix, true)` for solid; `is_wall(iy, ix)`; `set_boundary_walls()` sets the domain edge cells as walls (reflective box, no periodic wrap). Default all false (periodic). Reflective slip + adiabatic heat at walls. The SDL demo calls `set_boundary_walls()` so the run is a box, not wormholes.

**Porting to The Powder Toy (TPT)** — Moderate effort (roughly half a day to a day if you know both codebases):

| Task | Effort | What |
|------|--------|------|
| **Data mapping** | Easy | PT has `pv`, `hv`, `rho`, `vx`, `vy`. Convert to conserved U: ρ from `rho`, u/v from `vx`/`vy`, E from ρ, p (from `pv`), and e = p/((γ−1)ρ). Copy back: from U get ρ, u, v, p; set `pv`=p, `vx`=u, `vy`=v, `hv`=p/(ρ·R) if PT needs temperature. Use double in solver, copy to/from PT’s `float` grids. |
| **Grid** | Easy | One `State(ny, nx, dx)` with ny=YCELLS, nx=XCELLS; same shape as PT’s air grid. |
| **Walls** | Medium | Solver is currently periodic. Add a wall mask (e.g. from `bmap_blockair`): at cell faces next to a wall, use zero normal flux (no mass/energy through wall; momentum flux = pressure from fluid side). Either pass a `is_wall(iy,ix)` into `step()` and branch in the flux loop, or precompute ghost cells. |
| **Call site** | Easy | In `Air::update_air()`: create or reuse a gas `State`, copy in `pv`/`vx`/`vy`/`rho` → U; call `step()` then `apply_heat_diffusion(dt)` (or skip heat for air-only); copy U → primitives → `pv`/`vx`/`vy`/`rho`; set `hv` from `temperature(iy,ix)` (Kelvin). Link `heat.cpp` if you use heat. |

So: **not hard** — one file, no deps, same grid and physics. The only non-trivial part is implementing wall BCs in the flux logic (stop using periodic wrap at walls; use zero flux or mirror state).

## Will it work? Trade-offs

**Yes.** First-order + Rusanov is a standard, provably stable choice for compressible Euler. Under CFL and with step rejection you get:

- **Stability** — no overshoots (no reconstruction), no blow-ups; invalid steps are rejected and state is rolled back.
- **Conservation** — mass and energy are conserved to machine precision with periodic BCs (your headless run shows ~1e-14 error).
- **Real physics** — same equations (compressible Euler, ideal gas); no artificial clamps or floors on the conserved state.

**What we trade off:**

| We give up | We get |
|------------|--------|
| **Sharp shocks** | First-order + Rusanov is dissipative: shocks and contacts smear over several cells. You see smooth pressure waves, not crisp discontinuities. |
| **High-order accuracy** | No reconstruction → no second-order in space. Good enough for “air feels right” in a game; not for high-fidelity CFD. |
| **Fancy Riemann solver** | Rusanov is simple and robust; HLLC/Roe would be sharper but more fragile and more code. |

So: it will run stably and conserve; you trade **sharpness/accuracy** for **stability and simplicity**.

## Physics constants (in code)

- **Gas (air_solver_rusanov.cpp):** γ = 1.4, R = 287 J/(kg·K), c_v, c_p, μ = 1.8e-5 Pa·s. ρ_min, e_min, CFL = 0.35; dt limited by CFL and viscous diffusion.
- **Heat (heat.hpp/cpp):** R, γ, c_v, c_p, Pr = 0.71, k_thermal = μ·c_p/Pr. `max_heat_dt()` for thermal diffusion limit.

## Tuning

- **CFL** — in code, `const double CFL = 0.35`. Lower = more stable, smaller dt. Raise to 0.4–0.45 if you never reject.
- **Step rejection** — if `step()` returns 0, the update would have produced negative ρ or e; state is unchanged. Retry with smaller dt (e.g. halve) or keep current state for one frame.
- **Viscosity / heat** — change `mu` or `k_thermal` (or `Pr`) to tune dissipation and thermal diffusion.

## How to improve (while keeping it quick)

- **Resolution** — In the SDL `main`, increase `ny`/`nx` (e.g. 120×180) and/or reduce `scale` for more cells and sharper waves; cost scales with cell count.
- **CFL** — Bump to 0.4 if you see no step rejects; you advance more per step and stay fast.
- **Second-order (MUSCL + minmod)** — Add limited linear reconstruction (minmod slope) and evaluate flux from face-left/face-right states. Sharper shocks, still stable; ~20–30% more cost per step. Keep Rusanov flux and step rejection.
- **Better colormap** — e.g. show density or speed instead of pressure, or overlay velocity arrows every N cells for flow direction.
- **Walls / boundaries** — For PT: zero normal flux at solid faces (copy cell state into ghost, flip normal velocity), so you can do rooms and obstacles without changing the core solver.
- **Multiple fluids** — Same solver; add a passive scalar (e.g. smoke concentration) advected with the flow if you need tracers.
