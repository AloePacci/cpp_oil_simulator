## oilspillsim v0.1.0 (2026-11-17)

Python package `oilspillsim` 0.0.3 → 0.1.0 (`setup.py`). The core physics now behaves consistently when the time step
`_dt` changes, and several optional dynamics were added to make the spill evolve during an episode. The `reset`
signature of 0.0.3 (given sources, currents, wind and particles) is unchanged and works with every new option.

### BREAKING CHANGE

- **Physics: time-step-consistent defaults.** For the same seed and parameters, trajectories and density maps differ
  from 0.0.x (only the initial wind and source draws are equal). At the reference `_dt = 0.3` the *statistics* are the
  same as before (e.g. 803 vs 806 particles and the same spill area after ~160 steps); at other `_dt` they change on
  purpose. The old behaviour is still available, bit for bit (see *Legacy physics* below).
  - **Random walk ∝ √dt** (`diffusion_dt_ref = 0.3`). The Brownian term was `U(-2,2)·gamma` per step whatever `_dt`,
    so the spread after a given physical time shrank like 1/√dt. It is now multiplied by `sqrt(dt / diffusion_dt_ref)`:
    `gamma` keeps its meaning at `_dt = 0.3`, and the spread no longer depends on `_dt` (instantaneous release, 60 time
    units, `_dt` 0.15 / 0.3 / 0.6 / 1.2: before 6.3 / 4.6 / 3.3 / 2.3 cells, now 4.5 / 4.6 / 4.6 / 4.5).
  - **Emission per unit of time** (`emission_dt_ref = 0.3`). The sources released ~`_flow` particles per *step*
    whatever `_dt`. Now each step releases `_flow · dt / emission_dt_ref` particles, i.e. `_flow` is a rate per unit of
    time, equal to the old per-step rate at `_dt = 0.3` (60 time units: before 1994 / 997 / 499 / 250 particles for the
    four `_dt` above, now 1000 for all of them).
  - **Bug fix — skipped particle after a removal** (`fix_removal_skip = true`). In `step()`, when a particle left the
    grid it was erased inside the `for (i…; i++)` loop, so the particle that shifted into its column skipped its move
    that step. Now it is processed. On an open map with strong wind, 39.8 of 39.8 particles per step skipped their
    move before, 0 now. Maps whose edge is land are not affected (particles stop on land and never leave the grid).

### Feat

- **`set_dynamics(...)`** — one call to configure the dynamics; options persist across `reset()`. Calling it resets
  every option that is not passed to its default. All the *optional* dynamics are off by default and draw from their
  own random generator (seeded from the episode seed), so turning one on does not change the random sequence of the
  particles.
  - **Variable wind** (`wind_variability`, `wind_relaxation = 1/30`, `wind_max = √2`): Ornstein–Uhlenbeck process
    around the prevailing wind (the wind drawn at `reset`, readable as `wind_mean`), integrated with Euler–Maruyama.
    Stationary std per component = `wind_variability / sqrt(2·wind_relaxation)`, correlation time
    `1 / wind_relaxation`. |wind| is capped at `wind_max`, the envelope of the original `U(-1,1)²` draw, so the wind
    stays within the values the model always produced. Example: `wind_variability = 0.08` gives a std of ~0.27 per
    component and ~30° of direction change per 100 steps.
  - **Moving plume origin**:
    - `source_speed` (cells per time unit) and `source_turn` (heading random walk, rad/√time): a self-propelled source
      (e.g. a leaking vessel) that bounces back when its next position is land or outside the map.
    - `source_drift`: the source also drifts with the same forcing as the oil, `kw·wind + kc·current`
      (0 = none, 1 = like a particle).
    - The two can be combined. The source keeps a continuous position (`get_source_position()`) and `source_points`
      is its cell.
  - **Decay / evaporation** (`decay_rate`, `evaporable_fraction = 1`): each particle disappears with probability
    `1 − exp(−decay_rate·dt)` per step (exponential mass loss, half-life `ln 2 / decay_rate` time units).
    `evaporable_fraction < 1` marks only that fraction of the new particles as volatile; the rest is a persistent
    residue. Checked: 2001 → 1012 particles after one half-life.
  - The three time-step corrections above are also options of `set_dynamics`: `diffusion_dt_ref`,
    `emission_dt_ref` (≤ 0 = legacy) and `fix_removal_skip`.
- **`warmup(n_steps)`** — runs `n_steps` steps in C++ with the GIL released, to start the spill before an episode.
  Identical to calling `step()` `n_steps` times, and ~3× faster than stepping from Python through an environment
  wrapper.
- **`set_wind(wind, also_mean=True)`** — set the current wind; with `also_mean` it also becomes the prevailing wind
  of the variable-wind process.
- **`set_source_fuel(fuel)`** — finite release: emission stops when the fuel runs out (with a large `_flow`, an
  instantaneous puff).
- **New read-only attributes and properties**: `wind_mean`, `step_count`, `remaining_fuel`, `particle_count`, and the
  method `get_source_position()`.

### Fix

- **Seed reproducibility independent of the compiler.** `v << rand(gen), rand(gen)` (Eigen comma initialiser) has an
  evaluation order that depends on the compiler and C++ standard: the published 0.0.1 wheel (GCC 10) put the first draw
  in the second component, a GCC 11 build did the opposite, so the same seed gave a different spill. The draws are
  now explicit, in the order of the published wheel (wind at `reset` and both random-walk terms in `step`).
- **CMake build of the Python module.** `PYBIND11_MODULE(oilspillsim, …)` was built by the CMake target `SIMULATOR`,
  so the resulting `SIMULATOR*.so` could not be imported (`PyInit_SIMULATOR` missing). The module name is now the
  macro `OSS_MODULE_NAME`: `setup.py` defines it as `oilspillsim` (pip package, unchanged API) and the CMake build
  gets a module called `SIMULATOR`.
- With `_random_seed = -1` the generator of the optional dynamics is also seeded randomly.

### Legacy physics

To reproduce data generated with 0.0.x exactly:

```python
sim.reset(_seed=seed)
sim.set_dynamics(diffusion_dt_ref=-1, emission_dt_ref=-1, fix_removal_skip=False)
```

With these options the build is **bit-identical** to the published `oilspillsim` 0.0.1 wheel (density maps, particle
positions and wind; checked on 3 seeds × 2 configurations, including `_dt` 0.6, `_flow` 15 and 3 sources).

### Tests

- `tests/test.py` (run by CI) now has 10 assert-based tests instead of only importing the module (~1 s):
  - basics and seed reproducibility, including `warmup(n)` = n × `step()`;
  - the optional dynamics leave the particle random sequence untouched;
  - emission per unit of time and √dt diffusion, each checked against legacy mode;
  - no skipped particle after a removal (open map);
  - variable wind (cap and variability);
  - moving source (never on land);
  - decay half-life and residual fraction;
  - finite release, and `reset` with a given source, wind and particles combined with the dynamics.

## v2.0.0 (2024-12-19)

### BREAKING CHANGE

- VAE paper version out

### Feat

- **Simulation**: released paper version

## v1.0.0 (2024-09-26)
