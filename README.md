# cpp_py_oil_simulator

A C++ oil spill simulator with Python bindings (pybind11). This repository provides a performant core simulator written in C++ and a small pybind11 wrapper so the simulator can be used from Python as the `oilspillsim` package.

Highlights
- Lagrangian particle model: particles are advected by wind and currents, follow a random walk, are emitted by sources and are rendered as a normalized density map.
- Physics consistent with the time step `_dt` (since 0.1.0): diffusion ∝ √dt and emission per unit of time.
- Optional dynamics (since 0.1.0): time-varying wind, moving plume origin, decay / evaporation, finite releases, and a fast `warmup`.
- Core numerical code uses Eigen for linear algebra (bundled under `src/Eigen`).
- Python bindings implemented in `src/utils/pywrap/pywrap.cpp` using pybind11.
- Minimal example and utilities live under `src/` and `utils/`.

See [CHANGELOG.md](CHANGELOG.md) for the full list of changes.

## Quick install

Install using pip from the project root:

```bash
pip install .
```

This builds the extension module (`oilspillsim`) with pybind11. Ensure you have a C++ toolchain and CMake available on your platform. Inside an existing environment that already has `pybind11` and `setuptools`, `pip install --no-build-isolation .` avoids downloading them again.

## Build (developer)

Build the native extension directly using the provided CMake and setuptools integration:

```bash
# from project root
pip install -e .
```

Or build with CMake in `src/` (native build):

```bash
mkdir -p build && cd build
cmake ../src
cmake --build .
```

The CMake build produces the C++ executable `BasicSimulator` and a Python module named `SIMULATOR` (`import SIMULATOR; SIMULATOR.SIMULATOR(...)`), with the same API as the pip package. The pip build names the module and the class `oilspillsim`. The name is set by the `OSS_MODULE_NAME` macro in `pywrap.cpp`.

## Usage (Python)

After installation the package exposes the `oilspillsim` extension module:

```py
import numpy as np
from oilspillsim import oilspillsim

base = np.ones((100, 100), np.int32)          # 1 = water, 0 = land / obstacle
base[:8, :] = base[-8:, :] = base[:, :8] = base[:, -8:] = 0

sim = oilspillsim(base, _dt=0.3, _kw=1.0, _kc=1.0, _gamma=0.3, _flow=5,
                  _number_of_sources=1, _max_contamination_value=5,
                  _source_fuel=100000, _random_seed=42)
sim.reset(_seed=42)

# optional dynamics (all off by default)
sim.set_dynamics(wind_variability=0.08,   # time-varying wind
                 source_speed=0.3,        # moving source (cells / time unit)
                 decay_rate=0.01)         # evaporation (1 / time unit)

sim.warmup(200)                           # start the spill: 200 steps in C++
for t in range(300):
    sim.step()
    density = sim.get_normalized_density()   # (rows, cols) in [0, 1]
    wind = sim.wind_speed                     # changes every step with wind_variability > 0
```

A map can also be loaded from a file: `oilspillsim("src/mapas/defuniak_map.csv", ...)`.

## Accessible functions (API)

The `oilspillsim` extension exposes a class named `oilspillsim` (C++ type `SIMULATOR`). Below is a concise summary of constructors, methods and read-only attributes provided by the pybind11 wrapper. For exact function signatures and default values see `src/utils/pywrap/pywrap.cpp`.

- Class: `oilspillsim` (C++: `SIMULATOR`)
  - Constructor:
    - `oilspillsim(_filepath: str, **kwargs)`
      - Create the simulator by loading a map from a file (CSV or supported format by the core).
    - `oilspillsim(base_matrix: numpy.ndarray or Eigen::MatrixXi, **kwargs)`
      - Create the simulator from an in-memory base matrix (map, usually a 2D array, 0 for obstacle 1 for valid position).
  - Contructor kwargs:
    - `_dt`: (float=10.0) internal time variable. Each time oilspillsim.step() is called, particles will move relative to this.
    - `_kw`: (float=0.5) wind gain. wind component affecting particles is weighted by this value
    - `_kc`: (float=1.0) current gain. current components affecting particles are weighted by this value.
    - `_gamma`: (float=1.0) brownian movement gain. Random movement due to particle intereffects are weighet by this value. Since 0.1.0 the step is scaled by `sqrt(_dt / diffusion_dt_ref)` (see *Physics*).
    - `_flow`: (float=50.0) emission rate of each source. Since 0.1.0 it is a rate per unit of time: each step releases `_flow · _dt / emission_dt_ref` particles per source (= `_flow` per step at `_dt = 0.3`).
    - `_number_of_sources`: (int=3) ammount of sources that liberate oil to the environment
    - `_max_contamination_value`: (int=5) value to dictate maximun concentration of contamination per node. No more than "_max_contamination_value" particles can be present in the same node.
    - `_source_fuel`: (int=1000) maximum fuel a source can liberate to the environment, once a source liberates "_source_fuel" particles, it stops liberating
    - `_random_seed`: (int=-1) seed to relive experiments, leave to -1 to allways use a different seed
    - `_triangular`: (bool=False) special condition promoving a different particle behaviour where particles moves following a cross shape.

  - Main methods:
    - `step()`
      - Advance the simulation by one time step. The GIL is released during execution for performance.
    - `reset(_seed=-1, _source_points_pos, _v, _u, _wind_speed, _contamination_position)` (all but `_seed` optional: omit them, `None` is not accepted)
      - Reset the simulation state. If `_seed` is provided the RNG is seeded for reproducibility. -1 used for random seed.
      - If `_source_points_pos` (2 × number_of_sources, int) is provided, source points will take chosen position, otherwise seed source points will be taken.
      - `_v`, `_u` (current field, same shape as the map), `_wind_speed` (2) and `_contamination_position` (2 × N initial particles) override the random draws.
    - `get_normalized_density(gaussian: bool = True) -> numpy.ndarray`
      - Return the normalized density matrix (optionally smoothed using a Gaussian kernel if `gaussian=True`).
    - `warmup(n_steps: int)` *(0.1.0)*
      - Run `n_steps` steps in C++ (GIL released) to start the spill before an episode. Identical to calling `step()` `n_steps` times.
    - `set_dynamics(...)` *(0.1.0)*
      - Configure the physics and the optional dynamics (table below). Call it after `reset()`; the options persist across later resets. Every option not passed is set to its default.
    - `set_wind(wind, also_mean=True)` *(0.1.0)*
      - Set the current wind (2,). With `also_mean` it also becomes the prevailing wind of the variable-wind process.
    - `set_source_fuel(fuel: int)` *(0.1.0)*
      - Remaining fuel of the sources. A finite release stops when it runs out; with a large `_flow`, an instantaneous puff.
    - `get_source_position() -> numpy.ndarray` *(0.1.0)*
      - Continuous position (2 × number_of_sources) of the sources (they can move, see `source_speed` / `source_drift`).

  - Read-only attributes (accessible from Python):
    - `source_points` — coordinates of contamination sources (list/array).
    - `contamination_position` — array containing all the discrete particles positions (list/array).
    - `density` — raw hydrocarbon density matrix (not normalized).
    - `x`, `y` — spatial coordinate vectors associated with the grid.
    - `u`, `v` — velocity fields (x/y components) used by the simulation to simulate currents (array).
    - `wind_speed` — wind speed value or field used in the model (array(2)). Changes every step when `wind_variability > 0`.
    - `wind_mean` *(0.1.0)* — prevailing wind of the variable-wind process (array(2)).
    - `step_count`, `remaining_fuel`, `particle_count` *(0.1.0)* — steps since `reset`, fuel left, particles alive.

### `set_dynamics` options

Time is measured in simulator time units (`steps · _dt`) and positions in cells.

| Option | Default | Effect |
|---|---|---|
| `wind_variability` | 0 | Noise intensity σ of an Ornstein–Uhlenbeck wind around the prevailing wind (0 = constant wind, as before). Stationary std per component = σ / √(2θ). |
| `wind_relaxation` | 1/30 | Reversion rate θ to the prevailing wind; correlation time 1/θ time units. |
| `wind_max` | √2 | Cap on \|wind\|: the envelope of the original `U(-1,1)²` draw, so the wind stays in the usual range. |
| `source_speed` | 0 | Self-propelled source speed (cells / time unit). The source bounces back from land and from the map edge. |
| `source_turn` | 0 | Heading random walk of the moving source (rad / √time). |
| `source_drift` | 0 | The source drifts with `kw·wind + kc·current` (1 = like an oil particle). Can be combined with `source_speed`. |
| `decay_rate` | 0 | Decay / evaporation hazard λ (1 / time unit): each particle disappears with probability `1 − exp(−λ·dt)` per step (half-life ln 2 / λ). |
| `evaporable_fraction` | 1 | Fraction of new particles that can decay; the rest is a persistent residue. |
| `diffusion_dt_ref` | 0.3 | Random walk × `sqrt(_dt / diffusion_dt_ref)`. ≤ 0: legacy (per-step random walk). |
| `emission_dt_ref` | 0.3 | Emission of `_flow · _dt / emission_dt_ref` particles per step. ≤ 0: legacy (`_flow` per step). |
| `fix_removal_skip` | True | A particle leaving the grid no longer makes the next one skip its move. False: legacy. |

The optional dynamics (wind, source, decay) use their own random generator, seeded from the episode seed. Turning them on or off does not change the random sequence of the particles.

### Physics

Each step, every particle moves by `dt · (kw · wind + kc · current(cell)) + U(-2, 2) · gamma · sqrt(dt / dt_ref)`. Particles that hit land stay in place, particles that leave the grid are removed, and a cell holding `_max_contamination_value` particles pushes newcomers to a free neighbour. Each source releases `_flow · dt / dt_ref` particles per step while it has fuel. The density map is the particle count per cell divided by `_max_contamination_value`, by default smoothed with a 5×5 Gaussian kernel. The currents are a static random cellular field and the initial wind is drawn `U(-1, 1)` per component, both from the seed.

### Legacy physics (reproduce 0.0.x data)

Since 0.1.0 the defaults make the physics independent of `_dt`. To reproduce data generated with 0.0.x exactly (bit for bit, same seed):

```py
sim.reset(_seed=seed)
sim.set_dynamics(diffusion_dt_ref=-1, emission_dt_ref=-1, fix_removal_skip=False)
```

Notes:
- For exact signatures and default values check `src/utils/pywrap/pywrap.cpp`.
- Exposed structures (for example Eigen matrices) are converted to/from NumPy arrays via `pybind11/eigen.h`.

## Tests

The CI test script is `tests/test.py` (also referenced in `pyproject.toml`), accesible from [github](https://github.com/AloePacci/cpp_oil_simulator). It checks reproducibility, the time-step-consistent physics against legacy mode, and every optional dynamic (about 1 s). Run it with:

```bash
python tests/test.py
```

## Contributing

Contributions and bug reports are welcome. Please open issues or pull requests on the GitHub repository: https://github.com/AloePacci/cpp_oil_simulator

## Authors

Alejandro Casado Pérez — acasado4@us.es

Acknowledgements: Samuel Yanes Luis (help and guidance).

## License

This project is provided under the terms in `LICENSE`.
