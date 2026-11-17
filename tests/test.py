"""Smoke + regression tests of the oilspillsim extension (run by CI: `python tests/test.py`).

Only numpy is needed. Every check is an assert; the script exits with an error on the first failure.
"""
import math

import numpy as np

import oilspillsim
from oilspillsim import oilspillsim as SIMULATOR

H = W = 60
BASE = np.ones((H, W), np.int32)
BASE[:5, :] = BASE[-5:, :] = BASE[:, :5] = BASE[:, -5:] = 0          # land border
OPEN = np.ones((H, W), np.int32)                                       # no land: particles can leave the grid
LEGACY = dict(diffusion_dt_ref=-1, emission_dt_ref=-1, fix_removal_skip=False)


def make(seed=3, base=BASE, **kw):
    args = dict(_dt=0.3, _kw=1.0, _kc=1.0, _gamma=0.3, _flow=5, _number_of_sources=1,
                _max_contamination_value=5, _source_fuel=100000, _random_seed=seed)
    args.update(kw)
    s = SIMULATOR(base, **args)
    s.reset(_seed=seed)
    return s


def test_basic():
    s = make()
    for _ in range(50):
        s.step()
    d = np.array(s.get_normalized_density())
    assert d.shape == (H, W) and np.isfinite(d).all() and d.max() > 0
    assert s.step_count == 50 and s.particle_count == np.array(s.contamination_position).shape[1]
    assert isinstance(oilspillsim.__version__, str)


def test_seed_reproducible():
    a, b = make(seed=11), make(seed=11)
    a.warmup(100)
    for _ in range(100):
        b.step()
    assert np.array_equal(np.array(a.contamination_position), np.array(b.contamination_position)), "warmup(n) != n x step()"


def test_optional_dynamics_do_not_touch_main_rng():
    """With the optional dynamics at 0, set_dynamics() leaves the trajectory unchanged."""
    a, b = make(seed=5), make(seed=5)
    b.set_dynamics(wind_variability=0.0, source_speed=0.0, decay_rate=0.0)
    a.warmup(80)
    b.warmup(80)
    assert np.array_equal(np.array(a.contamination_position), np.array(b.contamination_position))


def test_emission_per_unit_time():
    emitted = {}
    for dt in (0.15, 0.3, 0.6, 1.2):
        s = make(_dt=dt, _kw=0, _kc=0, _gamma=0, _source_fuel=10**7)
        s.warmup(int(round(60 / dt)))
        emitted[dt] = 10**7 - s.remaining_fuel
    assert max(emitted.values()) - min(emitted.values()) <= 2, emitted          # ~1000 whatever dt
    legacy = {}
    for dt in (0.3, 1.2):
        s = make(_dt=dt, _kw=0, _kc=0, _gamma=0, _source_fuel=10**7)
        s.set_dynamics(**LEGACY)
        s.warmup(int(round(60 / dt)))
        legacy[dt] = 10**7 - s.remaining_fuel
    assert legacy[0.3] > 3 * legacy[1.2], legacy                                 # legacy: per STEP


def test_diffusion_scales_with_sqrt_dt():
    spread = {}
    for dt in (0.15, 1.2):
        v = []
        for seed in range(3):
            s = make(seed=seed, _dt=dt, _kw=0, _kc=0, _flow=1e6, _source_fuel=300, _max_contamination_value=10**6)
            s.warmup(int(round(30 / dt)))
            v.append(np.array(s.contamination_position).std(1).mean())
        spread[dt] = np.mean(v)
    assert abs(spread[1.2] / spread[0.15] - 1) < 0.15, spread                  # legacy ratio would be ~0.35


def test_no_skipped_particle_after_removal():
    s = make(seed=4, base=OPEN, _kw=3, _kc=0, _gamma=0.5, _flow=40, _max_contamination_value=10**6, _source_fuel=10**7)
    s.set_wind(np.array([0.0, 1.0]))
    s.warmup(80)
    for _ in range(10):
        before = {tuple(c) for c in np.array(s.contamination_position).T}
        s.step()
        after = np.array(s.contamination_position)
        assert not any(tuple(c) in before for c in after.T), "a particle skipped its move"


def test_variable_wind():
    s = make(seed=7)
    s.set_dynamics(wind_variability=0.08)
    w = []
    for _ in range(600):
        s.step()
        w.append(np.array(s.wind_speed))
    w = np.array(w)
    assert np.linalg.norm(w, axis=1).max() <= math.sqrt(2) + 1e-9
    assert w.std(0).min() > 0.05
    assert np.allclose(np.array(s.wind_mean), w.mean(0), atol=0.4)


def test_moving_source():
    for kw in (dict(source_speed=0.5, source_turn=0.3), dict(source_drift=1.0)):
        s = make(seed=7)
        s.set_dynamics(**kw)
        p0 = np.array(s.get_source_position())[:, 0].copy()
        for _ in range(200):
            s.step()
            r, c = np.array(s.source_points)[:, 0]
            assert BASE[r, c] == 1, "source on land"
        assert np.linalg.norm(np.array(s.get_source_position())[:, 0] - p0) > 1.0, kw


def test_decay_half_life():
    s = make(seed=2, _kw=0, _kc=0, _flow=1e6, _source_fuel=2000, _max_contamination_value=10**6)
    s.step()
    n0 = s.particle_count
    s.set_dynamics(decay_rate=0.05)
    s.warmup(int(round(math.log(2) / 0.05 / 0.3)))
    assert abs(s.particle_count / n0 - 0.5) < 0.06, (n0, s.particle_count)
    s = make(seed=2, _kw=0, _kc=0, _flow=1e6, _source_fuel=2000, _max_contamination_value=10**6)
    s.set_dynamics(decay_rate=1.0, evaporable_fraction=0.3)
    s.warmup(200)
    assert abs(s.particle_count / 2000 - 0.7) < 0.05, s.particle_count            # persistent residue


def test_finite_release_and_reset_args():
    s = make()
    s.set_source_fuel(100)
    s.warmup(100)
    assert s.remaining_fuel == 0 and s.particle_count <= 101
    s = SIMULATOR(BASE, _number_of_sources=1, _random_seed=1)
    s.reset(_seed=1, _source_points_pos=np.array([[30], [30]], np.int32), _wind_speed=np.array([0.5, -0.2]),
            _contamination_position=np.array([[30.0, 31.0], [30.0, 30.0]]))
    s.set_dynamics(decay_rate=0.05, evaporable_fraction=0.5, wind_variability=0.05)
    s.warmup(50)
    assert np.allclose(np.array(s.wind_mean), [0.5, -0.2])


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"all {len(tests)} tests passed (oilspillsim {oilspillsim.__version__})")
