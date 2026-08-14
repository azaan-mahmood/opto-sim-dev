"""Tests for the SS-TDM DFB laser model (Kim, Chung & Lee 2000, JQE 36(7)).

PARKED -- this file is rewritten once the DFB feature work is finished.

It was written against an earlier API that no longer exists.  It imports
``LaserParams`` and ``_iir_lowpass``, which the module never had, plus a
module ``src.lasers.drive`` that is now ``src.lasers.laser_driver``; it
calls ``_coupling_matrix(kappa)`` where the module has ``_coupling()``
returning scalars, and reads ``res.power_right``/``res.field_right`` where
``SimResult`` has ``P_right``/``E_right``.  Five of its tests exercise a
finite gain-bandwidth IIR filter that is deliberately not implemented, so
there is nothing to port there at all.

The bodies are left untouched so the rewrite can pick over them.  Without
the skip the file fails at import, and a collection failure makes pytest
abandon the whole run, so a bare ``pytest`` at the repository root collects
nothing -- which reads like a clean run and is not one.
"""

import numpy as np
import pytest
import warnings

pytest.skip("rewritten after the DFB feature work; see the module docstring",
            allow_module_level=True)

from src.lasers.dfblaser import (DFBLaser, LaserParams, _iir_lowpass)
from src.lasers import drive

DNU_TEST = 0.2e12          # Hz, filter bandwidth used in the frequency tests
DT_TEST = 0.185e-12        # s, N=40 device step at 600 um
POLE_TEST = float(np.exp(-np.pi * DNU_TEST * DT_TEST))


def run_steady(laser, current, t_end=20e-9, record_every=4):
    """Steady right-facet power over the last 25% of a run (W)."""
    res = laser.simulate(current, t_end=t_end, record_every=record_every)
    n = max(1, len(res.power_right) // 4)
    return res, float(np.mean(res.power_right[-n:]))


# --- gain-bandwidth filter ---------------------------------------------------

def test_gain_filter_dc_gain_one():
    x = np.full(8, 3.0 + 2.0j)
    state = np.zeros(8, dtype=complex)
    for _ in range(200):
        y = _iir_lowpass(x, state, POLE_TEST)
        state = y
    assert np.allclose(y, x, rtol=1e-6)


def test_gain_filter_three_db_bandwidth():
    n = 4000
    t = np.arange(n) * DT_TEST
    x = np.exp(2j * np.pi * 100e9 * t)          # designed -3 dB point
    state = np.zeros(1, dtype=complex)
    out = []
    for xn in x:
        y = _iir_lowpass(xn, state, POLE_TEST)
        state = y
        out.append(y)
    out = np.array(out)
    gain = np.mean(np.abs(out[n // 2:]) ** 2)   # |x|^2 = 1
    assert 0.40 < gain < 0.60, f"gain at -3dB point: {gain:.3f}"


def test_gain_filter_flat_in_band():
    n = 2000
    t = np.arange(n) * DT_TEST
    x = np.exp(2j * np.pi * 10e9 * t)           # 1/20 of the -3 dB point
    state = np.zeros(1, dtype=complex)
    out = []
    for xn in x:
        y = _iir_lowpass(xn, state, POLE_TEST)
        state = y
        out.append(y)
    out = np.array(out)
    gain = np.mean(np.abs(out[n // 2:]) ** 2)
    assert 0.95 < gain < 1.05, f"in-band gain: {gain:.3f}"


# --- coupling matrix power conservation (paper Eq. 17) -----------------------

def test_coupling_matrix_power_conservation():
    rng = np.random.default_rng(7)
    laser = DFBLaser(LaserParams(), n_sections=24, seed=7)
    for _ in range(50):
        kappa = (rng.normal(0, 60, 24) + 1j * rng.normal(0, 60, 24)).astype(complex)
        F = rng.standard_normal(24) + 1j * rng.standard_normal(24)
        R = rng.standard_normal(24) + 1j * rng.standard_normal(24)
        sech, cf, cr = laser._coupling_matrix(kappa)
        F_out = sech * F + cf * R
        R_out = cr * F + sech * R
        pin = np.sum(np.abs(F) ** 2 + np.abs(R) ** 2)
        pout = np.sum(np.abs(F_out) ** 2 + np.abs(R_out) ** 2)
        assert pout == pytest.approx(pin, rel=1e-10)


# --- convergence criterion (paper Fig. 5: kappa*dz < 0.2) ---------------------

def test_convergence_warning_fires_below_5_kappa_L():
    p = LaserParams(kappa_i0=50.0, kappa_g0=0.0)   # kappa*L = 3.0, need N>=15
    with pytest.warns(UserWarning, match="5\\*kappa\\*L"):
        DFBLaser(p, n_sections=10)


def test_no_warning_above_5_kappa_L():
    p = LaserParams(kappa_i0=50.0, kappa_g0=0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        DFBLaser(p, n_sections=40)


# --- end-to-end behaviour ----------------------------------------------------

def test_laser_smoke_with_filter():
    p = LaserParams(kappa_i0=30.0, kappa_g0=20.0)
    laser = DFBLaser(p, n_sections=20, seed=42)
    res = laser.simulate(0.05, t_end=4e-9, record_every=2)
    assert np.all(np.isfinite(res.power_right))
    assert np.all(res.power_right >= 0.0)


def test_seed_determinism():
    p = LaserParams(kappa_i0=30.0, kappa_g0=20.0)
    r1 = DFBLaser(p, n_sections=20, seed=42).simulate(0.05, t_end=2e-9)
    r2 = DFBLaser(p, n_sections=20, seed=42).simulate(0.05, t_end=2e-9)
    assert np.array_equal(r1.power_right, r2.power_right)
    assert np.array_equal(r1.field_right, r2.field_right)


def test_40nm_filter_equivalent_to_flat_in_band():
    # The paper's ~40 nm gain profile (5 THz) is wider than the simulation
    # band (2.7 THz Nyquist at N=40); per Sec. III it must leave the single
    # mode essentially untouched.  The gain-excess split keeps the carrier
    # path exact, so the filter is a near-null on the steady power.
    p_flat = LaserParams(kappa_i0=30.0, kappa_g0=20.0, gain_bw_nm=None)
    p_filt = LaserParams(kappa_i0=30.0, kappa_g0=20.0, gain_bw_nm=40.0)
    _, p_flat = run_steady(DFBLaser(p_flat, n_sections=20, seed=42), 0.1)
    _, p_filt = run_steady(DFBLaser(p_filt, n_sections=20, seed=42), 0.1)
    rel = abs(p_filt - p_flat) / p_flat
    assert rel < 0.10, f"filter changed steady power by {rel * 100:.1f}%"


def test_narrow_gain_filter_reduces_power():
    # A representable bandwidth (0.8 nm ~ 100 GHz at 1550 nm) attenuates the
    # DFB modes off line center, raising threshold and lowering power.
    p_flat = LaserParams(kappa_i0=30.0, kappa_g0=20.0, gain_bw_nm=None)
    p_filt = LaserParams(kappa_i0=30.0, kappa_g0=20.0, gain_bw_nm=0.8)
    _, p_flat = run_steady(DFBLaser(p_flat, n_sections=20, seed=42), 0.1)
    _, p_filt = run_steady(DFBLaser(p_filt, n_sections=20, seed=42), 0.1)
    assert p_filt < 0.99 * p_flat, f"narrow filter did not reduce power: {p_filt} vs {p_flat}"
    assert p_filt > 1e-3 * p_flat, f"narrow filter killed lasing entirely: {p_filt} vs {p_flat}"


# --- drive demo module -------------------------------------------------------

def test_drive_module_imports():
    assert callable(drive.cw(0.1))
    assert callable(drive.gain_switched(0.04, 0.22, 500e-12, 150e-12))
    gs = drive.gain_switched(0.04, 0.22, 100e-9, 30e-9)
    assert gs(5e-9) == pytest.approx(0.22)     # inside the pulse
    assert gs(60e-9) == pytest.approx(0.04)    # between pulses
