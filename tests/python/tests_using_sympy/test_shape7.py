"""
Test spectral_wave_data shape 7 against the analytic reference in shape_7.py.

Covers:
  - All supported kinematic quantities: phi, phi_t, grad_phi, pressure,
    elev, elev_t, grad_elev, grad_elev_2nd.
  - phi continuity across sigma-layer boundaries.
  - phi continuity at zref (between layered and below-zref branches).
  - Sigma capping: phi(z > zeta) == phi(z = zeta).
  - get() accessor for nlayers.
"""

import sys
import os
import math

import pytest

import numpy as np

import shape_7
import corsys
from tfun import Tfun

from spectral_wave_data import SpectralWaveData
from test_utils import should_run_quick_tests

assert sys.version_info > (3, 4)

# ---------------------------------------------------------------------------
# Shape-7 parameters used for all tests
# ---------------------------------------------------------------------------
NLAYERS = 3
SIG = [0.0, 0.5, 1.0]
# zref well below wave troughs; d=3.0 → d_eff = 3.0 + (-1.5) = 1.5 > 0
ZREF = -1.5
# Scale hfuns small so zeta << |zref| (prevents degenerate sigma columns)
HSCALE = 0.2

# Shape 7 reconstructs time derivatives via 4-step 2nd-order finite differences.
# These are exact ONLY for polynomials of degree <= 2.  Use poly_order=2 for
# ALL Tfun objects so that h_cur / ht_cur / c_cur / ct_cur are numerically exact.
POLY_ORDER = 2


def _create_tfuns(n_swd, seed):
    """Create n_swd+1 degree-2 Tfun objects (exact with 2nd-order FD)."""
    np.random.seed(seed)
    return [Tfun(POLY_ORDER) for _ in range(n_swd + 1)]

# ---------------------------------------------------------------------------
# Test matrix (subset of shape-2 matrix for speed)
# ---------------------------------------------------------------------------
impls = [1]
dks = [0.3]
ipols = [0, 1, 2]
nswds = [3]
seeds = [0]
ds = [3.0]
appsys = [
    (0.0,  0.0,  0.0,   0.0),
    (13.8, 0.0,  0.0,   0.0),
    (0.0,  0.0,  0.34,  0.0),
    (-5.3, 18.2, 0.0,  37.4),
]  # 4 coordinate systems (vs 8 in test_shape2)

waves_param = [
    {'impl': impl, 'dk': dk, 'd': d, 'ipol': ipol, 'nswd': nswd,
     'seed': seed, 'x0': x0, 'y0': y0, 't0': t0, 'beta': beta}
    for impl in impls
    for dk in dks
    for d in ds
    for ipol in ipols
    for nswd in nswds
    for seed in seeds
    for x0, y0, t0, beta in appsys
]


def dicts2ids(ds):
    return ['_'.join(f'{k}:{v}' for k, v in d.items()) for d in ds]


@pytest.fixture(scope='function', params=waves_param, ids=dicts2ids(waves_param))
def make_waves(request, tmp_path_factory):
    impl  = request.param['impl']
    dk    = request.param['dk']
    d     = request.param['d']
    ipol  = request.param['ipol']
    nswd  = request.param['nswd']
    seed  = request.param['seed']
    x0    = request.param['x0']
    y0    = request.param['y0']
    t0    = request.param['t0']
    beta  = request.param['beta']

    mysys = corsys.CorSys(x0, y0, t0, beta)

    # Elevation: degree-2 Tfuns scaled small to keep zeta well above zref=-1.5.
    # Degree 2 ensures the 2nd-order finite-difference reconstruction in the
    # shape-7 Fortran is numerically exact.
    hfuns_raw = _create_tfuns(nswd, seed)
    hfuns = [h.scale(HSCALE) for h in hfuns_raw]

    # Independent degree-2 potential functions per sigma-layer
    cfuns_layers = []
    for m in range(NLAYERS):
        cfuns_m = _create_tfuns(nswd, seed * 100 + m + 1)
        cfuns_layers.append(cfuns_m)

    swd_anal = shape_7.Shape7(dk, nswd, d, ZREF, NLAYERS, SIG,
                              cfuns_layers, hfuns, mysys)

    tmpdir = str(tmp_path_factory.mktemp('swd7'))
    fname = (
        f'shp7_impl{impl}_dk{dk}_d{d}_ipol{ipol}'
        f'_nswd{nswd}_seed{seed}_x0{x0}_y0{y0}_t0{t0}_beta{beta}.swd'
    )
    file_swd = os.path.join(tmpdir, fname)
    swd_anal.write_swd(file_swd, dt=0.1, nsteps=11)
    swd_anal.check_swd_meta(file_swd, nswd)

    swd_num = SpectralWaveData(file_swd, x0, y0, t0, beta,
                               rho=1025.0, impl=impl, dc_bias=True)
    yield swd_anal, swd_num
    swd_num.close()


# ---------------------------------------------------------------------------
# Evaluation points
# ---------------------------------------------------------------------------
xs = [0.0, 8.4]
ys = [0.0, -5.2]
# z_below < ZREF: shape-2 extrapolation path
# z_low, z_mid:  inside sigma column (different sigma values)
z_below = ZREF - 1.5   # = -3.0
z_low   = ZREF + 0.3   # = -1.2  (low sigma ~0.14 for typical H)
z_mid   = 0.0           # middle of sigma column (sigma ~0.7 for typical H)
zs_all  = [z_below, z_low, z_mid]
ts = [0.2, 0.23, 0.52]


# ---------------------------------------------------------------------------
# Main kinematic test
# ---------------------------------------------------------------------------

def test_waves(make_waves):
    """All supported quantities match the analytic reference within tolerance."""
    swd_anal, swd_num = make_waves
    quick_check = should_run_quick_tests()

    for t in ts:
        swd_num.update_time(t)

        for x in xs:
            for y in ys:
                # --- Surface elevation (z-independent) ---
                elev_anal = swd_anal.elev(x, y, t)
                elev_num  = swd_num.elev(x, y)
                assert math.isclose(elev_num, elev_anal, rel_tol=1e-4, abs_tol=1e-4), \
                    f'elev mismatch at x={x}, y={y}, t={t}'

                if not quick_check:
                    elev_t_anal = swd_anal.elev_t(x, y, t)
                    elev_t_num  = swd_num.elev_t(x, y)
                    assert math.isclose(elev_t_num, elev_t_anal,
                                        rel_tol=1e-4, abs_tol=1e-4), \
                        f'elev_t mismatch at x={x}, y={y}, t={t}'

                    ge_anal = swd_anal.grad_elev(x, y, t)
                    ge_num  = swd_num.grad_elev(x, y)
                    assert math.isclose(ge_num.x, ge_anal['x'], rel_tol=1e-4, abs_tol=1e-4)
                    assert math.isclose(ge_num.y, ge_anal['y'], rel_tol=1e-4, abs_tol=1e-4)
                    assert math.isclose(ge_num.z, ge_anal['z'], rel_tol=1e-4, abs_tol=1e-4)

                    ge2_anal = swd_anal.grad_elev_2nd(x, y, t)
                    ge2_num  = swd_num.grad_elev_2nd(x, y)
                    assert math.isclose(ge2_num.xx, ge2_anal['xx'], rel_tol=1e-4, abs_tol=1e-4)
                    assert math.isclose(ge2_num.xy, ge2_anal['xy'], rel_tol=1e-4, abs_tol=1e-4)
                    assert math.isclose(ge2_num.yy, ge2_anal['yy'], rel_tol=1e-4, abs_tol=1e-4)

                for z in zs_all:
                    # --- Potential ---
                    phi_anal = swd_anal.phi(x, y, z, t)
                    phi_num  = swd_num.phi(x, y, z)
                    assert math.isclose(phi_num, phi_anal, rel_tol=1e-4, abs_tol=1e-4), \
                        f'phi mismatch at x={x}, y={y}, z={z}, t={t}: ' \
                        f'anal={phi_anal:.6g}, num={phi_num:.6g}'

                    if quick_check:
                        continue

                    # --- Euler time derivative (chain-rule corrected) ---
                    phi_t_anal = swd_anal.phi_t(x, y, z, t)
                    phi_t_num  = swd_num.phi_t(x, y, z)
                    assert math.isclose(phi_t_num, phi_t_anal, rel_tol=1e-4, abs_tol=1e-4), \
                        f'phi_t mismatch at x={x}, y={y}, z={z}, t={t}: ' \
                        f'anal={phi_t_anal:.6g}, num={phi_t_num:.6g}'

                    # --- Particle velocity ---
                    gp_anal = swd_anal.grad_phi(x, y, z, t)
                    gp_num  = swd_num.grad_phi(x, y, z)
                    assert math.isclose(gp_num.x, gp_anal['x'], rel_tol=1e-4, abs_tol=1e-4), \
                        f'grad_phi.x mismatch at z={z}, t={t}'
                    assert math.isclose(gp_num.y, gp_anal['y'], rel_tol=1e-4, abs_tol=1e-4), \
                        f'grad_phi.y mismatch at z={z}, t={t}'
                    assert math.isclose(gp_num.z, gp_anal['z'], rel_tol=1e-4, abs_tol=1e-4), \
                        f'grad_phi.z mismatch at z={z}, t={t}'

                    # --- Pressure ---
                    prs_anal = swd_anal.pressure(x, y, z, t, rho=1025.0, grav=9.81)
                    prs_num  = swd_num.pressure(x, y, z)
                    assert math.isclose(prs_num, prs_anal, rel_tol=1e-3), \
                        f'pressure mismatch at x={x}, y={y}, z={z}, t={t}: ' \
                        f'anal={prs_anal:.6g}, num={prs_num:.6g}'


# ---------------------------------------------------------------------------
# Shape-7-specific structural tests
# ---------------------------------------------------------------------------

def test_layer_continuity(make_waves):
    """phi must be continuous across every interior sigma-layer boundary.

    At each interior layer position z_m = zref + sig[m] * H the value
    obtained just below and just above must agree to within a loose tolerance
    (slight numerical difference is expected from the piecewise-linear
    formulation at the exact layer boundary).
    """
    swd_num = make_waves[1]
    eps = 1e-4

    swd_num.update_time(ts[0])
    for x in xs[:1]:
        for y in ys[:1]:
            eta = swd_num.elev(x, y)
            H = eta - ZREF
            for m in range(1, NLAYERS - 1):  # interior layers only
                z_layer = ZREF + SIG[m] * H
                phi_below = swd_num.phi(x, y, z_layer - eps)
                phi_above = swd_num.phi(x, y, z_layer + eps)
                assert math.isclose(phi_below, phi_above, rel_tol=1e-3, abs_tol=1e-5), \
                    (f'phi discontinuity at interior layer m={m}, '
                     f'z_layer={z_layer:.4f}: '
                     f'phi(below)={phi_below:.6g}, phi(above)={phi_above:.6g}')


def test_zref_continuity(make_waves):
    """phi must be continuous at z = zref (layered branch vs shape-2 branch).

    By construction, the shape-2 extrapolation at z'=0 gives exactly the
    bottom-layer value (Z_j(0)=1), so phi must agree on both sides of zref.
    """
    swd_num = make_waves[1]
    eps = 1e-4

    swd_num.update_time(ts[0])
    for x in xs[:1]:
        for y in ys[:1]:
            phi_above = swd_num.phi(x, y, ZREF + eps)   # sigma-layer branch, sigma ≈ 0
            phi_below = swd_num.phi(x, y, ZREF - eps)   # shape-2 extrapolation branch
            assert math.isclose(phi_above, phi_below, rel_tol=1e-3, abs_tol=1e-5), \
                (f'phi discontinuity at zref={ZREF}: '
                 f'above={phi_above:.6g}, below={phi_below:.6g}')


def test_sigma_cap(make_waves):
    """phi and grad_phi at z > zeta must equal the surface (sigma=1) values.

    When sigma > 1 it is clamped to 1, so any z above the free surface
    must return the same kinematics as z = zeta.
    """
    swd_num = make_waves[1]

    swd_num.update_time(ts[0])
    for x in xs[:1]:
        for y in ys[:1]:
            eta = swd_num.elev(x, y)
            z_at    = eta           # sigma = 1.0 exactly
            z_above = eta + 2.0     # sigma > 1 → clamped to 1.0

            phi_at    = swd_num.phi(x, y, z_at)
            phi_above = swd_num.phi(x, y, z_above)
            assert math.isclose(phi_above, phi_at, rel_tol=1e-7, abs_tol=1e-9), \
                f'sigma cap (phi): eta={eta:.4f}, phi_at={phi_at:.6g}, phi_above={phi_above:.6g}'

            gp_at    = swd_num.grad_phi(x, y, z_at)
            gp_above = swd_num.grad_phi(x, y, z_above)
            assert math.isclose(gp_above.x, gp_at.x, rel_tol=1e-7, abs_tol=1e-9)
            assert math.isclose(gp_above.z, gp_at.z, rel_tol=1e-7, abs_tol=1e-9)


def test_below_zref_matches_shape2_at_layer0(make_waves):
    """Just below zref, the potential must match the analytic reference formula.

    This exercises the shape-2-style branch with the bottom-layer coefficients
    and effective depth d_eff = d + zref.
    """
    swd_anal, swd_num = make_waves

    swd_num.update_time(ts[0])
    for x in xs[:1]:
        for y in ys[:1]:
            for z in [ZREF - 0.5, ZREF - 1.0, ZREF - 1.4]:
                phi_anal = swd_anal.phi(x, y, z, ts[0])
                phi_num  = swd_num.phi(x, y, z)
                assert math.isclose(phi_num, phi_anal, rel_tol=1e-4, abs_tol=1e-4), \
                    f'phi below zref mismatch at z={z}: anal={phi_anal:.6g}, num={phi_num:.6g}'


def test_get_nlayers(make_waves):
    """get('nlayers') must return the correct number of sigma layers."""
    swd_num = make_waves[1]
    assert swd_num.get('nlayers') == NLAYERS
