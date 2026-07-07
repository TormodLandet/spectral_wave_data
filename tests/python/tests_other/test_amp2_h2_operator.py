"""
End-to-end test for amp=2 lazy H2 implementation (spectral_wave_data_shape_1_or_2_impl_7).

Strategy:
  1. Write a synthetic amp=2 shape-1 SWD file containing a single Airy wave
     (deep-water linear wave theory, small amplitude so H2 ≈ linear).
  2. Open the file with SpectralWaveData (impl=0 → auto-selects H2 impl).
  3. Verify:
     a. Wave elevation matches analytical input exactly (uses h-spectral directly).
     b. Surface potential phi at z=0 is close to linear theory.
     c. Kinematics decay with depth (deep-water exponential profile).
     d. Sequential stepping gives same result as seeking.
     e. Shape 2 (finite depth) also works.
     f. Non-power-of-two nx raises an error.
     g. get() accessors return correct values.

Note: the H2 operator is a nonlinear method, so for large-amplitude waves it
      will differ from linear theory.  These tests use A/lambda ≈ 0.003 (very
      small steepness) so the H2 correction is negligible and linear theory is
      a valid reference.

SWD spectral coefficient convention used here (matches fft_real_to_swd in
swd_fft.f90):
    h_swd(j) = (2/n) * conj(rfft(eta)[j])   for j in 1..n//2-1
    h_swd(0) = (1/n) * rfft(eta)[0]          (DC)
    h_swd(n//2) = (1/n) * rfft(eta)[n//2]    (Nyquist, real)
"""

import os
import struct
import datetime
import math

import numpy as np
import pytest

from spectral_wave_data import SpectralWaveData

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
GRAV = 9.81


# ---------------------------------------------------------------------------
# SWD file writer (amp=2 shape 1 or 2)
# ---------------------------------------------------------------------------


def _to_complex64_bytes(arr):
    """Pack a complex128 ndarray as interleaved float32 little-endian (SWD format)."""
    arr_c64 = arr.astype(np.complex64)
    buf = np.empty(len(arr_c64) * 2, dtype=np.float32)
    buf[0::2] = arr_c64.real
    buf[1::2] = arr_c64.imag
    return buf.astype("<f4").tobytes()


def _real_to_swd_coeff(eta_real):
    """Convert a real spatial array to SWD complex spectral coefficients."""
    n_real = len(eta_real)
    rfft_out = np.fft.rfft(eta_real)  # length n_real//2 + 1
    n_swd = n_real // 2  # highest index = n
    h = np.zeros(n_swd + 1, dtype=complex)
    h[0] = rfft_out[0] / n_real
    h[1:-1] = 2.0 * np.conj(rfft_out[1:-1]) / n_real
    h[n_swd] = rfft_out[n_swd] / n_real
    return h


def write_amp2_swd(
    filename, n_swd, dk, nsteps, dt, h_list, ht_list, c_list, ct_list, depth=None, grav=GRAV
):
    """
    Write a minimal amp=2 SWD file (shape 1 for deep water, shape 2 for finite depth).

    Parameters
    ----------
    n_swd   : int   — highest spectral index (nx = 2*n_swd, must be power of two)
    dk      : float — wavenumber spacing
    nsteps  : int   — number of time steps
    dt      : float — time step size
    h_list  : list of ndarray(n_swd+1, complex)  — elevation coefficients per step
    ht_list : list of ndarray(n_swd+1, complex)  — d/dt elevation per step
    c_list  : list of ndarray(n_swd+1, complex)  — surface potential per step
    ct_list : list of ndarray(n_swd+1, complex)  — d/dt potential per step
    depth   : float or None  — water depth; None → deep water (shape 1)
    """
    shp = 1 if depth is None else 2

    prog = b"SWD_test_amp2".ljust(30)[:30]
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S").encode()[:20].ljust(20)
    cid = b"Amp=2 test file\x00"

    with open(filename, "wb") as f:
        # --- SWD header ---
        f.write(struct.pack("<f", 37.0221))  # magic
        f.write(struct.pack("<i", 100))  # fmt
        f.write(struct.pack("<i", shp))  # shp
        f.write(struct.pack("<i", 2))  # amp = 2
        f.write(prog)  # cprog (30 bytes)
        f.write(now)  # cdate (20 bytes)
        f.write(struct.pack("<i", len(cid)))  # nid
        f.write(cid)  # cid
        f.write(struct.pack("<f", grav))  # grav
        f.write(struct.pack("<f", 1.0))  # lscale
        f.write(struct.pack("<i", 0))  # nstrip
        f.write(struct.pack("<i", nsteps))  # nsteps
        f.write(struct.pack("<f", dt))  # dt
        f.write(struct.pack("<i", -1))  # order (fully nonlinear)
        f.write(struct.pack("<i", n_swd))  # n
        f.write(struct.pack("<f", dk))  # dk
        if shp == 2:
            f.write(struct.pack("<f", depth))  # d (shape 2 only)

        # --- Temporal data ---
        for istep in range(nsteps):
            f.write(_to_complex64_bytes(h_list[istep]))
            f.write(_to_complex64_bytes(ht_list[istep]))
            f.write(_to_complex64_bytes(c_list[istep]))
            f.write(_to_complex64_bytes(ct_list[istep]))


# ---------------------------------------------------------------------------
# Analytical Airy wave helpers
# ---------------------------------------------------------------------------


def _airy_coefficients(n_swd, dk, A, nsteps, dt, depth=None):
    """
    Generate h/ht/c/ct SWD coefficient arrays for a single mode-1 Airy wave.

    Deep water:  omega = sqrt(g * k1), phi_surface = (omega/k1)*A*sin(kx-omega*t)
    Finite depth: omega = sqrt(g * k1 * tanh(k1*depth))
    """
    k1 = dk  # fundamental wavenumber
    if depth is None:
        omega = math.sqrt(GRAV * k1)
    else:
        omega = math.sqrt(GRAV * k1 * math.tanh(k1 * depth))

    h_list, ht_list, c_list, ct_list = [], [], [], []
    for istep in range(nsteps):
        t = istep * dt
        # h_swd(1) = A * exp(+i*omega*t) so that
        #   Re{h_swd(1)*exp(-i*k1*x)} = A*cos(k1*x - omega*t)
        h = np.zeros(n_swd + 1, dtype=complex)
        ht = np.zeros(n_swd + 1, dtype=complex)
        c = np.zeros(n_swd + 1, dtype=complex)
        ct = np.zeros(n_swd + 1, dtype=complex)

        h[1] = A * np.exp(1j * omega * t)
        ht[1] = 1j * omega * h[1]  # d(eta)/dt = A*omega*sin(kx - omega*t)

        # c_swd(1) s.t. Re{c_swd(1)*exp(-i*k1*x)} = (omega/k1)*A*sin(k1*x - omega*t)
        # Derivation: sin(k1*x - omega*t) = Re{-i*exp(-i*(k1*x - omega*t))}
        #   = Re{-i*exp(i*omega*t)*exp(-i*k1*x)} -> c_swd(1) = +i*(omega/k1)*h[1]
        c[1] = 1j * (omega / k1) * h[1]
        ct[1] = 1j * omega * c[1]  # d(c)/dt

        h_list.append(h)
        ht_list.append(ht)
        c_list.append(c)
        ct_list.append(ct)

    return h_list, ht_list, c_list, ct_list, omega, k1


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def airy_deep_water_swd(tmp_path_factory):
    """Shape-1 amp=2 SWD file with a single deep-water Airy wave."""
    tmp = str(tmp_path_factory.mktemp("amp2_deep"))
    fname = os.path.join(tmp, "airy_deep.swd")
    n_swd = 512  # nx = 1024 (power of two)
    dk = 0.05
    A = 0.15  # small amplitude for linear validity
    nsteps, dt = 8, 0.5
    h_list, ht_list, c_list, ct_list, omega, k1 = _airy_coefficients(
        n_swd, dk, A, nsteps, dt, depth=None
    )
    write_amp2_swd(fname, n_swd, dk, nsteps, dt, h_list, ht_list, c_list, ct_list, depth=None)
    return fname, n_swd, dk, A, nsteps, dt, omega, k1


@pytest.fixture(scope="module")
def airy_finite_depth_swd(tmp_path_factory):
    """Shape-2 amp=2 SWD file with a single finite-depth Airy wave."""
    tmp = str(tmp_path_factory.mktemp("amp2_finite"))
    fname = os.path.join(tmp, "airy_finite.swd")
    n_swd = 512
    dk = 0.05
    A = 0.10
    depth = 20.0
    nsteps, dt = 8, 0.5
    h_list, ht_list, c_list, ct_list, omega, k1 = _airy_coefficients(
        n_swd, dk, A, nsteps, dt, depth=depth
    )
    write_amp2_swd(fname, n_swd, dk, nsteps, dt, h_list, ht_list, c_list, ct_list, depth=depth)
    return fname, n_swd, dk, A, nsteps, dt, omega, k1, depth


# ---------------------------------------------------------------------------
# Test 1: wave elevation matches analytical input
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("x_app", [0.0, 5.0, 12.3])
def test_elev_matches_analytical_deep(airy_deep_water_swd, x_app):
    """
    elev() must match A*cos(k*x - omega*t) — elevation comes directly from h-spectral.
    """
    fname, n_swd, dk, A, nsteps, dt, omega, k1 = airy_deep_water_swd
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        for istep in range(1, nsteps - 1):
            t = istep * dt
            swd.update_time(t)
            eta_num = swd.elev(x_app, 0.0)
            eta_anal = A * math.cos(k1 * x_app - omega * t)
            assert math.isclose(eta_num, eta_anal, rel_tol=1e-4, abs_tol=1e-5), (
                f"elev mismatch at x={x_app}, t={t}: num={eta_num:.6g}, anal={eta_anal:.6g}"
            )
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 2: phi decreases with depth (deep-water exponential decay)
# ---------------------------------------------------------------------------


def test_phi_decays_with_depth(airy_deep_water_swd):
    """
    |phi| must be larger near the surface than deeper — basic deep-water physics.
    """
    fname, n_swd, dk, A, nsteps, dt, omega, k1 = airy_deep_water_swd
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        t = 2 * dt
        swd.update_time(t)
        x_app, y_app = 3.0, 0.0
        eta = swd.elev(x_app, y_app)

        # Sample phi at z just below surface, z=-3, z=-7
        z_levels = [eta - 0.3, -3.0, -7.0]
        phi_vals = [abs(swd.phi(x_app, y_app, z)) for z in z_levels]

        # |phi| must strictly decrease going deeper
        for i in range(len(phi_vals) - 1):
            assert phi_vals[i] > phi_vals[i + 1], (
                f"phi did not decrease: phi({z_levels[i]:.1f})={phi_vals[i]:.4g}, "
                f"phi({z_levels[i + 1]:.1f})={phi_vals[i + 1]:.4g}"
            )
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 3: phi close to linear theory at small amplitude (deep water)
# ---------------------------------------------------------------------------


def test_phi_linear_theory(airy_deep_water_swd):
    """
    For small steepness A*k << 1, H2 should produce phi close to linear theory.
    Linear: phi(x, z, t) = (A*omega/k1) * exp(k1*z) * sin(k1*x - omega*t)
    """
    fname, n_swd, dk, A, nsteps, dt, omega, k1 = airy_deep_water_swd
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        t = 1.5 * dt
        swd.update_time(t)
        x_app, y_app = 4.0, 0.0
        z = -2.0

        phi_num = swd.phi(x_app, y_app, z)
        # Linear theory value
        phi_lin = (A * omega / k1) * math.exp(k1 * z) * math.sin(k1 * x_app - omega * t)

        # Allow generous tolerance because H2 uses finite steps and the field
        # differs slightly even for small amplitudes
        assert math.isclose(phi_num, phi_lin, rel_tol=0.10, abs_tol=1e-3 * A), (
            f"phi too far from linear theory: H2={phi_num:.4g}, linear={phi_lin:.4g}"
        )
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 4: sequential stepping == seeking gives same kinematics
# ---------------------------------------------------------------------------


def test_sequential_equals_seek(airy_deep_water_swd):
    """
    update_time in order t0 < t1 < t2 must give same result as
    opening fresh and jumping directly to t2.
    """
    fname, n_swd, dk, A, nsteps, dt, omega, k1 = airy_deep_water_swd
    x_app, y_app, z = 1.7, 0.0, -1.5
    t_target = 3 * dt

    # Sequential
    swd_seq = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        for istep in range(1, 4):
            swd_seq.update_time(istep * dt)
        phi_seq = swd_seq.phi(x_app, y_app, z)
        eta_seq = swd_seq.elev(x_app, y_app)
    finally:
        swd_seq.close()

    # Direct seek
    swd_seek = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        swd_seek.update_time(t_target)
        phi_seek = swd_seek.phi(x_app, y_app, z)
        eta_seek = swd_seek.elev(x_app, y_app)
    finally:
        swd_seek.close()

    assert math.isclose(phi_seq, phi_seek, rel_tol=1e-5, abs_tol=1e-8), (
        f"phi: sequential={phi_seq:.8g} vs seek={phi_seek:.8g}"
    )
    assert math.isclose(eta_seq, eta_seek, rel_tol=1e-6, abs_tol=1e-9), (
        f"elev: sequential={eta_seq:.8g} vs seek={eta_seek:.8g}"
    )


# ---------------------------------------------------------------------------
# Test 5: get() accessors
# ---------------------------------------------------------------------------


def test_get_accessors_deep(airy_deep_water_swd):
    """get() accessors must return values matching the file header."""
    fname, n_swd, dk, A, nsteps, dt, omega, k1 = airy_deep_water_swd
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        assert swd["shp"] == 1, f"shp={swd['shp']}"
        assert swd["amp"] == 2, f"amp={swd['amp']}"
        assert swd["n"] == n_swd, f"n={swd['n']}"
        assert math.isclose(swd["dk"], dk, rel_tol=1e-5), f"dk={swd['dk']}"
        assert math.isclose(swd["dt"], dt, rel_tol=1e-5), f"dt={swd['dt']}"
        assert swd["d"] < 0, "d should be negative for deep water"
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 6: shape 2 (finite depth) — constructor works and eta matches
# ---------------------------------------------------------------------------


def test_shape2_amp2_elev(airy_finite_depth_swd):
    """Shape-2 amp=2 file must open without error and return correct elevation."""
    fname, n_swd, dk, A, nsteps, dt, omega, k1, depth = airy_finite_depth_swd
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        assert swd["shp"] == 2
        assert swd["amp"] == 2
        assert math.isclose(swd["d"], depth, rel_tol=1e-4)

        for istep in range(1, nsteps - 1):
            t = istep * dt
            swd.update_time(t)
            x_app = 6.0
            eta_num = swd.elev(x_app, 0.0)
            eta_anal = A * math.cos(k1 * x_app - omega * t)
            assert math.isclose(eta_num, eta_anal, rel_tol=1e-4, abs_tol=1e-5), (
                f"elev mismatch at t={t}: {eta_num:.6g} vs {eta_anal:.6g}"
            )
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 7: non-power-of-two nx is now ALLOWED (performance warning only)
# ---------------------------------------------------------------------------


def test_non_power_of_two_nx_opens_ok(tmp_path):
    """
    nx = 2*n no longer needs to be a power of two.  Opening a file with
    n_swd=100 (nx=200) must succeed without raising an error.
    """
    fname = str(tmp_path / "non_pow2_nx.swd")
    n_swd = 100  # nx = 200, not a power of two
    dk = 0.1
    nsteps, dt = 4, 1.0
    A = 0.05
    k1 = dk
    omega = math.sqrt(GRAV * k1)
    h_list, ht_list, c_list, ct_list = [], [], [], []
    for istep in range(nsteps):
        t = istep * dt
        h = np.zeros(n_swd + 1, dtype=complex)
        ht = np.zeros(n_swd + 1, dtype=complex)
        c = np.zeros(n_swd + 1, dtype=complex)
        ct = np.zeros(n_swd + 1, dtype=complex)
        h[1] = A * np.exp(1j * omega * t)
        ht[1] = 1j * omega * h[1]
        c[1] = 1j * (omega / k1) * h[1]
        ct[1] = 1j * omega * c[1]
        h_list.append(h)
        ht_list.append(ht)
        c_list.append(c)
        ct_list.append(ct)
    write_amp2_swd(fname, n_swd, dk, nsteps, dt, h_list, ht_list, c_list, ct_list, depth=None)
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        swd.update_time(dt)
        eta = swd.elev(0.0, 0.0)
        assert math.isfinite(eta)
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 8: env-var SWD_NUM_H2_STEPS override
# ---------------------------------------------------------------------------


def test_env_var_nsteps(airy_deep_water_swd, monkeypatch):
    """
    Setting SWD_NUM_H2_STEPS=5 should result in nlayers = 5+1 = 6.
    """
    fname, n_swd, dk, A, nsteps, dt, omega, k1 = airy_deep_water_swd
    monkeypatch.setenv("SWD_NUM_H2_STEPS", "5")
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        nlayers = swd["nlayers"]
        assert nlayers == 6, f"Expected nlayers=6 (nstep+1=5+1), got {nlayers}"
    finally:
        swd.close()
    monkeypatch.delenv("SWD_NUM_H2_STEPS", raising=False)


# ---------------------------------------------------------------------------
# Test 9: grad_phi components are finite and non-zero (smoke test)
# ---------------------------------------------------------------------------


def test_grad_phi_finite(airy_deep_water_swd):
    """All components of grad_phi must be finite and at least one non-zero."""
    fname, n_swd, dk, A, nsteps, dt, omega, k1 = airy_deep_water_swd
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        swd.update_time(dt)
        gp = swd.grad_phi(2.5, 0.0, -1.0)
        assert math.isfinite(gp.x), f"grad_phi.x not finite: {gp.x}"
        assert math.isfinite(gp.z), f"grad_phi.z not finite: {gp.z}"
        assert abs(gp.x) + abs(gp.z) > 1e-6, "grad_phi is all zeros"
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 10: nlayers = nstep+1 with exact sigma positions
# ---------------------------------------------------------------------------


def test_default_nlayers_is_nstep_plus_one(airy_deep_water_swd):
    """Default nlayers must equal nstep+1 = 21."""
    fname = airy_deep_water_swd[0]
    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        assert swd["nlayers"] == 21, f"Expected 21, got {swd['nlayers']}"
    finally:
        swd.close()


# ---------------------------------------------------------------------------
# Test 11: nonlinear kinematics validated against raschii Stokes wave
# ---------------------------------------------------------------------------


def test_nonlinear_kinematics_vs_raschii(tmp_path):
    """
    Use raschii to write an amp=2 SWD file for a moderately steep Stokes wave
    (5th order), then compare H2 kinematics at depth against raschii's analytic
    velocity field.
    """
    raschii = pytest.importorskip("raschii", reason="raschii >= 2.0.0 required for this test")

    depth = 20.0  # m
    height = 2.0  # wave height  H/d = 0.10  (moderate steepness)
    length = 40.0  # wavelength
    N_order = 5

    wave = raschii.StokesWave(height=height, depth=depth, length=length, N=N_order)

    fname = str(tmp_path / "stokes_amp2.swd")
    # raschii 2.0.0 supports amp=2 for shape=2
    wave.write_swd(fname, tmax=wave.period * 4, dt=wave.period / 20, amp=2)

    swd = SpectralWaveData(fname, 0.0, 0.0, 0.0, 0.0, rho=1025.0, impl=0)
    try:
        t_eval = wave.period * 0.5
        swd.update_time(t_eval)

        x_points = [0.0, length / 4, length / 2]
        z_points = [-1.0, -3.0, -depth * 0.5]

        for x_app in x_points:
            for z_app in z_points:
                # raschii uses z=0 at bottom, SWD uses z=0 at still water
                z_raschii = z_app + wave.depth

                gp = swd.grad_phi(x_app, 0.0, z_app)
                u_swd = float(gp.x)
                w_swd = float(gp.z)

                vel = wave.velocity(x=x_app, z=z_raschii, t=t_eval)
                u_raschii = float(vel[0])  # horizontal
                w_raschii = float(vel[1])  # vertical

                # H2 should recover analytic kinematics within some percentage error (allowing for
                # some numerical error due to linear interpolation, H2-op is not perfect etc)
                u_scale = max(abs(u_raschii), 1e-4)
                w_scale = max(abs(w_raschii), 1e-4)
                err_u = abs(u_swd - u_raschii) / u_scale
                err_w = abs(w_swd - w_raschii) / w_scale
                assert err_u < 0.02, (
                    f"u at x={x_app}, z={z_app}: SWD={u_swd:.4g}, "
                    f"raschii={u_raschii:.4g} => err={err_u:.4g}"
                )
                assert err_w < 0.02, (
                    f"w at x={x_app}, z={z_app}: SWD={w_swd:.4g}, "
                    f"raschii={w_raschii:.4g} => err={err_w:.4g}"
                )
    finally:
        swd.close()
