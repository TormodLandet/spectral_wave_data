"""
Analytic reference class and binary SWD writer for shape class 7.

Shape 7: long-crested waves in a wave-following sigma-coordinate.
The velocity potential is stored on nlayers sigma-layers deforming with the
free surface.  No time derivatives are stored in the file; the Fortran
reader reconstructs them numerically from the four-step window.

The analytic kinematics implemented here mirror the Fortran implementation
exactly, so comparisons between this reference and the compiled library
validate correctness.
"""

import math
from struct import pack

import numpy as np

from test_utils import run_swd_meta_check


# ---------------------------------------------------------------------------
# Tfun helpers (evaluate SymPy polynomial Tfun objects numerically)
# ---------------------------------------------------------------------------

def _tfun_eval(tfun, t):
    """Evaluate a Tfun object at time t (returns Python complex)."""
    result = complex(0.0)
    t_pow = 1.0
    for c in tfun.c:
        result += complex(c) * t_pow
        t_pow *= t
    return result


def _tfun_eval_der(tfun, t):
    """Evaluate the time derivative of a Tfun object at time t."""
    result = complex(0.0)
    t_pow = 1.0
    for i in range(1, len(tfun.c)):
        result += i * complex(tfun.c[i]) * t_pow
        t_pow *= t
    return result


# ---------------------------------------------------------------------------
# Reference class
# ---------------------------------------------------------------------------

class Shape7:
    """Analytic reference for shape-7 sigma-coordinate wave kinematics.

    Parameters
    ----------
    dk : float
        Wave-number spacing.
    n : int
        Highest spectral index (n+1 components, j=0..n).
    d : float
        Water depth (> 0) or -1 for infinite depth.
    zref : float
        z-position of the bottom sigma-layer (must be <= 0,
        and d+zref > 0 for finite depth).
    nlayers : int
        Number of sigma layers (>= 2).
    sig : sequence of float, length nlayers
        Sigma positions: sig[0]=0, sig[-1]=1, strictly increasing.
    cfuns_layers : list of nlayers lists, each of length n+1 (Tfun objects)
        Spectral amplitudes per sigma-layer. cfuns_layers[m][j] is the
        Tfun for layer index m (0=bottom) and spectral component j.
    hfuns : list of n+1 Tfun objects
        Spectral amplitudes for surface elevation h_j(t).
    sys : corsys.CorSys
        Application-to-SWD coordinate transform (x0, y0, t0, beta).
    """

    def __init__(self, dk, n, d, zref, nlayers, sig, cfuns_layers, hfuns, sys):
        assert len(cfuns_layers) == nlayers, "cfuns_layers must have nlayers entries"
        for m in range(nlayers):
            assert len(cfuns_layers[m]) == n + 1
        assert len(hfuns) == n + 1
        assert abs(sig[0]) < 1.0e-9, "sig[0] must be 0"
        assert abs(sig[-1] - 1.0) < 1.0e-9, "sig[-1] must be 1"

        self.dk = dk
        self.n = n
        self.d = d
        self.zref = zref
        self.nlayers = nlayers
        self.sig = np.asarray(sig, dtype=float)
        self.cfuns = cfuns_layers
        self.hfuns = hfuns
        self.sys = sys

        self.cbeta = math.cos(sys.beta * math.pi / 180.0)
        self.sbeta = math.sin(sys.beta * math.pi / 180.0)

        if d > 0.0:
            self.d_eff = d + zref
            assert self.d_eff > 0.0, f"d_eff = d + zref = {self.d_eff} must be > 0"
        else:
            self.d_eff = None  # infinite depth

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _app2swd(self, x_app, y_app, t_app):
        xswd = self.sys.x0 + x_app * self.cbeta + y_app * self.sbeta
        tswd = self.sys.t0 + t_app
        return xswd, tswd

    def _Xj(self, xswd, j):
        """Spatial basis function  X_j = exp(-i k_j x_swd)."""
        return math.cos(j * self.dk * xswd) - 1j * math.sin(j * self.dk * xswd)

    def _Zfun_below(self, kj, z_prime):
        """Vertical basis function and z-derivative for below-zref extrapolation.

        z_prime = z - zref  (<= 0 when below zref).
        Matches the Fortran depth_Zfun subroutine.
        """
        if self.d_eff is None:
            # Infinite depth: Z_j = exp(k_j z')
            ez = math.exp(kj * z_prime)
            return ez, kj * ez
        else:
            # Finite depth: Z_j = cosh(kj*(z'+d_eff)) / cosh(kj*d_eff)
            sh = math.cosh(kj * (z_prime + self.d_eff))
            ch = math.cosh(kj * self.d_eff)
            Zfun = sh / ch
            Zfun_z = kj * math.sinh(kj * (z_prime + self.d_eff)) / ch
            return Zfun, Zfun_z

    def _eval_elev(self, xswd, tswd):
        """Return (zeta, zeta_x, zeta_t) at (xswd, tswd)."""
        zeta = 0.0
        zeta_x = 0.0
        zeta_t = 0.0
        for j in range(self.n + 1):
            kj = j * self.dk
            hj = _tfun_eval(self.hfuns[j], tswd)
            htj = _tfun_eval_der(self.hfuns[j], tswd)
            Xj = self._Xj(xswd, j)
            zeta += (hj * Xj).real
            if j > 0:
                # d/dx_swd Re(h_j X_j) = Re(h_j * (-ik_j) * X_j) = k_j Im(h_j X_j)
                zeta_x += kj * (hj * Xj).imag
            zeta_t += (htj * Xj).real
        return zeta, zeta_x, zeta_t

    def _eval_phi_layers(self, xswd, tswd):
        """Return (phi_m, u_m, Phi_t_m) numpy arrays of length nlayers.

        phi_m  : potential on each layer
        u_m    : d(phi_m)/d(x_swd) on each layer
        Phi_t_m: d(phi_m)/dt|sigma on each layer (time deriv at fixed sigma)
        """
        phi_m = np.zeros(self.nlayers)
        u_m = np.zeros(self.nlayers)
        Phi_t_m = np.zeros(self.nlayers)
        for j in range(self.n + 1):
            kj = j * self.dk
            Xj = self._Xj(xswd, j)
            for m in range(self.nlayers):
                cjm = _tfun_eval(self.cfuns[m][j], tswd)
                ctjm = _tfun_eval_der(self.cfuns[m][j], tswd)
                cx = cjm * Xj
                phi_m[m] += cx.real
                if j > 0:
                    u_m[m] += kj * cx.imag
                Phi_t_m[m] += (ctjm * Xj).real
        return phi_m, u_m, Phi_t_m

    def _interp_sigma(self, sigma_eval, phi_m, u_m, Phi_t_m=None):
        """Linear interpolation in sigma.  sigma_eval must already be clamped to [0,1].

        Returns (phi, u, phi_sigma, Phi_t_sigma).
        phi_sigma  = d(phi)/d(sigma) in the bracketing interval.
        """
        # Find bracket [m1, m2] containing sigma_eval
        m1 = self.nlayers - 2  # fallback: last interval
        for m in range(self.nlayers - 1):
            if sigma_eval <= self.sig[m + 1]:
                m1 = m
                break
        m2 = m1 + 1

        dsig = self.sig[m2] - self.sig[m1]
        s = (sigma_eval - self.sig[m1]) / dsig

        phi = float(phi_m[m1] + s * (phi_m[m2] - phi_m[m1]))
        u = float(u_m[m1] + s * (u_m[m2] - u_m[m1]))
        phi_sigma = float((phi_m[m2] - phi_m[m1]) / dsig)
        Phi_t_sigma = 0.0
        if Phi_t_m is not None:
            Phi_t_sigma = float(Phi_t_m[m1] + s * (Phi_t_m[m2] - Phi_t_m[m1]))
        return phi, u, phi_sigma, Phi_t_sigma

    # ------------------------------------------------------------------
    # Public kinematic methods (matching the Fortran API contract)
    # ------------------------------------------------------------------

    def phi(self, x_app, y_app, z_app, t_app):
        xswd, tswd = self._app2swd(x_app, y_app, t_app)
        if z_app >= self.zref:
            zeta, zeta_x, zeta_t = self._eval_elev(xswd, tswd)
            H = zeta - self.zref
            sigma = max(0.0, min(1.0, (z_app - self.zref) / H))
            phi_m, u_m, Phi_t_m = self._eval_phi_layers(xswd, tswd)
            phi, u, phi_sigma, _ = self._interp_sigma(sigma, phi_m, u_m)
            return phi
        else:
            z_prime = z_app - self.zref
            result = _tfun_eval(self.cfuns[0][0], tswd).real   # j=0, Z_0=1
            for j in range(1, self.n + 1):
                kj = j * self.dk
                Xj = self._Xj(xswd, j)
                cjm = _tfun_eval(self.cfuns[0][j], tswd)
                Zfun, _ = self._Zfun_below(kj, z_prime)
                result += (cjm * Xj).real * Zfun
            return result

    def phi_t(self, x_app, y_app, z_app, t_app):
        """Euler time derivative of phi at fixed physical z (chain-rule corrected)."""
        xswd, tswd = self._app2swd(x_app, y_app, t_app)
        if z_app >= self.zref:
            zeta, zeta_x, zeta_t = self._eval_elev(xswd, tswd)
            H = zeta - self.zref
            sigma = max(0.0, min(1.0, (z_app - self.zref) / H))
            phi_m, u_m, Phi_t_m = self._eval_phi_layers(xswd, tswd)
            phi, u, phi_sigma, Phi_t_sigma = self._interp_sigma(sigma, phi_m, u_m, Phi_t_m)
            # phi_t|z = Phi_t|sigma - sigma * zeta_t / H * phi_sigma
            return Phi_t_sigma - sigma * zeta_t / H * phi_sigma
        else:
            z_prime = z_app - self.zref
            result = _tfun_eval_der(self.cfuns[0][0], tswd).real
            for j in range(1, self.n + 1):
                kj = j * self.dk
                Xj = self._Xj(xswd, j)
                ctjm = _tfun_eval_der(self.cfuns[0][j], tswd)
                Zfun, _ = self._Zfun_below(kj, z_prime)
                result += (ctjm * Xj).real * Zfun
            return result

    def grad_phi(self, x_app, y_app, z_app, t_app):
        """Particle velocity (gradient of phi)."""
        xswd, tswd = self._app2swd(x_app, y_app, t_app)
        if z_app >= self.zref:
            zeta, zeta_x, zeta_t = self._eval_elev(xswd, tswd)
            H = zeta - self.zref
            sigma = max(0.0, min(1.0, (z_app - self.zref) / H))
            phi_m, u_m, Phi_t_m = self._eval_phi_layers(xswd, tswd)
            phi, u_swd, phi_sigma, _ = self._interp_sigma(sigma, phi_m, u_m)
            # Sigma chain-rule: d(phi)/d(x) = u_swd - sigma*zeta_x/H * phi_sigma
            phi_xswd = u_swd - sigma * zeta_x / H * phi_sigma
            phi_z = phi_sigma / H
            return {'x': phi_xswd * self.cbeta,
                    'y': phi_xswd * self.sbeta,
                    'z': phi_z}
        else:
            z_prime = z_app - self.zref
            phi_xswd = 0.0
            phi_z = 0.0
            for j in range(1, self.n + 1):
                kj = j * self.dk
                Xj = self._Xj(xswd, j)
                cjm = _tfun_eval(self.cfuns[0][j], tswd)
                Zfun, Zfun_z = self._Zfun_below(kj, z_prime)
                cx = cjm * Xj
                phi_xswd += kj * cx.imag * Zfun
                phi_z += cx.real * Zfun_z
            return {'x': phi_xswd * self.cbeta,
                    'y': phi_xswd * self.sbeta,
                    'z': phi_z}

    def pressure(self, x_app, y_app, z_app, t_app, rho, grav):
        """Fully nonlinear Bernoulli pressure."""
        phi_t_val = self.phi_t(x_app, y_app, z_app, t_app)
        vel = self.grad_phi(x_app, y_app, z_app, t_app)
        kin_e = 0.5 * (vel['x']**2 + vel['y']**2 + vel['z']**2)
        return -rho * phi_t_val - rho * grav * z_app - rho * kin_e

    def elev(self, x_app, y_app, t_app):
        xswd, tswd = self._app2swd(x_app, y_app, t_app)
        zeta, _, _ = self._eval_elev(xswd, tswd)
        return zeta

    def elev_t(self, x_app, y_app, t_app):
        xswd, tswd = self._app2swd(x_app, y_app, t_app)
        _, _, zeta_t = self._eval_elev(xswd, tswd)
        return zeta_t

    def grad_elev(self, x_app, y_app, t_app):
        xswd, tswd = self._app2swd(x_app, y_app, t_app)
        elev_x_swd = 0.0
        for j in range(1, self.n + 1):
            kj = j * self.dk
            hj = _tfun_eval(self.hfuns[j], tswd)
            Xj = self._Xj(xswd, j)
            elev_x_swd += kj * (hj * Xj).imag
        return {'x': elev_x_swd * self.cbeta,
                'y': elev_x_swd * self.sbeta,
                'z': 0.0}

    def grad_elev_2nd(self, x_app, y_app, t_app):
        xswd, tswd = self._app2swd(x_app, y_app, t_app)
        elev_xx_swd = 0.0
        for j in range(1, self.n + 1):
            kj = j * self.dk
            hj = _tfun_eval(self.hfuns[j], tswd)
            Xj = self._Xj(xswd, j)
            elev_xx_swd -= kj**2 * (hj * Xj).real
        return {'xx': elev_xx_swd * self.cbeta**2,
                'xy': elev_xx_swd * self.sbeta * self.cbeta,
                'yy': elev_xx_swd * self.sbeta**2}

    # ------------------------------------------------------------------
    # SWD binary writer
    # ------------------------------------------------------------------

    def write_swd(self, file_swd, dt, nsteps, too_short_file=False):
        """Write a shape-7 SWD binary file consumable by the Fortran/Python API.

        Per time step the file contains:
          h[0..n]              (elevation amplitudes, complex float32)
          c_m[0..n]  m=0..nlayers-1  (potential amplitudes, no derivatives)
        """
        with open(file_swd, 'wb') as out:
            # ---- Common SWD header ----
            out.write(pack('<f', 37.0221))   # magic
            out.write(pack('<i', 100))        # fmt
            out.write(pack('<i', 7))          # shp = 7
            out.write(pack('<i', 1))          # amp = 1  (only supported value)
            out.write(pack('<30s', b'shape7_sympy_test'.ljust(30)))
            out.write(pack('<20s', b'yyyy:mm:dd hh:ss'.ljust(20)))
            cid = b"{'shape':7,'test':True}"
            nid = len(cid)
            out.write(pack('<i', nid))
            out.write(pack(f'<{nid}s', cid))
            out.write(pack('<f', 9.81))       # grav
            out.write(pack('<f', 1.0))        # lscale
            out.write(pack('<i', 0))          # nstrip
            out.write(pack('<i', nsteps))
            out.write(pack('<f', dt))
            out.write(pack('<i', -1))         # order (-1 = fully nonlinear)
            # ---- Shape-7 header ----
            out.write(pack('<i', self.n))
            out.write(pack('<f', self.dk))
            out.write(pack('<f', self.d))
            out.write(pack('<f', self.zref))
            out.write(pack('<i', self.nlayers))
            for m in range(self.nlayers):
                out.write(pack('<f', float(self.sig[m])))
            # ---- Time steps ----
            nout = nsteps // 2 if too_short_file else nsteps
            for i in range(nout):
                t_swd = i * dt
                # Elevation h[0..n]
                for j in range(self.n + 1):
                    hj = _tfun_eval(self.hfuns[j], t_swd)
                    out.write(pack('<f', float(hj.real)))
                    out.write(pack('<f', float(hj.imag)))
                # Potential per layer c_m[0..n]
                for m in range(self.nlayers):
                    for j in range(self.n + 1):
                        cjm = _tfun_eval(self.cfuns[m][j], t_swd)
                        out.write(pack('<f', float(cjm.real)))
                        out.write(pack('<f', float(cjm.imag)))

    def check_swd_meta(self, file_swd, n):
        # swd_meta does not expose shape-7 specific fields (n, dk, nlayers, …)
        # in its output, so we only verify shp=7 is reported correctly.
        run_swd_meta_check(file_swd, shp=7)
