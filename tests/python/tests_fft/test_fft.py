import os
import random
from itertools import product
from pathlib import Path

import numpy as np
import pytest
from spectral_wave_data import SpectralWaveData

# Test type can be normal, quick, or all
# If SWD_TEST_TYPE is set to "quick", then only a small random sample are run
SWD_TEST_TYPE = os.environ.get("SWD_TEST_TYPE", "normal")

# How many test cases to run for the tests. Running all these tests takes 5-6 minute
# All test cases passed in July 2026 (150528 passed when using `-k fft` with pytest)
NUM_CASES_QUICK_TESTS = 5_000  # Set SWD_TEST_TYPE="quick"

TEST_FILES_DIR = Path(__file__).parent / "inputfiles"

# Define the HUGE matrix of test cases to run (150 528 tests!!!)
ALL_FUNCS = ["elev", "grad_phi"]
Z_INDEPENDENT_FUNCS = ["elev"]
ALL_SHP_IMPL = [(1, 1), (2, 1), (4, 1), (4, 2), (5, 1)]
ALL_X0_Y0_BETA = [(0.0, 0.0, 0.0), (-23.67, 56.4, -23.65)]
ALL_NX = [32, 33, 34, 35]
ALL_NY = [1, 32, 33, 34, 35]
ALL_NSUMX = [-1, 13, 14]
ALL_NSUMY = [-1, 13, 14]
ALL_DC_BIAS = [True, False]
ALL_NORDER = [-1, 0, 2, 4, 7]
ALL_NX_FFT = [-1, -2, 50, 51]
ALL_NY_FFT = [-1, -2, 50, 51]


def enumerate_all_fft_test_cases(quick: bool = False):
    """
    make a list with the combination of the parameters to test
    """
    case_params = []
    test_ids = []
    for shape, impl in ALL_SHP_IMPL:
        for x0, y0, beta in ALL_X0_Y0_BETA:
            for nx in ALL_NX:  # hosm nx
                for ny in ALL_NY:  # hosm ny
                    swd_name = f"shape{shape}_impl{impl}_nx{nx}_ny{ny}.swd"
                    swd_path = TEST_FILES_DIR / swd_name
                    if not swd_path.is_file():
                        # This swd doesn't exist ==> illegal combination of shp, impl, nx, ny
                        continue

                    for nsumx, nsumy, dc_bias, norder, nx_fft, ny_fft in product(
                        ALL_NSUMX, ALL_NSUMY, ALL_DC_BIAS, ALL_NORDER, ALL_NX_FFT, ALL_NY_FFT
                    ):
                        # skip values for nsumy and ny_fft not valid for 1D
                        if shape in (1, 2) and (nsumy > 1 or abs(ny_fft) != 1):
                            continue
                        for func in ALL_FUNCS:
                            for z in [0, -2.1, 2.1]:
                                # skip z-dependent for fields independent of z
                                if func in Z_INDEPENDENT_FUNCS and (norder != 0 or z != 0):
                                    continue
                                # skip norder for z <= 0
                                if z <= 0 and norder != 0:
                                    continue

                                test_ids.append(
                                    f"{shape=}, {impl=}, {x0=}, {y0=}, {beta=}, {nx=}, {ny=},"
                                    f" {nsumx=}, {nsumy=}, {dc_bias=}, {norder=}, {nx_fft=},"
                                    f" {ny_fft=}, {func=}, {z=}"
                                )
                                case_params.append(
                                    (
                                        shape,
                                        impl,
                                        x0,
                                        y0,
                                        beta,
                                        nx,
                                        ny,
                                        nsumx,
                                        nsumy,
                                        dc_bias,
                                        norder,
                                        nx_fft,
                                        ny_fft,
                                        func,
                                        z,
                                    )
                                )
    # How many of the test cases to run
    NUM_CASES = len(case_params)
    if SWD_TEST_TYPE == "quick":
        NUM_CASES = NUM_CASES_QUICK_TESTS

    if NUM_CASES < len(case_params):
        # shuffle the list and take a random sample of NUM_CASES test cases
        random.seed(42)
        combined = list(zip(test_ids, case_params))
        random.shuffle(combined)
        test_ids, case_params = zip(*combined)
        test_ids = list(test_ids)[:NUM_CASES]
        case_params = list(case_params)[:NUM_CASES]

    return case_params, test_ids


case_params, test_ids = enumerate_all_fft_test_cases()


@pytest.mark.parametrize("test_case", case_params, ids=test_ids)  # tests to run through
def test_fft(test_case):
    shape, impl, x0, y0, beta, nx, ny, nsumx, nsumy, dc_bias, norder, nx_fft, ny_fft, func, z = (
        test_case
    )
    swd_name = f"shape{shape}_impl{impl}_nx{nx}_ny{ny}.swd"
    swd_path = TEST_FILES_DIR / swd_name
    with SpectralWaveData(
        swd_path,
        x0=x0,
        y0=y0,
        beta=beta,
        nsumx=nsumx,
        nsumy=nsumy,
        norder=norder,
        dc_bias=dc_bias,
    ) as swd:
        swd.update_time(0.0)

        x_fft, y_fft = swd.xy_fft(nx_fft=nx_fft, ny_fft=ny_fft)
        assert x_fft.shape == y_fft.shape

        # test only nxp x nyp gridpoints
        nxp, nyp = 5, 5
        ixs = np.unique(np.linspace(0, x_fft.shape[0] - 1, nxp).astype(int))
        iys = np.unique(np.linspace(0, y_fft.shape[1] - 1, nyp).astype(int))

        if func in Z_INDEPENDENT_FUNCS:
            arr_fft = getattr(swd, f"{func}_fft")(nx_fft=nx_fft, ny_fft=ny_fft)
        else:
            arr_fft = getattr(swd, f"{func}_fft")(z, nx_fft=nx_fft, ny_fft=ny_fft)

        # scalar output
        if arr_fft.ndim == 2:
            res = np.zeros((ixs.size, iys.size))
            res_fft = np.zeros((ixs.size, iys.size))
            assert x_fft.shape == arr_fft.shape
            for i, ix in enumerate(ixs):
                for j, iy in enumerate(iys):
                    if func in Z_INDEPENDENT_FUNCS:
                        res[i, j] = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy])
                    else:
                        res[i, j] = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy], z)
                    res_fft[i, j] = arr_fft[ix, iy]
            assert np.allclose(res, res_fft)
        elif arr_fft.ndim == 3:  # 3-component output
            res = np.zeros((3, ixs.size, iys.size))
            res_fft = np.zeros((3, ixs.size, iys.size))
            assert (3, x_fft.shape[0], x_fft.shape[1]) == arr_fft.shape
            for i, ix in enumerate(ixs):
                for j, iy in enumerate(iys):
                    if func in Z_INDEPENDENT_FUNCS:
                        resxyz = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy])
                    else:
                        resxyz = getattr(swd, func)(x_fft[ix, iy], y_fft[ix, iy], z)
                    res[:, i, j] = np.asarray([resxyz.x, resxyz.y, resxyz.z])
                    res_fft[:, i, j] = arr_fft[:, ix, iy]
            print(np.max(np.abs(res[0, :, :] - res_fft[0, :, :])))
            print(np.max(np.abs(res[1, :, :] - res_fft[1, :, :])))
            print(np.max(np.abs(res[2, :, :] - res_fft[2, :, :])))
            assert np.allclose(res[0, :, :], res_fft[0, :, :])
            assert np.allclose(res[1, :, :], res_fft[1, :, :])
            assert np.allclose(res[2, :, :], res_fft[2, :, :])
