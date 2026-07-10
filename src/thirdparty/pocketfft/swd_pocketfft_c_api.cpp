/*
 * Thin C wrapper around the pocketfft C++ header-only library.
 *
 * Exposes a plain C API so Fortran code can call it via iso_c_binding
 * without any C++ in the Fortran sources.  All C++ exceptions are caught
 * here and converted to integer error codes; no C++ exception may escape
 * the extern "C" boundary.
 *
 * Build requirements: C++17  (pocketfft_hdronly.h uses std::optional etc.)
 * Thread safety: pocketfft caches plans in thread-local storage; distinct
 * SwdRfftPlan1D objects are independent and safe to use from different threads
 * simultaneously.  A single object must not be used concurrently.
 */
/* POCKETFFT_NO_MULTITHREADING is passed via -D on the compiler command line (CMake).
 * Do not redefine here to avoid the "macro redefined" warning. */
#include "pocketfft_hdronly.h"
#include "swd_pocketfft_c_api.h"

#include <complex>
#include <new>           /* std::bad_alloc */
#include <stdexcept>

/* -------------------------------------------------------------------------
 * Internal plan structure — just stores n so callers can access it.
 * pocketfft caches its own internal twiddle tables the first time a given
 * length is requested, so we do not need to store anything else here.
 * -------------------------------------------------------------------------*/
struct SwdRfftPlan1D {
    int n;
};

/* -------------------------------------------------------------------------
 * Internal 2-D plan structure — records nx and ny.  PocketFFT caches its own
 * twiddle tables internally by size, so nothing else is needed here.
 * -------------------------------------------------------------------------*/
struct SwdRfft2Plan {
    int nx;
    int ny;
};

/* -------------------------------------------------------------------------
 * Public C API
 * -------------------------------------------------------------------------*/
extern "C" {

int swd_rfft_plan_create(int n, swd_rfft_plan_t *plan_out)
{
    if (!plan_out || n <= 0) return -1;
    try {
        *plan_out = new SwdRfftPlan1D{n};
        return 0;
    } catch (...) {
        *plan_out = nullptr;
        return -1;
    }
}

void swd_rfft_plan_destroy(swd_rfft_plan_t plan)
{
    delete plan;   /* delete nullptr is a no-op per C++ standard */
}

int swd_rfft_forward(swd_rfft_plan_t plan, const double *in, void *out)
{
    if (!plan || !in || !out) return -1;
    try {
        const int n = plan->n;
        pocketfft::shape_t  shape_in {static_cast<size_t>(n)};
        pocketfft::stride_t stride_in{static_cast<ptrdiff_t>(sizeof(double))};
        pocketfft::stride_t stride_out{static_cast<ptrdiff_t>(sizeof(std::complex<double>))};

        pocketfft::r2c(shape_in, stride_in, stride_out,
                       /*axis=*/0, /*forward=*/true,
                       in,
                       reinterpret_cast<std::complex<double>*>(out),
                       /*fct=*/1.0,
                       /*nthreads=*/1);
        return 0;
    } catch (...) {
        return -1;
    }
}

int swd_rfft_backward(swd_rfft_plan_t plan, const void *in, double *out)
{
    if (!plan || !in || !out) return -1;
    try {
        const int n = plan->n;
        pocketfft::shape_t  shape_out{static_cast<size_t>(n)};
        pocketfft::stride_t stride_in {static_cast<ptrdiff_t>(sizeof(std::complex<double>))};
        pocketfft::stride_t stride_out{static_cast<ptrdiff_t>(sizeof(double))};

        pocketfft::c2r(shape_out, stride_in, stride_out,
                       /*axis=*/0, /*forward=*/false,
                       reinterpret_cast<const std::complex<double>*>(in),
                       out,
                       /*fct=*/1.0 / static_cast<double>(n),
                       /*nthreads=*/1);
        return 0;
    } catch (...) {
        return -1;
    }
}

int swd_rfft2_plan_create(int nx, int ny, swd_rfft2_plan_t *plan_out)
{
    if (!plan_out || nx <= 0 || ny <= 0) return -1;
    try {
        *plan_out = new SwdRfft2Plan{nx, ny};
        return 0;
    } catch (...) {
        *plan_out = nullptr;
        return -1;
    }
}

void swd_rfft2_plan_destroy(swd_rfft2_plan_t plan)
{
    delete plan;   /* delete nullptr is a no-op per C++ standard */
}

int swd_rfft2_c2r(swd_rfft2_plan_t plan, const void *in, double *out)
{
    if (!plan || !in || !out) return -1;
    try {
        using C = std::complex<double>;
        const int    nx   = plan->nx;
        const int    ny   = plan->ny;
        const size_t snx  = static_cast<size_t>(nx);
        const size_t sny  = static_cast<size_t>(ny);
        const size_t snxh = snx / 2 + 1;

        // Fortran column-major layout: axis 0 is x (contiguous, fastest varying),
        // axis 1 is y.  The half-spectrum has shape (nx/2+1, ny).
        //
        // PocketFFT c2r with multiple axes:
        //   "first carry out a c2c transform along all axes except the last one,
        //    then a c2r transform on the last axis in axes"
        // We want c2r on x (axis 0) and c2c on y (axis 1).
        // => axes = {1, 0}: c2c on y first, then c2r on x.
        // For ny = 1 the c2c on axis 1 is a trivial no-op (single element).
        pocketfft::shape_t  shape_out {snx,            sny};
        pocketfft::stride_t stride_in {static_cast<ptrdiff_t>(sizeof(C)),
                                       static_cast<ptrdiff_t>(sizeof(C) * snxh)};
        pocketfft::stride_t stride_out{static_cast<ptrdiff_t>(sizeof(double)),
                                       static_cast<ptrdiff_t>(sizeof(double) * snx)};
        pocketfft::shape_t  axes      {1, 0};

        pocketfft::c2r(shape_out, stride_in, stride_out, axes,
                       /*forward=*/false,
                       reinterpret_cast<const C *>(in), out,
                       /*fct=*/1.0,   // unnormalized: scale=1 (caller pre-scales coefficients)
                       /*nthreads=*/1);
        return 0;
    } catch (...) {
        return -1;
    }
}

} /* extern "C" */
