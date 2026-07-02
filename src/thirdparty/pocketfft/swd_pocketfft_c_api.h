/*
 * Pure C API for the pocketfft C++ header-only library.
 *
 * This header is included by C and Fortran (via iso_c_binding) callers.
 * The implementation is in pocketfft_c_api.cpp (C++17).
 *
 * All functions return 0 on success, -1 on error.
 * The complex data format is interleaved (re0, im0, re1, im1, ...) doubles,
 * matching both C's double _Complex and Fortran's complex(c_double) layout.
 */
#ifndef POCKETFFT_C_API_H
#define POCKETFFT_C_API_H

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque handle to a 1-D real FFT plan. */
typedef struct SwdRfftPlan1D* swd_rfft_plan_t;

/*
 * Create a plan for 1-D real FFTs of length n.
 * Writes the new plan into *plan_out.  Returns 0 on success, -1 on failure.
 */
int swd_rfft_plan_create(int n, swd_rfft_plan_t *plan_out);

/* Destroy a plan previously created by swd_rfft_plan_create. */
void swd_rfft_plan_destroy(swd_rfft_plan_t plan);

/*
 * Forward real-to-complex 1-D FFT (matches numpy.fft.rfft scaling: scale=1).
 *   in  : real[n]         (const double*)
 *   out : complex[n/2+1]  (void* to interleaved double pairs)
 * Returns 0 on success, -1 on error.
 */
int swd_rfft_forward(swd_rfft_plan_t plan, const double *in, void *out);

/*
 * Backward complex-to-real 1-D FFT (matches numpy.fft.irfft: scale = 1/n).
 *   in  : complex[n/2+1]  (const void* to interleaved double pairs)
 *   out : real[n]          (double*)
 * Returns 0 on success, -1 on error.
 */
int swd_rfft_backward(swd_rfft_plan_t plan, const void *in, double *out);

#ifdef __cplusplus
}
#endif
#endif /* POCKETFFT_C_API_H */
