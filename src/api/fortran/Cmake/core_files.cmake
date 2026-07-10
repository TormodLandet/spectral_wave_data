# All the Fortran implementation source files except kind_values.f90
set(SRC_CORE
  ${DIR_SRC_API_F}/swd_fft_lib.f90
  ${DIR_SRC_API_F}/hosm_h2_operator.f90
  ${DIR_SRC_API_F}/multilayer_long_crested.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_1_or_2_impl_7.f90
  ${DIR_SRC_API_F}/open_swd_file.F90
  ${DIR_SRC_API_F}/spectral_interpolation.f90
  ${DIR_SRC_API_F}/spectral_wave_data.f90
  ${DIR_SRC_API_F}/spectral_wave_data_allocate.f90
  ${DIR_SRC_API_F}/spectral_wave_data_c.f90
  ${DIR_SRC_API_F}/spectral_wave_data_error.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_1_impl_1.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_2_impl_1.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_3_impl_1.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_4_impl_1.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_4_impl_2.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_5_impl_1.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_6_impl_1.f90
  ${DIR_SRC_API_F}/spectral_wave_data_shape_7_impl_1.f90
  ${DIR_SRC_API_F}/swd_write_shape_1_or_2.f90
  ${DIR_SRC_API_F}/swd_write_shape_3.f90
  ${DIR_SRC_API_F}/swd_write_shape_4_or_5.f90
  ${DIR_SRC_API_F}/swd_write_shape_6.f90
  ${DIR_SRC_API_F}/swd_write_shape_7.f90
  ${DIR_SRC_API_F}/swd_version.f90
)

# You need to define DIR_THIRDPARTY before including swd_fft.cmake
# so that it can locate the vendored PocketFFT sources.
set(DIR_THIRDPARTY ${DIR_SRC_API_F}/../../thirdparty)

# Include the shared FFT helper functions:
# * swd_fft_define_sources(source_list_variable)
# * swd_fft_enable_for_target(target)
# Call the last one on every library/executable using SWD's FFT lib.
include(${CMAKE_CURRENT_LIST_DIR}/swd_fft.cmake)

# Add the FFT backend sources to the core source list.
swd_fft_define_sources(SRC_FFT)
list(APPEND SRC_CORE ${SRC_FFT})

# Bundle the Intel compiler libraries when compiling shared libraries
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (UNIX)
        message("SWD: enabling -static-intel for .so files")
        set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -static-intel")
    endif()
    if (WIN32)
        message("SWD: enabling /static for .lib files")
        set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} /libs:static /threads")
    endif()
endif()


# Bundle the Intel compiler libraries when compiling shared libraries on Linux
#if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND UNIX)
#    message("SWD: enabling -static-intel for .so files")
#    set(CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_Fortran_FLAGS} -static-intel")
#endif()