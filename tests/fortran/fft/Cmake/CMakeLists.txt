# Build and run Fortran FFT unit tests for SpectralWaveData.
# Mirrors the structure of tests/fortran/temporal_interpolation/Cmake/.

cmake_minimum_required(VERSION 3.15)

set(TEST_NAME test_swd_fft)
project(${TEST_NAME} C CXX Fortran)

set(DIR_SRC_API_F ${CMAKE_CURRENT_SOURCE_DIR}/../../../src/api/fortran)
set(DIR_TESTS_FFT ${CMAKE_CURRENT_SOURCE_DIR}/..)

include(${DIR_SRC_API_F}/Cmake/core_files.cmake)

# Library sources + kind_values for the Fortran API
set(SRC_LIB
    ${SRC_CORE}
    ${DIR_SRC_API_F}/kind_values.f90
)

add_library(swd_fft_test_lib STATIC ${SRC_LIB})
swd_enable_fft(swd_fft_test_lib)

add_executable(${TEST_NAME} ${DIR_TESTS_FFT}/fft/test_swd_fft.f90)
target_link_libraries(${TEST_NAME} PRIVATE swd_fft_test_lib)

enable_testing()
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
