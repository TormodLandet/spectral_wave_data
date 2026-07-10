###############################################################################
# SWD FFT CMake helper utilities

###############################################################################
# FFT backend selection
#
# The developer compiling SWD must choose one of the FFT backends:
#
# POCKETFFT: vendored, zero extra dependencies, should be used for the Python'
#            wheels we are going to distribute.
#
#            POCKETFFT is the **DEFAULT** option! Use it unless you have a good
#            reason to use FFTW3 instead (or MKL's FFTW3 interface).
#
# FFTW     : requires system FFTW3 (dnf install fftw-devel / apt install
#            libfftw3-dev, OR SIMILAR); selected with -DSWD_FFT_BACKEND=FFTW.
#
#            IF you link to normal FFTW then the resulting SWD library will
#            have license **GPL**, not MIT! We distribute the wheel as MIT,
#            hence that CAN NOT linkt to fftw.
#
#            The Intel Math Kernel Library, MKL, also has a fftw3.f03 and is
#            or at least should be compatible with FFTW3. If you have MKL
#            installed you can possibly use that instead of FFTW3 (untested).
#
# Exactly one backend file is compiled; both expose the identical
# `module swd_fft_backend` interface so the rest of the build is unaffected.
#
set(SWD_FFT_BACKEND "POCKETFFT" CACHE STRING
    "FFT backend to use: POCKETFFT (default, vendored) or FFTW (fftw3.f03)")
set_property(CACHE SWD_FFT_BACKEND PROPERTY STRINGS POCKETFFT FFTW)


# #############################################################################
# Setup source files for the chosen FFT backend
#
# Use this function to get the list of source files for the chosen FFT backend.
# 
# Usage (in any CMakeLists.txt that includes core_files.cmake):
#
#     swd_fft_define_sources(SRC_FFT)
#     list(APPEND SRC_CORE ${SRC_FFT})
#
function(swd_fft_define_sources source_list_variable)
    if(SWD_FFT_BACKEND STREQUAL "POCKETFFT")
        set(SRC_FFT
            ${DIR_SRC_API_F}/swd_fft_backend_pocketfft.f90
            ${DIR_THIRDPARTY}/pocketfft/swd_pocketfft_c_api.cpp
        )
    elseif(SWD_FFT_BACKEND STREQUAL "FFTW")
        set(SRC_FFT
            ${DIR_SRC_API_F}/swd_fft_backend_fftw.f90
        )
    else()
        message(FATAL_ERROR "SWD_FFT_BACKEND must be POCKETFFT or FFTW (got: ${SWD_FFT_BACKEND})")
    endif()
    # Set the output variable to the list of FFT source files with the name specified by the caller.
    set(${source_list_variable} "${SRC_FFT}" PARENT_SCOPE)
    message(STATUS "SWD will use FFT backend: ${SWD_FFT_BACKEND}")
endfunction()


# #############################################################################
# Setup linking for FFT libraries
# 
# After adding the FFT sources (above) you must call swd_fft_enable_for_target(your_target)
# for every library or executable that is built from SRC_CORE (i.e. every target
# that includes source files swd_fft_lib.f90 and swd_fft_backend_*.f90).
#
# Usage (in any CMakeLists.txt that includes core_files.cmake):
#
#   include(<path>/swd_fft.cmake)    # or included transitively via core_files.cmake
#   add_library(mylib ...)
#   swd_fft_enable_for_target(mylib)
#
function(swd_fft_enable_for_target target)
    if(SWD_FFT_BACKEND STREQUAL "POCKETFFT")
        # pocketfft is a single-header C++17 library; add its directory to the
        # include path of the target.
        target_include_directories(${target} PRIVATE
            ${DIR_THIRDPARTY}/pocketfft)

        # Disable the pocketfft thread pool (SWD objects are not thread-safe anyway,
        # and the thread pool adds a libpthread dependency we do not want).
        target_compile_definitions(${target} PRIVATE
            POCKETFFT_NO_MULTITHREADING)

        # Require C++17 for the .cpp translation unit.
        target_compile_features(${target} PRIVATE cxx_std_17)

        # Linux / macOS (GCC or Clang): statically embed the C++ and GCC runtime
        # so the shared library loaded by Python via ctypes has no new runtime deps.
        if(UNIX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
            target_link_options(${target} PRIVATE
                -static-libstdc++
                -static-libgcc)
        endif()

        # Windows (MSVC): statically link the MSVC runtime (/MT) to avoid a
        # dependency on MSVCP140.dll in the DLL loaded by Python.
        if(WIN32 AND MSVC)
            set_property(TARGET ${target} PROPERTY
                MSVC_RUNTIME_LIBRARY "MultiThreaded")
        endif()

    elseif(SWD_FFT_BACKEND STREQUAL "FFTW")
        # Locate FFTW3 (double precision). Either Intel oneMKL or system FFTW3
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel"  # Intel Fortran Compiler Classic, ifort
            OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")  # Intel Fortran Compiler (ifx)
            # Use fftw3-compatible implementation from Intel MKL
            # Tested and working on Linux with FC=ifx from oneAPI 2024.1
            if(NOT DEFINED ENV{MKLROOT})
                # MKLROOT is not set, try to guess it from the compiler path.
                get_filename_component(COMP_PATH ${CMAKE_Fortran_COMPILER} DIRECTORY)
                set(MKLROOT "${COMP_PATH}/../mkl")
            else()
                # Use MKLROOT from the environment variable.
                set(MKLROOT $ENV{MKLROOT})
            endif()
            message(STATUS "SWD will use MKLROOT=${MKLROOT}")

            # Knowing MKLROOT we can guess the location of fftw3.f03
            set(FFTW3_F03_INCLUDE_DIR "${MKLROOT}/include/fftw")
            
            # We let CMake find the MKL libraries to link to using the built-in FindBLAS
            set(BLAS_VENDOR "Intel10_64_seq")  # Use the sequential Intel MKL library on x86_64. 
            find_package(BLAS REQUIRED)
            target_link_libraries(${target} PRIVATE BLAS::BLAS)
            get_target_property(FFTW3_LIB_TO_USE BLAS::BLAS INTERFACE_LINK_LIBRARIES)  # For message() only

        else()
            # Use "normal" FFTW3 from system libraries (GPL license!)
            # On RHEL/CentOS: dnf install fftw-devel, on Debian/Ubuntu: apt install libfftw3-dev
            find_package(FFTW3 QUIET)
            if(FFTW3_FOUND)
                target_link_libraries(${target} PRIVATE FFTW3::fftw3)
                get_target_property(FFTW3_LIB_TO_USE FFTW3::fftw3 LOCATION)  # For message() only
            else()
                # Fall back to manual library search if find_package is unavailable
                find_library(FFTW3_LIB fftw3 REQUIRED
                    DOC "Path to the FFTW3 double-precision shared/static library")
                target_link_libraries(${target} PRIVATE ${FFTW3_LIB})
                set(FFTW3_LIB_TO_USE ${FFTW3_LIB})  # For message() only
            endif()
        
            # The Fortran backend uses `include 'fftw3.f03'`, so the directory that
            # contains fftw3.f03 must be on the compiler include path.  find_package
            # does not always set this (e.g. when only find_library succeeded), so
            # search for the header explicitly.  Common locations: /usr/include.
            find_path(FFTW3_F03_INCLUDE_DIR fftw3.f03
                HINTS ${FFTW3_INCLUDE_DIRS}
                PATHS /usr/include /usr/local/include
                DOC "Directory containing the FFTW3 Fortran include file fftw3.f03")
            if(NOT FFTW3_F03_INCLUDE_DIR)
                message(FATAL_ERROR
                    "SWD_FFT_BACKEND=FFTW: could not find fftw3.f03. Install the FFTW "
                    "development package (dnf install fftw-devel / apt install "
                    "libfftw3-dev) or set -DFFTW3_F03_INCLUDE_DIR=<dir>.")
            endif()
        endif()
        target_include_directories(${target} PRIVATE ${FFTW3_F03_INCLUDE_DIR})
        message(STATUS "SWD will use fftw3.f03 from ${FFTW3_F03_INCLUDE_DIR}")
        message(STATUS "SWD will use fftw3 library  ${FFTW3_LIB_TO_USE}")
    endif()
endfunction()
