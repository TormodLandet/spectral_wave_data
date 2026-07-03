# swd_fft.cmake — helper to wire the PocketFFT C++ backend to a CMake target.
#
# Include this file and call swd_enable_fft(your_target) for every library or
# executable that is built from source lists containing swd_pocketfft_c_api.cpp.
# The function is a no-op when USE_POCKETFFT is OFF.
#
# Usage (in any CMakeLists.txt that includes core_files.cmake):
#
#   include(<path>/swd_fft.cmake)    # or included transitively via core_files.cmake
#   add_library(mylib ...)
#   swd_enable_fft(mylib)

function(swd_enable_fft target)
    if(NOT USE_POCKETFFT)
        return()
    endif()

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
endfunction()
