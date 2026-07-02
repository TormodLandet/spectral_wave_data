# Build SpectralWaveData library for Python package

# CMake version must be at least 3.15 for the CMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded to work
cmake_minimum_required(VERSION 3.15)

set(LIB SpectralWaveData)

project(${LIB} C CXX Fortran)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(WIN32)
    set_property (DIRECTORY PROPERTY VS_STARTUP_PROJECT ${LIB})
    # To static link Intel libraries into SWD (not yet recognized by VS, to be set manually in VS. see Readme file)
    set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS_RELEASE} /libs:static")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_BUILD_TYPE Release)

set(DIR_SRC_API_F ../../fortran)
set(DIR_SRC_API_C ../../c)
set(DIR_SRC_API_PY .)

include(${DIR_SRC_API_F}/Cmake/core_files.cmake)

# Core files plus kind values from the Python directory
set(SRC_API_F
    ${SRC_CORE}
    ${DIR_SRC_API_PY}/kind_values.f90
)

add_library(${LIB} SHARED ${SRC_API_F})
target_include_directories(${LIB} INTERFACE ${DIR_SRC_API_C})

# Setup for the vendored PocketFFT library.
if(USE_POCKETFFT)
    target_include_directories(${LIB} PRIVATE ${CMAKE_SOURCE_DIR}/../../../thirdparty/pocketfft)
    
    # Disable the pocketfft thread pool
    target_compile_definitions(${LIB} PRIVATE POCKETFFT_NO_MULTITHREADING)
    
    # Statically embed libstdc++ and libgcc so the .so has no new runtime dependencies even if
    # we now use C++ features that may require linking to routines in the C++ standard library.
    if(UNIX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_link_options(${LIB} PRIVATE -static-libstdc++ -static-libgcc)
    endif()
endif()
