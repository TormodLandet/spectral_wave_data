# Installing the spectral_wave_data Python-package from source

This Python package provides the spectral_wave_data ocean wave model.

Programmed by: Jens B. Helmers & Odin Gramstad, DNV

## Prerequisites for compiling the Spectral Wave Data library

### Windows

To build on Windows you will need to install:

* Intel Fortran compiler "ifort" from a recent (free) Intel OneAPI
  installation bundle.
* CMake
* Microsoft Visual Studio
* Python 3.x with the "wheel" package installed

If you do not have access to Visual Studio you may install `ninja` or
another build tool and modify the CMake command to use that as the build
system instead of creating a Visual Studio Solution. We have not tested
this, but it will probably work.

### Linux

To build on Linux you will need to install:

* A Fortran compiler, either gfortran 9 or later or a recent version of 
  ifort from the (free) Intel OneAPI installation bundle.
* CMake
* Make
* Python 3.x with the "wheel" package installed

## Building the Fortran library

### Windows

Go to the `Cmake` directory and call the provided bat file

```bat
cd src\api\python\Cmake
cmd /C CMakeBuild_Win64.bat
```

You may need to adjust the CMakeBuild_Win64.bat file to reflect your actual version of Visual Studio.

Details for building with Visual Studio:

1) Open the generated spectral_wave_data.sln with Visual Studio
2) Make sure the 'Release' flag and actual Binary configuration (64 vs 32bit) is selected in the tool bar of Visual Studio.
3) Build -> Build Solution
4) Check that the solution build without error messages and check that SpectralWaveData.dll is created in the Release sub-folder.
5) Put the `SpectralWaveData.dll` shared library file in the `Cmake/Build_Win64` sub-directory

### Linux

Go to the `Cmake` directory and call the provided bash file

```bash
cd src/api/python/Cmake
bash ./CMakeBuild_Linux.bash
```

## Building the Python package

You can build a Python package (`*.whl` file) by running pip (in this directory, the one containing `pyproject.toml`)

```bash
pip wheel .
```

You can now install the `spectral_wave_data` library from the generated wheel file
```bash
pip install spectral_wave_data-*.whl
```

You can also directly install the package without going via a `*.whl` file by running
```bash
pip install .
```

## Using the Python Wheel file on other computers

### Windows

Unless you have installed a recent Intel Fortran and Microsoft C++
compiler you may have to install some redistributable packages to be
downloaded from Intel and Microsoft:

* https://software.intel.com/en-us/articles/intel-compilers-redistributable-libraries-by-version
* https://www.microsoft.com/en-us/download/details.aspx?id=52685

### Linux

You may need to use the `auditwheel` tool to make the newly created Wheel
file into a "manylinux" Wheel which works across multiple Linux versions.

```bash
auditwheel repair spectral_wave_data*.whl
```
