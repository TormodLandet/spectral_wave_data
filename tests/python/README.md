# SWD Python unit tests


## Example using uv

First create a virtual environment for the tests and install the requirements.

These example commands are for Linux, similar commands will work on other operating systems like
Windows:

```bash
cd tests/python
# Create a virtual environment
uv venv .venv
# Install external dependencies
uv pip install --python=.venv/bin/python -r requirements_tests.txt
# Install the freshly compiled and packaged *.whl file
uv pip install --python=.venv/bin/python --reinstall [WHEEL_FILE]
# Run the quick unit tests (takes about 1 minute)
SWD_TEST_TYPE=quick .venv/bin/python -m pytest -v
```

Here `[WHEEL_FILE]` is likely `../../src/api/python/dist/spectral_wave_data-*.whl` if you built the
wheel using `uv build --wheel` in the `src/api/python` directory. Remember to compile the library
(`libSpectralWaveData.so` on Linux) before building the wheel. The wheel-builder just *copies* this
file which must *already* exist from a previous compilation using CMake etc,
**it does not compile the Fortran code as a part of the building of the wheel file**!


## Example using helper scripts

All these scripts should be executed in a terminal window having access to python 3

To install all required python packages run 'install_requirements_tests.bat' in the same environment.

Also install the actual python package spectral_wave_data to test. You need to install the wheel you
have just built.

Run e.g. run_all_tests.bat
