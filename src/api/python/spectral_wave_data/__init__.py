# During development, the version file may not exist yet
# It is generated while building the wheel/source distribution
try:
    from .version import version_short, version_full
except ImportError:
    version_short = "0.0"
    version_full = "0.0.0 (local development version)"

from .spectral_wave_data import (
    SpectralWaveData,
    SwdError,
    SwdFileCantOpenError,
    SwdFileBinaryError,
    SwdFileDataError,
    SwdInputValueError,
    SwdAllocateError,
    SwdIsClosedError,
)
