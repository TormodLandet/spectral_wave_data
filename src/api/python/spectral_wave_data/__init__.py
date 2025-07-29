# During development, the version file may not exist yet
# It is generated while building the wheel/source distribution
version_short: str = "0.0.0"
version_full: str = "0.0.0"
try:
    from .version import version_short, version_full
except ImportError:
    pass

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
