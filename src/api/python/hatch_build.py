"""
This file is part of the hatchling-based build system
See pyproject.toml
"""

import sys
import sysconfig
from pathlib import Path

from hatchling.metadata.plugin.interface import MetadataHookInterface
from hatchling.builders.hooks.plugin.interface import BuildHookInterface

class CustomMetadataHook(MetadataHookInterface):
    """
    Our custom code to get the library version
    """

    def update(self, metadata: dict):
        py_dir = Path(self.root).resolve(strict=True)
        version_file = py_dir.parent.parent.parent / "version_definition.py"

        locals = {}
        exec(open(version_file).read(), globals(), locals)

        # Update the [project] table metadata
        metadata["version"] = locals["version_full"]


class CustomBuildHook(BuildHookInterface):
    """
    Our custom code to include the compiled library in the wheel
    and to set the appropriate wheel tag
    """

    def initialize(self, version: str, build_data: dict):
        ################################
        # Bundle the compiled library
        ################################

        py_dir = Path(self.root).resolve(strict=True)
        if sys.platform.startswith("linux"):
            lib = py_dir / "Cmake/Build_Linux/libSpectralWaveData.so"
        elif sys.platform.startswith("win"):
            lib = py_dir / "Cmake/Build_Win64/SpectralWaveData.dll"
        else:
            raise RuntimeError(f"Unsupported build platform: {sys.platform}")

        if not lib.exists():
            self.app.abort(
                f"The required file from the Fortran build folder is missing: {lib}"
            )

        build_data["force_include"][str(lib.resolve(strict=True))] = (
            f"spectral_wave_data/{lib.name}"
        )

        ################################
        # Setup the Wheel metadata
        ################################

        # We include a compiled library so the wheel is not pure Python
        build_data["pure_python"] = False

        # Create the tag for the Wheel file
        # See https://peps.python.org/pep-0425/
        # Should end up with
        # - py3-none-linux_x86_64
        # - py3-none-win_amd64
        # The "none" is because we don't have any ABI requirements since we
        # use ctypes and not the Python C API
        python_requirement = "py3"
        abi_requirement = "none"
        platform_requirement = (
            sysconfig.get_platform().replace("-", "_").replace(".", "_")
        )
        build_data["tag"] = (
            f"{python_requirement}-{abi_requirement}-{platform_requirement}"
        )
