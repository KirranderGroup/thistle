"""Compiles the Fortran kernels when the package is built (``pip install .``).

The library is built with ``make`` (gfortran + OpenMP by default).  Intel oneAPI (load its
environment first, e.g. ``source /opt/intel/oneapi/setvars.sh`` or ``module load intel``):

    FC=ifx pip install .
    FC=ifx EXTRA=-xHost pip install .        # optimised for the build machine's CPU

Set THISTLE_SKIP_MAKE=1 to package an already built library as it is.
"""
import os
import shutil
import subprocess

from setuptools import Distribution, setup
from setuptools.command.build_py import build_py

HERE = os.path.dirname(os.path.abspath(__file__))


class BuildWithFortran(build_py):
    def run(self):
        if not os.environ.get("THISTLE_SKIP_MAKE"):
            if shutil.which("make") is None:
                raise RuntimeError("`make` is needed to compile the Fortran kernels")
            fc = os.environ.get("FC", "gfortran")
            if shutil.which(fc.split()[0]) is None:
                hint = ("load the Intel oneAPI environment first (`source /opt/intel/oneapi/setvars.sh` "
                        "or `module load intel`)" if os.path.basename(fc.split()[0]).split("-")[-1] in ("ifx", "ifort") else
                        "install gfortran (Linux: your package manager, e.g. `apt install gfortran`; "
                        "macOS: `brew install gcc`; clusters: `module load gcc`) or set FC to your compiler")
                raise RuntimeError(f"Fortran compiler `{fc}` not found: {hint}.")
            # always from scratch: a library left from another compiler or version must never be reused
            subprocess.check_call(["make", "-C", HERE, "clean"])
            subprocess.check_call(["make", "-C", HERE])
            _check_loadable()
        super().run()


def _check_loadable():
    """The library must load into this Python: on Apple silicon an x86_64 (Rosetta) Python cannot
    load the arm64 library that Homebrew's gfortran builds, and the error would only show later."""
    import ctypes
    import glob
    import platform
    libs = glob.glob(os.path.join(HERE, "python", "pyxscat", "libpyxscat.*"))
    if not libs:
        return
    try:
        ctypes.CDLL(libs[0])
    except OSError as e:
        msg = str(e)
        if "architecture" in msg or "mach-o" in msg.lower():
            raise RuntimeError(
                f"the compiled library does not match this Python ({platform.machine()}): {msg}\n"
                "On Apple silicon with an Intel (x86_64) conda, make a native environment instead:\n"
                "    CONDA_SUBDIR=osx-arm64 conda create -n thistle python=3.12 -y\n"
                "    conda activate thistle && conda config --env --set subdir osx-arm64\n"
                "and install again.") from e
        import warnings
        warnings.warn(f"the compiled library could not be loaded here: {msg}")


class BinaryDistribution(Distribution):
    """The wheel contains a compiled library, so it is platform specific."""

    def has_ext_modules(self):
        return True


setup(cmdclass={"build_py": BuildWithFortran}, distclass=BinaryDistribution)
