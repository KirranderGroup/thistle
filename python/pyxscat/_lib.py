# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Andrés Moreno Carrascosa
"""Loading of the compiled PyXSCAT library and the raw ctypes signatures."""
import ctypes
import glob
import os

import numpy as np
from numpy.ctypeslib import ndpointer

_here = os.path.dirname(os.path.abspath(__file__))


def _find_library():
    env = os.environ.get("PYXSCAT_LIB")
    if env:
        return env
    hits = sorted(glob.glob(os.path.join(_here, "libpyxscat.*")))
    if not hits:
        raise ImportError(
            "The compiled PyXSCAT library was not found next to the Python package. "
            "Run `make` in the repository root, or set PYXSCAT_LIB to its path."
        )
    return hits[0]


_dbl = ndpointer(dtype=np.float64, flags="F_CONTIGUOUS")
_int = ndpointer(dtype=np.int32, flags="F_CONTIGUOUS")

try:
    lib = ctypes.CDLL(_find_library())
except OSError as exc:                       # typically a missing compiler runtime library
    raise ImportError(
        f"The PyXSCAT library could not be loaded ({exc}). If it was built with Intel ifx/ifort, "
        "load the oneAPI environment (`source /opt/intel/oneapi/setvars.sh` or `module load intel`; with "
        "conda's ifx: `export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`) "
        "so that the OpenMP runtime libiomp5 is found; with gfortran, make sure its runtime "
        "(libgfortran, libgomp) is on the library path.") from exc

lib.pyxscat_two_electron_j0.restype = None
lib.pyxscat_two_electron_j0.argtypes = [
    ctypes.c_int, ctypes.c_int, ctypes.c_int,             # nprim, ncontr, nq
    _dbl, _dbl, _dbl, _dbl,                               # x, y, z, ga
    _int, _int, _int, _int, _int,                         # l, m, n, group, contr
    _dbl, _dbl,                                           # coeff, q
    ctypes.c_int,                                         # mode
    ctypes.c_void_p,                                      # dm2 (or NULL)
    _dbl, _dbl,                                           # da, db
    ctypes.c_double, ctypes.c_double,                     # wc, wx
    ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,       # mode 3: nkc, Ac(k,i,j), Bc(k,i,j)  (k fastest)
    ctypes.c_int, ctypes.c_void_p, ctypes.c_void_p,       #         nkx, Ax, Bx
    ctypes.c_void_p,                                      #         dgemm function pointer (or NULL)
    ctypes.c_double, ctypes.c_double, ctypes.c_double,    # cutoffcentre, cutoffz, cutoffmd
    ctypes.c_double,                                      # screen_tol
    ctypes.c_int,                                         # lorder (0: j0, 2: j2 P2)
    ctypes.c_int,                                         # nthreads
    _dbl, _dbl,                                           # res, stats
    ctypes.POINTER(ctypes.c_int),                         # ierr
]

lib.pyxscat_hermite_table.restype = None
if hasattr(lib, "pyxscat_rdm12_bitwise_groups"):
    lib.pyxscat_rdm12_bitwise_groups.restype = None
    lib.pyxscat_rdm12_bitwise_groups.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int,                     # ndet, norb, nw
        ndpointer(dtype=np.int64, flags="C_CONTIGUOUS"),              # alpha words (ndet, nw)
        ndpointer(dtype=np.int64, flags="C_CONTIGUOUS"),              # beta words
        ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),            # coefficients
        ctypes.c_int, ctypes.c_int,                                   # nint, ng
        ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),              # gstart (ng+1), 1-based
        ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),              # grp (ndet), 1-based
        ctypes.c_int,                                                 # nnb
        ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),              # nbstart (ng+1), 1-based
        ndpointer(dtype=np.int32, flags="C_CONTIGUOUS"),              # nb (nnb), 1-based group numbers
        ctypes.c_int,                                                 # nthreads
        ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),            # dm1
        ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),            # dm2
        ctypes.POINTER(ctypes.c_int)]

if hasattr(lib, "pyxscat_rdm12_bitwise"):
    lib.pyxscat_rdm12_bitwise.restype = None
    lib.pyxscat_rdm12_bitwise.argtypes = [
        ctypes.c_int, ctypes.c_int, ctypes.c_int,                     # ndet, norb, nw (words per string)
        ndpointer(dtype=np.int64, flags="C_CONTIGUOUS"),              # alpha words (ndet, nw)
        ndpointer(dtype=np.int64, flags="C_CONTIGUOUS"),              # beta words
        ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),            # coefficients
        ctypes.c_int,                                                 # nthreads
        ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),            # dm1[p, q]
        ndpointer(dtype=np.float64, flags="C_CONTIGUOUS"),            # dm2[p, q, r, s]
        ctypes.POINTER(ctypes.c_int),                                 # ierr
    ]
lib.pyxscat_hermite_table.argtypes = [ctypes.c_int, ctypes.c_int, _int, _dbl, _dbl, _dbl]

lib.pyxscat_pair_transform.restype = None
lib.pyxscat_pair_transform.argtypes = [
    ctypes.c_int, ctypes.c_int, ctypes.c_int,             # nprim, nao, nq
    _dbl, _dbl, _dbl, _dbl,                               # x, y, z, ga
    _int, _int, _int, _int, _int,                         # l, m, n, group, contr
    _dbl, _dbl,                                           # coeff, qv(3, nq)
    ctypes.c_double, ctypes.c_int,                        # screen, nthreads
    ndpointer(dtype=np.complex128, flags="F_CONTIGUOUS"),  # f(nq, nao, nao)
    ctypes.POINTER(ctypes.c_int),                         # ierr
]

lib.pyxscat_electron_nucleus.restype = None
lib.pyxscat_electron_nucleus.argtypes = [
    ctypes.c_int, ctypes.c_int, ctypes.c_int,             # nprim, nao, nq
    _dbl, _dbl, _dbl, _dbl,                               # x, y, z, ga
    _int, _int, _int, _int, _int,                         # l, m, n, group, contr
    _dbl, _dbl,                                           # coeff, q
    ctypes.c_int, _dbl, _dbl,                             # natom, charges, positions (3, natom)
    _dbl,                                                 # dm1 (nao, nao)
    ctypes.c_int, ctypes.c_double, ctypes.c_int,          # lorder, screen, nthreads
    _dbl, ctypes.POINTER(ctypes.c_int),                   # res, ierr
]

lib.pyxscat_elastic_pairs.restype = None
lib.pyxscat_elastic_pairs.argtypes = [
    ctypes.c_int, ctypes.c_int, ctypes.c_int,             # nprim, nao, nq
    _dbl, _dbl, _dbl, _dbl,                               # x, y, z, ga
    _int, _int, _int, _int, _int,                         # l, m, n, group, contr
    _dbl, _dbl,                                           # coeff, q
    _dbl, _dbl,                                           # da, db (nao, nao)
    ctypes.c_int, ctypes.c_double, ctypes.c_int,          # lorder, screen, nthreads
    _dbl, _dbl, ctypes.POINTER(ctypes.c_int),             # res, stats(2), ierr
]

lib.pyxscat_radial_series.restype = None
lib.pyxscat_radial_series.argtypes = [ctypes.c_int, _dbl, ctypes.c_double, ctypes.c_int, _dbl, _dbl]

try:
    lib.pyxscat_build_info.restype = None
    lib.pyxscat_build_info.argtypes = [ctypes.c_char_p, ctypes.c_int]
except AttributeError:                       # library built before this entry point existed
    pass

lib.pyxscat_max_threads.restype = ctypes.c_int
lib.pyxscat_max_threads.argtypes = []

ERRORS = {
    1: "empty basis or q grid",
    2: "a two-particle density matrix is required for this calculation",
    3: "unknown density provider",
    4: "angular momentum too high for the quartet kernel (4 l_max + 3 > 32; tested up to h functions)",
    5: "primitive shells (groups) must be numbered 1..ng and be contiguous",
    6: "contracted-function index out of range",
    7: "all members of a primitive shell must share centre, exponent and angular momentum",
    8: "only the Legendre orders l = 0 and l = 2 are implemented",
}


def blas_dgemm_pointer():
    """Address of the BLAS ``dgemm`` SciPy ships (Fortran calling convention, 32-bit integers),
    for the factored-density kernel; None when unavailable or when a quick check fails."""
    try:
        from scipy.linalg import cython_blas
        cap = cython_blas.__pyx_capi__["dgemm"]
        api = ctypes.pythonapi
        api.PyCapsule_GetName.restype = ctypes.c_char_p
        api.PyCapsule_GetName.argtypes = [ctypes.py_object]
        api.PyCapsule_GetPointer.restype = ctypes.c_void_p
        api.PyCapsule_GetPointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
        ptr = api.PyCapsule_GetPointer(cap, api.PyCapsule_GetName(cap))
        if not ptr:
            return None
        # check the ABI on a small product
        proto = ctypes.CFUNCTYPE(None, ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int),
                                 ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
                                 ctypes.POINTER(ctypes.c_double), ctypes.c_void_p, ctypes.POINTER(ctypes.c_int),
                                 ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
                                 ctypes.c_void_p, ctypes.POINTER(ctypes.c_int))
        fn = proto(ptr)
        a = np.asfortranarray(np.arange(6.0).reshape(3, 2))
        b = np.asfortranarray(np.arange(12.0).reshape(3, 4) + 1)
        c = np.zeros((2, 4), order="F")
        i = lambda v: ctypes.byref(ctypes.c_int(v))
        d = lambda v: ctypes.byref(ctypes.c_double(v))
        fn(b"T", b"N", i(2), i(4), i(3), d(1.0), a.ctypes.data, i(3), b.ctypes.data, i(3), d(0.0), c.ctypes.data, i(2))
        if not np.allclose(c, a.T @ b):
            return None
        return ptr
    except Exception:
        return None


DGEMM = blas_dgemm_pointer()
