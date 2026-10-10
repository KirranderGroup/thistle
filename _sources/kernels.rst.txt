Fortran kernels
===============

The numerical work of THISTLE is done by ten Fortran modules in ``src/fortran``, compiled by
the top-level ``Makefile`` into ``python/pyxscat/libpyxscat.so`` (``.dylib`` on macOS).
Python reaches them only through the ``bind(C)`` entry points of module ``xs_capi``, whose
ctypes signatures are declared in :mod:`pyxscat._lib`.

The modules, their routines, arguments, source listings and the comments attached to them are
documented by `FORD <https://forddocs.readthedocs.io/>`_, built into this site:

* `Fortran API reference (FORD) <fortran/index.html>`_

Modules
-------

================  ===============================================================================
Module            Contents
================  ===============================================================================
``xs_hermite``    Hermite (McMurchie–Davidson) expansion coefficients for any angular momentum
``xs_bessel_j0``  Spherical Bessel functions, derivatives and the recurrences of the j0 kernel
``xs_p0``         Closed forms at coincident centres (h = 0) for j0, j1 and j2
``xs_radial``     Radial functions j_n(x)/x^n and the Hermite radial coefficients, stable at any x
``xs_kernel_j0``  Screened shell-quartet two-electron kernel; density providers (dense, mean-field, factored)
``xs_pairft``     Analytic Fourier transforms of basis-function pairs at fixed q vectors
``xs_elastic``    Elastic scattering by density-weighted pairs
``xs_nuclear``    Electron–nucleus integrals for electron scattering
``xs_rdm``        Bitwise 1-RDM and 2-RDM builders from determinant expansions
``xs_capi``       The C-compatible entry points called from Python
================  ===============================================================================

Building the Fortran reference
------------------------------

FORD is run automatically by ``conf.py`` at the start of every Sphinx build (``make html``
in ``docs/``) and writes to ``docs/_ford/fortran``; Sphinx then copies that tree into the HTML
output as ``fortran/``.  Set ``THISTLE_DOCS_SKIP_FORD=1`` to skip it.  To run FORD on its own:

.. code-block:: bash

   cd docs && ford ford.md

Writing the Fortran comments
----------------------------

FORD only collects comments written with its markers.  The kernels use ``!>`` on the lines
*before* a module, procedure or variable; ``!!`` on the lines *after* a declaration also works
(for arguments, for instance).  Plain ``!`` comments are shown in the source listings only, and
ruler lines (``! -----``) between a ``!>`` block and the declaration are ignored.  The text is
Markdown, which has three consequences worth knowing:

* The first paragraph is the summary shown in the lists.
* Equations, tables and option lists are indented by at least four spaces after ``!>``
  (``!>      I2(q) = ...``) with an empty ``!>`` line before and after, so they are rendered
  verbatim as code blocks.  Prose at one space (``!> text``) is wrapped as a paragraph.
* In prose, anything of the form ``<name ...>`` is taken for an HTML tag and disappears, and a
  first line of the form ``word: text`` is taken for FORD metadata.  Put bra-kets in backticks
  (```<I| ... |I>```) and start the block with a sentence.

The project file is ``docs/ford.md``.  FORD cannot extract the source listing of three
one-statement functions with a kind-prefixed result (``ncum``, ``dfact``,
``upward_threshold``); their descriptions and interfaces are documented all the same.
