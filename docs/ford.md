---
project: PyXSCAT kernels
summary: The Fortran kernels of THISTLE, analytic X-ray and electron scattering integrals over Gaussian basis functions.
author: Andrés Moreno Carrascosa, Mats Simmermacher, Patrick Yuheng Wang, Eirik Moen Liane, Joseph Cooper, Nikola Zotev, Thomas Northey and Adam Kirrander
project_github: https://github.com/KirranderGroup/thistle
version: 0.2.0.dev0
src_dir: ../src/fortran
output_dir: ./_ford/fortran
extensions: f90
preprocess: false
display: public
         protected
         private
source: true
proc_internals: true
search: true
graph: false
print_creation_date: true
predocmark: >
docmark: !
---

The Fortran kernels of [THISTLE](../index.html): the shared library `libpyxscat` that the
Python package `pyxscat` calls through the C entry points of module `xs_capi` (ctypes;
see [the Python API](../api/pyxscat.html)).  Apache License 2.0.

Documentation comments use the FORD conventions: `!>` in the line(s) before a module,
procedure or variable, or `!!` on the line(s) after its declaration.  Plain `!` comments
are not collected.  `graph: false` above is because the call graphs need Graphviz (`dot`);
switch it on where that is installed.
