# SPDX-License-Identifier: Apache-2.0
"""Sphinx configuration for the THISTLE documentation.

Two generators feed one site:

* Sphinx (autodoc + MyST) documents the Python packages ``thistle`` and ``pyxscat`` and
  renders the Markdown pages in this directory.
* FORD documents the Fortran kernels in ``src/fortran``.  It is run from ``builder-inited``
  (see :func:`build_ford`) into ``docs/_ford/fortran`` and that tree is copied into the HTML
  output as ``fortran/`` through ``html_extra_path``.  Set ``THISTLE_DOCS_SKIP_FORD=1`` to
  skip the FORD step (e.g. while iterating on the Python pages).

Only ``docs/`` is scanned by Sphinx and only ``src/fortran`` by FORD: the virtual environment
at the repository root, ``build/`` and ``python/thistle.egg-info`` are never read.
"""
import glob
import os
import subprocess
import sys
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "python"))

with open(os.path.join(ROOT, "pyproject.toml"), "rb") as f:
    _meta = tomllib.load(f)["project"]

# ------------------------------------------------------------------------------ project
project = "THISTLE"
author = ", ".join(a["name"] for a in _meta["authors"])
copyright = "2026, Andrés Moreno Carrascosa, Patrick Y. Wang"
release = _meta["version"]
version = release

# ------------------------------------------------------------------------------ general
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_parser",
]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
master_doc = "index"
# _build and _ford are outputs; ford.md is FORD's project file, not a Sphinx page
exclude_patterns = ["_build", "_ford", "ford.md", "Thumbs.db", ".DS_Store", "**/__pycache__"]
templates_path = []
nitpicky = False

# MyST: the existing Markdown pages (QUICKSTART, THEORY, ...) are included as they are
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist", "attrs_inline"]
myst_heading_anchors = 3

# autodoc / napoleon
autodoc_member_order = "bysource"
autodoc_default_options = {"members": True, "undoc-members": True, "show-inheritance": True}
autodoc_typehints = "description"
autosummary_generate = False
napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = True
# optional run-time dependencies are not needed to render the docstrings
autodoc_mock_imports = ["pyscf", "scipy", "xraydb", "h5py"]
# the compiled kernel library: mocked only when it has not been built
if not os.environ.get("PYXSCAT_LIB") and not glob.glob(os.path.join(ROOT, "python", "pyxscat", "libpyxscat.*")):
    autodoc_mock_imports.append("pyxscat._lib")

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

# ------------------------------------------------------------------------------ HTML
from better import better_theme_path  # noqa: E402  (sphinx-better-theme)

html_theme = "better"
html_theme_path = [better_theme_path]
html_baseurl = "https://thistle.kirrander.group/"
html_title = f"THISTLE {release}"
html_short_title = "THISTLE"
html_theme_options = {
    "rightsidebar": False,
    "showheader": True,
    "showrelbartop": True,
    "showrelbarbottom": True,
    "linktotheme": False,
    "sidebarwidth": "16rem",
}
html_sidebars = {"**": ["globaltoc.html", "localtoc.html", "searchbox.html"]}
html_show_sourcelink = False
html_last_updated_fmt = "%Y-%m-%d"

# ------------------------------------------------------------------------------ FORD
FORD_PROJECT = os.path.join(HERE, "ford.md")
FORD_OUT = os.path.join(HERE, "_ford")          # FORD writes to _ford/fortran (see ford.md)
# html_extra_path is validated when the configuration is read, before the build starts,
# so the directory has to exist already; its content is copied to <outdir>/fortran/.
os.makedirs(os.path.join(FORD_OUT, "fortran"), exist_ok=True)
html_extra_path = [FORD_OUT]


def build_ford(app):
    """Run FORD on src/fortran so that its site is ready to be copied into the HTML output."""
    from sphinx.util import logging

    logger = logging.getLogger(__name__)
    if os.environ.get("THISTLE_DOCS_SKIP_FORD"):
        logger.info("FORD step skipped (THISTLE_DOCS_SKIP_FORD is set)")
        return
    logger.info("running FORD on src/fortran -> docs/_ford/fortran")
    subprocess.run([sys.executable, "-m", "ford", "--quiet", FORD_PROJECT], cwd=HERE, check=True)


def setup(app):
    app.connect("builder-inited", build_ford)
