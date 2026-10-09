# Sphinx documentation.  `make html` also runs FORD on src/fortran (from conf.py) and
# places its site at _build/html/fortran/.  THISTLE_DOCS_SKIP_FORD=1 make html skips FORD.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

ford:
	ford ford.md

clean:
	rm -rf "$(BUILDDIR)" _ford

.PHONY: help ford clean Makefile

%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
