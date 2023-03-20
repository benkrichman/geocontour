# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= ~/Library/Python/3.9/bin/sphinx-build #-3.9
SPHINXAPI     ?= ~/Library/Python/3.9/bin/sphinx-apidoc
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

#API documentation build
#options:
#	-M put module docs before sub-docs
#	-f force-overwrite existing files
#	-e make separate pages for each module
buildapi:
	@$(SPHINXAPI) -Mfe -o source/ ../geocontour
	@echo "Auto-generated API documentation done"
