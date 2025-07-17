SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Build HTML documentation
html:
	@$(SPHINXBUILD) $(SPHINXOPTS) -b html $(SOURCEDIR) $(BUILDDIR)/html

# Build Doxygen documentation
doxygen:
	@doxygen Doxyfile

# Build both Doxygen and Sphinx documentation
docs: doxygen html

# Clean build artifacts
clean:
	@rm -rf $(BUILDDIR)

# Serve documentation locally
serve:
	@echo "Starting local documentation server..."
	@echo "Open http://localhost:8000 in your browser"
	@cd $(BUILDDIR)/html && python -m http.server 8000

# Format code using clang-format
format:
	@../scripts/format_code.sh

# Build documentation using build script
build-docs:
	@../scripts/build_docs.sh

# Build and serve documentation
build-serve:
	@../scripts/build_docs.sh --serve

.PHONY: html doxygen docs clean serve format build-docs build-serve