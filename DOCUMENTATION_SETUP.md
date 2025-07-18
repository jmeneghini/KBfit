# KBfit Documentation and Code Quality Setup

## Overview

This document summarizes the comprehensive documentation and code quality system that has been set up for the KBfit project.

## What Was Implemented

### 1. Enhanced GitHub Actions Workflow (.github/workflows/docs.yml)

**Features:**
- **Automatic Documentation Generation**: Builds Doxygen and Sphinx documentation on every push
- **GitHub Pages Deployment**: Automatically deploys documentation to GitHub Pages
- **Code Formatting Checks**: Validates code formatting with clang-format
- **Auto-formatting**: Automatically formats code in pull requests
- **Multi-branch Support**: Deploys main branch to root, development branch to `/dev`

**Jobs:**
- `format-check`: Validates code formatting
- `build-docs`: Builds and deploys documentation
- `auto-format`: Automatically formats code in pull requests

### 2. Comprehensive Documentation System

**Sphinx Configuration (docs/conf.py):**
- Enhanced with modern extensions
- MathJax support for mathematical formulas
- Intersphinx linking to external documentation
- Read the Docs theme with custom styling
- Conditional myst-parser support

**Doxygen Configuration (docs/Doxyfile):**
- Comprehensive API documentation extraction
- Class diagrams and call graphs
- Source code browsing
- XML output for Sphinx integration

**Documentation Files:**
- `index.rst`: Main documentation homepage
- `installation.rst`: Comprehensive installation guide
- `getting_started.rst`: User guide and examples
- `physics_background.rst`: Theoretical foundation
- `contributing.rst`: Developer contribution guidelines
- `changelog.rst`: Version history and changes
- `api.rst`: API reference integration

### 3. Development Scripts

**scripts/format_code.sh:**
- Formats all source code with clang-format
- Shows formatting changes
- Provides commit guidance

**scripts/build_docs.sh:**
- Builds both Doxygen and Sphinx documentation
- Dependency checking
- Optional local server for testing

**scripts/README.md:**
- Usage instructions for all scripts
- Development workflow guidance
- Troubleshooting tips

### 4. Enhanced Makefile (docs/Makefile)

**Targets:**
- `docs`: Build both Doxygen and Sphinx documentation
- `html`: Build HTML documentation
- `doxygen`: Build Doxygen documentation
- `serve`: Serve documentation locally
- `format`: Format code
- `clean`: Clean build artifacts

## Key Optimizations Documented

### 1. Spectrum Fitting Performance Optimization

**Root Finding in CM Frame:**
- Eliminates redundant Ecm ↔ Elab coordinate transformations
- ~50% reduction in computational overhead
- Maintains numerical accuracy
- Uses `BoxQuantization::getDeltaEnergyPredictionsOptimized()`

**Performance Improvements:**
- Pre-allocated temporary vectors
- Cache-friendly memory access patterns
- Efficient parameter passing
- Minimized function call overhead

### 2. Comprehensive Code Documentation

**Added Documentation:**
- Detailed method documentation with performance notes
- Implementation details for optimization
- Physics background with mathematical formulas
- Call chain documentation for critical paths

## Usage Instructions

### Local Development

1. **Format Code:**
   ```bash
   ./scripts/format_code.sh
   ```

2. **Build Documentation:**
   ```bash
   ./scripts/build_docs.sh
   ```

3. **Serve Documentation Locally:**
   ```bash
   ./scripts/build_docs.sh --serve
   ```

4. **Using Make:**
   ```bash
   cd docs
   make docs        # Build all documentation
   make serve       # Serve locally
   make format      # Format code
   ```

### GitHub Actions

**Automatic Triggers:**
- **Push to main**: Deploys documentation to GitHub Pages
- **Push to spectrum_fits**: Deploys to `/dev` subdirectory
- **Pull Request**: Checks formatting and optionally auto-formats

**Manual Triggers:**
- Documentation can be manually triggered from GitHub Actions

## File Structure

```
KBfit/
├── .github/
│   └── workflows/
│       └── docs.yml                 # Enhanced CI/CD workflow
├── docs/
│   ├── conf.py                      # Sphinx configuration
│   ├── Doxyfile                     # Doxygen configuration
│   ├── Makefile                     # Enhanced build targets
│   ├── index.rst                    # Main documentation
│   ├── installation.rst             # Installation guide
│   ├── getting_started.rst          # User guide
│   ├── physics_background.rst       # Physics documentation
│   ├── contributing.rst             # Contribution guidelines
│   ├── changelog.rst                # Version history
│   ├── api.rst                      # API reference
│   └── readme.rst                   # README integration
├── scripts/
│   ├── format_code.sh               # Code formatting script
│   ├── build_docs.sh                # Documentation build script
│   └── README.md                    # Script documentation
└── DOCUMENTATION_SETUP.md           # This file
```

## Benefits

### For Developers

1. **Consistent Code Style**: Automatic formatting ensures consistent code style
2. **Easy Documentation**: Simple scripts for building and serving documentation
3. **Fast Feedback**: Immediate feedback on code formatting in PRs
4. **Comprehensive Guides**: Clear contribution and development guidelines

### For Users

1. **Professional Documentation**: High-quality, searchable documentation
2. **Always Up-to-Date**: Automatically updated on every change
3. **Multiple Formats**: Web-based documentation with mathematical formulas
4. **Easy Access**: Available on GitHub Pages

### For Maintainers

1. **Automated Quality Control**: Automatic code formatting and documentation builds
2. **Reduced Manual Work**: Automatic deployment and formatting
3. **Better Documentation**: Comprehensive API and theoretical documentation
4. **Version Control**: Automatic versioning and changelog

## Mathematical Documentation

The physics background documentation includes properly formatted mathematical formulas:

- **Quantization Condition**: det[1 - B^(P)(E_cm) K̃] = 0
- **Box Matrix Definition**: Full mathematical formulation
- **K-matrix Relationships**: Detailed mathematical relationships
- **Cayley Transformations**: Complete mathematical definitions

## Performance Documentation

Comprehensive documentation of recent optimizations:

- **Before/After Comparisons**: Clear performance improvement metrics
- **Implementation Details**: Step-by-step optimization explanations
- **Call Chain Documentation**: Complete function call sequences
- **Memory Access Patterns**: Cache-friendly optimization descriptions

## Next Steps

1. **Test the GitHub Actions**: Push changes to test the automated workflow
2. **Verify GitHub Pages**: Check that documentation deploys correctly
3. **Team Training**: Introduce team members to the new documentation system
4. **Continuous Improvement**: Gather feedback and enhance documentation

## Troubleshooting

### Common Issues

1. **Sphinx Build Failures**: Check conda environment and dependencies
2. **Formatting Issues**: Ensure clang-format is installed and accessible
3. **GitHub Pages**: Verify repository settings for GitHub Pages
4. **Mathematical Formulas**: Ensure MathJax is properly configured

### Getting Help

- Check script output for specific error messages
- Verify conda environment activation
- Test documentation builds locally before pushing
- Review GitHub Actions logs for deployment issues

This comprehensive setup provides a professional, maintainable documentation and code quality system for the KBfit project.