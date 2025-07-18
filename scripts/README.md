# KBfit Development Scripts

This directory contains utility scripts for KBfit development.

## Available Scripts

### `format_code.sh`
Formats all source code using clang-format.

**Usage:**
```bash
./scripts/format_code.sh
```

**Requirements:**
- clang-format installed
- Must be run from KBfit root directory

**Features:**
- Automatically finds and formats all source files
- Shows formatting changes with git diff
- Provides guidance on committing changes

### `build_docs.sh`
Builds both Doxygen and Sphinx documentation locally.

**Usage:**
```bash
# Build documentation
./scripts/build_docs.sh

# Build and serve documentation locally
./scripts/build_docs.sh --serve
```

**Requirements:**
- doxygen
- sphinx, sphinx-rtd-theme, breathe, myst-parser

**Features:**
- Checks for required dependencies
- Builds Doxygen XML for API documentation
- Builds Sphinx HTML documentation
- Optional local server for viewing results

## Development Workflow

### Before Committing
1. Format your code:
   ```bash
   ./scripts/format_code.sh
   ```

2. Build and test documentation:
   ```bash
   ./scripts/build_docs.sh
   ```

3. Run tests:
   ```bash
   cd build
   ctest
   ```

### Local Documentation Development
1. Build documentation with live server:
   ```bash
   ./scripts/build_docs.sh --serve
   ```

2. Open http://localhost:8000 in your browser

3. Make changes to documentation files

4. Rebuild documentation:
   ```bash
   ./scripts/build_docs.sh
   ```

## CI/CD Integration

The GitHub Actions workflow automatically:
- Checks code formatting
- Builds documentation
- Deploys to GitHub Pages
- Auto-formats pull requests

These scripts mirror the CI/CD process for local development.

## Troubleshooting

### Common Issues

**clang-format not found:**
```bash
sudo apt-get install clang-format
# or
conda install clang-format
```

**Sphinx not found:**
```bash
pip install sphinx sphinx-rtd-theme breathe myst-parser
```

**Documentation build fails:**
- Check that you're in the KBfit root directory
- Ensure all dependencies are installed
- Check for syntax errors in documentation files

### Getting Help

- Check the error messages for specific issues
- Ensure you're in the correct conda environment: `conda activate spectrum`
- Verify all dependencies are installed
- Check that you're in the KBfit root directory

## Adding New Scripts

When adding new development scripts:

1. Place them in the `scripts/` directory
2. Make them executable: `chmod +x script_name.sh`
3. Add appropriate error checking and user feedback
4. Update this README with usage instructions
5. Consider adding them to the CI/CD workflow if needed