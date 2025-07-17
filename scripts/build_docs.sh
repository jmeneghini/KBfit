#!/bin/bash

# Documentation building script for KBfit
# This script builds both Doxygen and Sphinx documentation locally

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if we're in the KBfit root directory
if [ ! -f "docs/conf.py" ]; then
    echo -e "${RED}Error: docs/conf.py not found${NC}"
    echo "Please run this script from the KBfit root directory"
    exit 1
fi

echo -e "${BLUE}Building KBfit documentation...${NC}"

# Create docs directory if it doesn't exist
mkdir -p docs/_build

# Check dependencies
echo -e "${YELLOW}Checking dependencies...${NC}"

# Check for doxygen
if ! command -v doxygen &> /dev/null; then
    echo -e "${RED}Error: doxygen is not installed${NC}"
    echo "Please install doxygen: conda install doxygen"
    exit 1
fi

# Check for sphinx
if ! command -v sphinx-build &> /dev/null; then
    echo -e "${RED}Error: sphinx-build is not installed${NC}"
    echo "Please install sphinx: pip install sphinx sphinx-rtd-theme breathe myst-parser"
    exit 1
fi

# Build Doxygen documentation
echo -e "${YELLOW}Building Doxygen documentation...${NC}"
cd docs
doxygen Doxyfile

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Doxygen documentation built successfully${NC}"
else
    echo -e "${RED}Error: Doxygen build failed${NC}"
    exit 1
fi

# Build Sphinx documentation
echo -e "${YELLOW}Building Sphinx documentation...${NC}"
sphinx-build -b html . _build/html

if [ $? -eq 0 ]; then
    echo -e "${GREEN}Sphinx documentation built successfully${NC}"
else
    echo -e "${RED}Error: Sphinx build failed${NC}"
    exit 1
fi

cd ..

# Check if documentation was built
if [ -f "docs/_build/html/index.html" ]; then
    echo -e "${GREEN}Documentation build complete!${NC}"
    echo -e "${BLUE}Open the documentation:${NC}"
    echo "  Browser: file://$(pwd)/docs/_build/html/index.html"
    echo "  Or run: python -m http.server 8000 --directory docs/_build/html"
else
    echo -e "${RED}Error: Documentation build failed${NC}"
    exit 1
fi

# Optional: serve the documentation locally
if [ "$1" = "--serve" ]; then
    echo -e "${BLUE}Starting local documentation server...${NC}"
    echo "Open http://localhost:8000 in your browser"
    echo "Press Ctrl+C to stop the server"
    cd docs/_build/html
    python -m http.server 8000
fi