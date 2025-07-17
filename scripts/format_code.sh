#!/bin/bash

# Code formatting script for KBfit
# This script formats all source files using clang-format

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if clang-format is available
if ! command -v clang-format &> /dev/null; then
    echo -e "${RED}Error: clang-format is not installed${NC}"
    echo "Please install clang-format to format code"
    exit 1
fi

# Check if we're in the KBfit root directory
if [ ! -f ".clang-format" ]; then
    echo -e "${RED}Error: .clang-format file not found${NC}"
    echo "Please run this script from the KBfit root directory"
    exit 1
fi

echo -e "${GREEN}Formatting KBfit source code...${NC}"

# Find all source files and format them
find source -name "*.h" -o -name "*.cc" -o -name "*.cpp" -o -name "*.hpp" -o -name "*.c" -o -name "*.cxx" -o -name "*.hxx" | while read file; do
    echo "Formatting: $file"
    clang-format -i "$file"
done

echo -e "${GREEN}Code formatting complete!${NC}"

# Check if there are any changes
if git diff --quiet --exit-code; then
    echo -e "${GREEN}No formatting changes needed${NC}"
else
    echo -e "${YELLOW}Formatting changes detected:${NC}"
    git diff --stat
    echo ""
    echo "Review the changes and commit if they look good:"
    echo "  git add -A"
    echo "  git commit -m 'Format code with clang-format'"
fi