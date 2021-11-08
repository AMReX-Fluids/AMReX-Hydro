#!/bin/bash -v
set -e # Exit with nonzero exit code if anything fails

# Doxygen
echo "Build the Doxygen documentation"
doxygen ./Docs/Doxygen/doxygen.conf &> doxygen.out


# Sphinx
#
# Sphinx with Breathe ( Doxygen doc integration ):
# Requires: Breathe to be installed and parent directory added to PATH
# see https://github.com/michaeljones/breathe

echo "Build the Sphinx documentation"
cd ./Docs
make html

