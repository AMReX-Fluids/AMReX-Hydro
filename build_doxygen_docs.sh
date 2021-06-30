#!/bin/bash -v
set -e # Exit with nonzero exit code if anything fails

# Doxygen
echo "Build the Doxygen documentation"
doxygen ./Docs/Doxygen/doxygen.conf &> doxygen.out


 
