#!/bin/bash -v
set -e # Exit with nonzero exit code if anything fails

# Doxygen
echo "Building the Doxygen documentation"

# Get AMReX Doxygen Tagfile
wget https://amrex-codes.github.io/amrex/docs_xml/doxygen/amrex-doxygen-web.tag.xml -O amrex-doxygen-web.tag.xml

# Render HTML Docs
doxygen ./Docs/Doxygen/doxygen.conf &> doxygen.out



