#!/bin/bash

set -e 

##### External Dependencies (downloadable)

### Comet v2023.01.2 Download
# Delete previously downloaded Comet
rm -rf bin/comet.linux.exe
# Download Comet
wget -O bin/comet.linux.exe https://github.com/UWPR/Comet/releases/download/v2023.01.2/comet.linux.exe 


### ThermoRawFileParser 1.4.2 Download and extraction
# Delete previously downloaded ThermoRawFileParser
rm -rf bin/ThermoRawFileParser
# Download ThermoRawFileParser
wget -O trfp.zip https://github.com/compomics/ThermoRawFileParser/releases/download/v1.4.2/ThermoRawFileParser1.4.2.zip
# Extract archive
unzip trfp.zip -d bin/ThermoRawFileParser
# Delete downloaded archive
rm trfp.zip


### Percolator v3.0.6 Download and extraction
# Delete previously downloaded OPercolatorpenMS
rm -rf bin/percolator/
# Download percolator
wget -O percolator.deb https://github.com/percolator/percolator/releases/download/rel-3-06/percolator-v3-06-linux-amd64.deb
# Extract only the data-part of the deb package
ar x percolator.deb data.tar.gz
# Extract the data into the percolator-Folder
mkdir -p ./bin/percolator
tar -xf data.tar.gz -C ./bin/percolator
# Delete downloaded and extracted deb package
rm -rf percolator.deb data.tar.gz


##### ProtGraphCPP Traversal Compilation

### Build ProtGraphTraverseInt Runtime Estimator (dryrun)
rm -rf bin/ProtGraphTraverseIntSourceDryRun/build
cmake -B bin/ProtGraphTraverseIntSourceDryRun/build -S bin/ProtGraphTraverseIntSourceDryRun
# Some machines (like older ubuntus) need to specify the new gcc compiler via:  "-D CMAKE_C_COMPILER=gcc-11 -D CMAKE_CXX_COMPILER=g++-11"
OLDDIR=$(pwd)
cd bin/ProtGraphTraverseIntSourceDryRun/
cmake --build ./build
cd $OLDDIR


### Build ProtGraphTraverseFloat Runtime Estimator (dryrun)
rm -rf bin/ProtGraphTraverseFloatSourceDryRun/build
cmake -B bin/ProtGraphTraverseFloatSourceDryRun/build -S bin/ProtGraphTraverseFloatSourceDryRun
# Some machines (like older ubuntus) need to specify the new gcc compiler via:  "-D CMAKE_C_COMPILER=gcc-11 -D CMAKE_CXX_COMPILER=g++-11"
OLDDIR=$(pwd)
cd bin/ProtGraphTraverseFloatSourceDryRun/
cmake --build ./build
cd $OLDDIR


# Build ProtGraphTraverseInt FASTA-Generation (using limits, which are retrieved by the binary search)
rm -rf bin/ProtGraphTraverseIntSourceVarLimitter/build
cmake -B bin/ProtGraphTraverseIntSourceVarLimitter/build -S bin/ProtGraphTraverseIntSourceVarLimitter
# Some machines (like older ubuntus) need to specify the new gcc compiler via:  "-D CMAKE_C_COMPILER=gcc-11 -D CMAKE_CXX_COMPILER=g++-11"
OLDDIR=$(pwd)
cd bin/ProtGraphTraverseIntSourceVarLimitter/
cmake --build ./build
cd $OLDDIR


# Build ProtGraphTraverseInt FASTA-Generation (using limits, which are retrieved by the binary search)
rm -rf bin/ProtGraphTraverseFloatSourceVarLimitter/build
cmake -B bin/ProtGraphTraverseFloatSourceVarLimitter/build -S bin/ProtGraphTraverseFloatSourceVarLimitter
# Some machines (like older ubuntus) need to specify the new gcc compiler via:  "-D CMAKE_C_COMPILER=gcc-11 -D CMAKE_CXX_COMPILER=g++-11"
OLDDIR=$(pwd)
cd bin/ProtGraphTraverseFloatSourceVarLimitter/
cmake --build ./build
cd $OLDDIR


##### Python Dependencies (Python 3.9)

# Protein-Graph-Generation and other exports
pip install protgraph==0.3.9
# Needed for exporting SQLite via ProtGraph
pip install apsw==3.42.0.0
# Needed to generate compact Excel-Files (HeatMaps)
pip install XlsxWriter==3.0.3
# Used to load in comet-txt-files
pip install pandas
# Needed for generating various plots
pip install plotly
# Engine, needed to run all workflows
pip install nextflow
# Needed for parsing xmls (e.g. Percolator output)
pip install lxml


##### Make all files in bin executable (excluding sub-directories) to be visible by processes in Nextflow
chmod +x bin/*
