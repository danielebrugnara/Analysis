[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10462604.svg)](https://doi.org/10.5281/zenodo.10462604)

# Analysis
Analysis of my PhD project

Prerequisites ROOT and NPTOOL

Compile and run as follows:

``cd Analysis``

``mkdir build``

``cd build``

``cmake [options] ..``

``make install``

``./bin/analysis [n_threads]``

# Infos:


| Folder       | Sub-Folder                | Description                                                                    |
|--------------|---------------------------|--------------------------------------------------------------------------------|
| ./Lib        |                           | Main analysis libraries (data sorting, gates, calibrations and corrections)    |
| ./Utilities/ | SimulationAnalysis/Agata  | Gamma detection simulation. Tested and checked with real source data           |
| ./Utilities/ | SimulationAnalysis/Mugast | Silicon detector simulation. Needs both data and simulation to work            |
| ./Utilities/ | ./Utilities/Tester        | Utility for testing                                                            |
| ./Utilities/ | ./VamosFPOptimization     | Root macro to align charge states to improve reaction fragment identification  |
| ./Utilities/ | ./Scripts                 | Various random scripts                                                         |
| ./Utilities/ | ./AgataEfficiency         | Compute efficiency with Eu source using gamma spectra                          |
| ./Utilities/ | ./TargetThickness         | Mathematical model for target thickness and deformation (Mathematica)          |
| ./Utilities/ | ./MaxLikelihood           | Max likelihood calculation macro to compute angular distribution deconvolution |
| ./Utilities/ | ./Reactions               | Utilities for reaction calculations with FRESCO                                |
