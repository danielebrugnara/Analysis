# Analysis
Analysis of my PhD project

Prerequisites ROOT and NPTOOL

Compile and run as follows:

cd Analysis

mkdir Build

cd Build

cmake ..

make install

./bin/analysis [n_threads]

# Infos:

./                                      -> main analysis
./Utilities/SimulationAnalysis/Agata    -> cryogenic target simulation
./Utilities/SimulationAnalysis/Mugast   -> mugast ang distr simulation
./Utilities/Tester                      -> Tester for shared libraries
./VamosFPOptimization                   -> Find shift to align reaction prods. in focal plane
./Scripts                               -> Used for femul replays on multiple machines    
./AgataEfficiency                       -> Compute efficiency with Eu 
./TargetThickness                       -> Mathematical model for target thickness
