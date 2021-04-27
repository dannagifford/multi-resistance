How to run the simulations:
First compile simulation_code.cpp with a C++ compiler. Requires the boost library to be installed so that the .ini files describing the environments can be loaded. On Ubuntu you can install this with 'apt-get install libboost-all-dev' (though you might only want the specific ones for your system).

g++ simulation_code.cpp -o simulation

Then call the simulation program with the .ini environment file as an argument. It will loop over initial mutator frequencies specified as an array in the source code (currently 0, 0.05, 0.1, 0.3, 0.5) e.g.

simulation rif/rif.ini

It will segfault if you forget to specify the .ini file.

The output will be saved as a .txt file with the same prefix in the same folder as the input file. Data are saved for each 'hour' of the simulation, plus an extra row is generated for the initial population conditions. Hence, there are 133 (=22*6+1) rows per replicate of the simulation, and 133*replicate rows in the output file. There is one column for each of the 8 possible genotypes, arranged as S, R, N, D, S', R', N', D', where prime indicates the mutator genetic background.

Plots and analysis can then be performed in R with 'plots.R' and 'model-fitting-simulations.R'

Rscript plots.R

Rscript model-fitting-simulations.R
