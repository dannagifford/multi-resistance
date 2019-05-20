#! /bin/bash
# First run to get data into needed formats
Rscript munging.R

# Then run for statistical models (will take time)
Rscript model-fitting.R

# Then run for plots
Rscript plots.R
