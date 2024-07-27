We used the commands in commands_parameter_sweeps.sh to compare cell type proportion and divergence in gene expression.
Running driver.sh would create individual scripts for each command in commands_parameter_sweeps.sh.

Note: some runs (when controlling for tau or expression level for high clustering resolution for Sestan_DLPFC or Allen_MTG) did not finish in 7 days
In that case, we reran starting at the first unfinished iteration and combined the two files to make a complete file

All commands use the script parameter_sweep_general_new.py, which goes through 100 iterations across all combinations of parameters (specified in the script) and computes relationship between cell type divergence and cell type proportion, stratifying by something if desired.
