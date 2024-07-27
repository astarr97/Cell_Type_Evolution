#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cores 1
#SBATCH --mem-per-cpu 128GB

python pseudobulk.py chimp cross_species_cluster
python pseudobulk.py chimp subclass