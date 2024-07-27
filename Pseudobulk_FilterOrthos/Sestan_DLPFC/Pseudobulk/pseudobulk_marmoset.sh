#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cores 1
#SBATCH --mem-per-cpu 128GB

python pseudobulk_Sestan_2022_DLPFC.py marmoset subclass
python pseudobulk_Sestan_2022_DLPFC.py marmoset subtype