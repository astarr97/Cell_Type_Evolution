#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cores 1
#SBATCH --mem-per-cpu 128GB

mkdir Allen_M1_Subsamples_cross_species_cluster_label
mkdir Allen_M1_Subsamples_subclass_label

python process_for_ctp.py mouse cross_species_cluster_label Allen_M1_Subsamples_cross_species_cluster_label
python process_for_ctp.py human cross_species_cluster_label Allen_M1_Subsamples_cross_species_cluster_label
python process_for_ctp.py marmoset cross_species_cluster_label Allen_M1_Subsamples_cross_species_cluster_label

python process_for_ctp.py mouse subclass_label Allen_M1_Subsamples_subclass_label
python process_for_ctp.py human subclass_label Allen_M1_Subsamples_subclass_label
python process_for_ctp.py marmoset subclass_label Allen_M1_Subsamples_subclass_label

python convert_orthos_Allen_M1.py Allen_M1_NoSubsamples_subclass_label,Allen_M1_NoSubsamples_cross_species_cluster_label,Allen_M1_Subsamples_cross_species_cluster_label,Allen_M1_Subsamples_subclass_label Allen_M1_NoSubsamples_subclass_label

ml load system rclone

tar -cvzf Allen_M1_Subsamples_cross_species_cluster_label.tar.gz Allen_M1_Subsamples_cross_species_cluster_label
rclone copy Allen_M1_Subsamples_cross_species_cluster_label.tar.gz fraser:astarr/CellTypeProp_New_And_Pseudobulk/Allen_M1

tar -cvzf Allen_M1_Subsamples_subclass_label.tar.gz Allen_M1_Subsamples_subclass_label
rclone copy Allen_M1_Subsamples_subclass_label.tar.gz fraser:astarr/CellTypeProp_New_And_Pseudobulk/Allen_M1
