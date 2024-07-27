#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --ntasks=1
#SBATCH --cores 1
#SBATCH --mem-per-cpu 64GB

mkdir Allen_MTG_Subsamples_cross_species_cluster
mkdir Allen_MTG_Subsamples_subclass

python process_for_ctp.py rhesus cross_species_cluster Allen_MTG_Subsamples_cross_species_cluster
python process_for_ctp.py gorilla cross_species_cluster Allen_MTG_Subsamples_cross_species_cluster
python process_for_ctp.py chimp cross_species_cluster Allen_MTG_Subsamples_cross_species_cluster
python process_for_ctp.py human cross_species_cluster Allen_MTG_Subsamples_cross_species_cluster
python process_for_ctp.py marmoset cross_species_cluster Allen_MTG_Subsamples_cross_species_cluster

python process_for_ctp.py rhesus subclass Allen_MTG_Subsamples_subclass
python process_for_ctp.py gorilla subclass Allen_MTG_Subsamples_subclass
python process_for_ctp.py chimp subclass Allen_MTG_Subsamples_subclass
python process_for_ctp.py human subclass Allen_MTG_Subsamples_subclass
python process_for_ctp.py marmoset subclass Allen_MTG_Subsamples_subclass

python convert_orthos_Allen_MTG.py Allen_MTG_Subsamples_cross_species_cluster,Allen_MTG_Subsamples_subclass

ml load system rclone

tar -cvzf Allen_MTG_Subsamples_cross_species_cluster.tar.gz Allen_MTG_Subsamples_cross_species_cluster
rclone copy Allen_MTG_Subsamples_cross_species_cluster.tar.gz fraser:astarr/CellTypeProp_New_And_Pseudobulk/Allen_MTG

tar -cvzf Allen_MTG_Subsamples_subclass.tar.gz Allen_MTG_Subsamples_subclass
rclone copy Allen_MTG_Subsamples_subclass.tar.gz fraser:astarr/CellTypeProp_New_And_Pseudobulk/Allen_MTG
