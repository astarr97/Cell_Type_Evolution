Each folder contains code for the corresponding dataset.
process_for_ctp.sh is the primary shell script that contains the commands needed to downsample cells and pseudobulk.
process_for_ctp.py was used to downsample cells and pseudobulk.
convert_orthos_<data_set>.py as used to restrict to one-to-one non-mitochondrial orthologs
For Allen_MTG and Sestan_DLPFC, there is additional code to pseudobulk without doing any downsampling.
The resulting files were used for differential expression analysis.
