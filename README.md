# Cell_Type_Evolution
Code used to study evolution of cell types.  Each subfolder contains a README to describe the code therein.

First, we used Convert_RDS_h5ad.R to convert the RDS files to h5ad files for use with scanpy.
Next, we used the code in Pseudobulk_FilterOrthos to downsample cells across 100 iterations, pseudboulk counts across cells, and restrict to one-to-one non-mitochondrial orthologs.
Next, we used the code in Parameter_Sweeps to compare cell type-specific expression divergence and cell type proportion.
Finally, the code in Analysis_Figure_Making was used, in conjunction with powerpoint and biorender, to make figures anddo the remaining analysis based on the outputs of the above code and other publicly available data.

Note: Code to align, quantify, and do differential expression analysis for the human-chimpanzee hybrid cortical organoid data can be found at https://github.com/astarr97/Z_Score.
