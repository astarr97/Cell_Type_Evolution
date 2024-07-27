This contains the code that was used for differential expression between species after pseudobulking.
The zip file Gene_Lengths contains the gene lengths for each species.
The files within must be copied to the folder in which DE analysis will be done.
In general, we used the jupyter notebook to split each pseudobulked table by cell type.
We then created config files using the next cell in the jupyter notebook.  These config files point to the files that will be used for DE analysis and include any covariates and the species being compared.  For example, the file for a cell type might look like:
Sample1_Human_L2-3IT_Pseudobulked_Counts.txt	Human
Sample2_Human_L2-3IT_Pseudobulked_Counts.txt	Human
Sample3_Human_L2-3IT_Pseudobulked_Counts.txt	Human
Sample4_Human_L2-3IT_Pseudobulked_Counts.txt	Human
Sample1_Chimp_L2-3IT_Pseudobulked_Counts.txt	Chimp
Sample2_Chimp_L2-3IT_Pseudobulked_Counts.txt	Chimp
Sample3_Chimp_L2-3IT_Pseudobulked_Counts.txt	Chimp
Sample4_Chimp_L2-3IT_Pseudobulked_Counts.txt	Chimp

We then used the R script which iterates through all the Config files and uses DESeq2 and apeglm to do differential expression analysis.

After that, we switch the sign of log fold-change for rhesus and marmoset so that they are consistent with the other files when comparing to human.
