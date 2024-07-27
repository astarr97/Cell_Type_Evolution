import sys
import pandas as pd
import numpy as np
import os

#Folders containing subsamplings to be filtered for orthologs and renamed
folders = sys.argv[1].split(",")
#Base folder that contains one file from each species to get the list of genes to keep
base = sys.argv[2]
genes_keep = []
for file in os.listdir(base):
    genes = pd.read_csv(base + "/" + file, sep = "\t")["0"]
    if len(genes_keep) == 0:
        genes_keep = genes
    else:
        genes_keep = np.intersect1d(genes_keep, genes)

#Read in the orthologs and restrict to 4-way one-to-one orthologs
orthos = pd.read_csv("Orthologs_AllenSestan_HumChpGorRheMarm.txt", sep = "\t")
orthos = orthos[orthos["Gene name"].isin([x if type(x) is str and "MT-" not in x else "RBFOX3" for x in list(orthos["Gene name"])])]
orthos = orthos[["Gene name", "Chimpanzee gene name", "Chimpanzee homology type", "Macaque gene name", "Macaque homology type", "White-tufted-ear marmoset gene name", "White-tufted-ear marmoset homology type"]]
orthos = orthos.dropna()
orthos_121 = orthos[(orthos["Chimpanzee homology type"] == "ortholog_one2one") & (orthos["Macaque homology type"] == "ortholog_one2one") & (orthos["White-tufted-ear marmoset homology type"] == "ortholog_one2one")]

#Further restrict to genes in genes_keep
#Because this is Sestan, everything has human gene names so we will just always use human gene names
genes_keep = np.intersect1d(genes_keep, orthos_121["Gene name"])
print(len(genes_keep))
for folder in os.listdir():
    if folder in folders:
        for file in os.listdir(folder):
            if ".txt" in file and "All121" not in file:
                v = pd.read_csv(folder + "/" + file, sep = "\t")
                #In this case, all the gene names are identical so we will use the human orthologs
                v = v[v["0"].isin(genes_keep)]
                assert(".txt" in file)
                v.to_csv(folder + "/" + file.replace(".txt", "_All121.txt"), sep = "\t", index = False)

    