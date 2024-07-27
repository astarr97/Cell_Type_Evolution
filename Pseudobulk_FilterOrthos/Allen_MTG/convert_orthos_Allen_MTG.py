import sys
import pandas as pd
import numpy as np
import os

#Folders containing subsamplings to be filtered for orthologs and renamed
folders = sys.argv[1].split(",")
#Base folder that contains one file from each species to get the list of genes to keep
base = sys.argv[2]
genes_keep = []
orthos = pd.read_csv("Orthologs_AllenSestan_HumChpGorRheMarm.txt", sep = "\t")
orthos = orthos[orthos["Gene name"].isin([x if type(x) is str and "MT-" not in x else "RBFOX3" for x in list(orthos["Gene name"])])]
orthos = orthos[["Gene name", "Chimpanzee gene name", "Chimpanzee homology type", "Gorilla gene name", "Gorilla homology type", "Macaque gene name", "Macaque homology type", "White-tufted-ear marmoset gene name", "White-tufted-ear marmoset homology type"]]
orthos = orthos.dropna()
orthos_121 = orthos[(orthos["Chimpanzee homology type"] == "ortholog_one2one") & (orthos["Gorilla homology type"] == "ortholog_one2one") & (orthos["Macaque homology type"] == "ortholog_one2one") & (orthos["White-tufted-ear marmoset homology type"] == "ortholog_one2one")]

#NEED TO FINISH WRITING THIS TO FIX GENE NAMES!
genes_keep = []
for file in os.listdir(base):
    if "121" not in file:
        gene = pd.read_csv(base + "/" + file, sep = "\t").set_index("features")
        if "Chimp" in file:
            orthos_121_new = orthos_121.copy().set_index("Chimpanzee gene name")
            x = orthos_121_new.join(gene).dropna()
            genes = x["Gene name"]
        elif "Rhesus" in file:
            orthos_121_new = orthos_121.copy().set_index("Macaque gene name")
            x = orthos_121_new.join(gene).dropna()
            genes = x["Gene name"]
        elif "Marmoset" in file:
            orthos_121_new = orthos_121.copy().set_index("White-tufted-ear marmoset gene name")
            x = orthos_121_new.join(gene).dropna()
            genes = x["Gene name"]
        elif "Gorilla" in file:
            orthos_121_new = orthos_121.copy().set_index("Gorilla gene name")
            x = orthos_121_new.join(gene).dropna()
            genes = x["Gene name"]
        else:
            orthos_121_new = orthos_121.copy().set_index("Gene name")
            x = orthos_121_new.join(gene).dropna()
            genes = x.index
        if len(genes_keep) == 0:
            genes_keep = genes
        else:
            genes_keep = np.intersect1d(genes_keep, genes)

orthos_121 = orthos_121[orthos_121["Gene name"].isin(genes_keep)].drop_duplicates("Gene name")

def filt_mtg(v, species):
    if species != "Human":
        orthos_121_new = orthos_121.copy().set_index(species + " gene name")
        v = v.join(orthos_121_new).dropna()
        v["features"] = v["Gene name"]
    else:
        orthos_121_new = orthos_121.copy().set_index("Gene name")
        v = v.join(orthos_121_new).dropna()
        v["features"] = v.index
    
    v = v[np.setdiff1d(v.columns, orthos_121_new.columns)]
    v = v.set_index("features")
    return v
    

#Need to rename them now!
for folder in folders:
    for file in os.listdir(folder):
        if "121" not in file:
            v = pd.read_csv(folder + "/" + file, sep = "\t")
            v.index = v["features"]
            if "Chimp" in file:
                v = filt_mtg(v, "Chimpanzee").sort_values("features")
                v.to_csv(folder + "/" + file.replace(".txt", "_All121.txt"), sep = "\t", index = True)
            elif "Rhesus" in file:
                v = filt_mtg(v, "Macaque").sort_values("features")
                v.to_csv(folder + "/" + file.replace(".txt", "_All121.txt"), sep = "\t", index = True)
            elif "Marmoset" in file:
                v = filt_mtg(v, "White-tufted-ear marmoset").sort_values("features")
                v.to_csv(folder + "/" + file.replace(".txt", "_All121.txt"), sep = "\t", index = True)
            elif "Gorilla" in file:
                v = filt_mtg(v, "Gorilla").sort_values("features")
                v.to_csv(folder + "/" + file.replace(".txt", "_All121.txt"), sep = "\t", index = True)
            else:
                v = filt_mtg(v, "Human").sort_values("features")
                v.to_csv(folder + "/" + file.replace(".txt", "_All121.txt"), sep = "\t", index = True)


    