import scanpy as sc
import pandas as pd
import scipy
import copy
import numpy as np
import sys
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

species = sys.argv[1]
group_by = sys.argv[2]

d_caps = {"chimp":"Chimp", "gorilla":"Gorilla", "human":"Human", "rhesus":"Rhesus", "marmoset":"Marmoset", "chimpanzee":"Chimpanzee"}
l = sc.read_mtx("snRNA-seq_" + d_caps[species] + "_counts.mtx.gz")
p = pd.read_csv("snRNA-seq_" + d_caps[species] + "_cell_meta.txt", sep = "\t")
g = pd.read_csv("snRNA-seq_" + d_caps[species] + "_genes.txt", sep = "\t", header = None)
l.var = p
l.obs = g
l.var_names_make_unique()
l.var.index = l.var["cell_name"]

#Takes approximately 8 hours to run
j = copy.deepcopy(l.var)

sample_names = list(set(l.var["samplename"]))
print(sample_names)
dfs = list(np.repeat(0, len(sample_names)))
indicator = 1
c = 0
total = len(list(l.obs.index))
for index, row in l.obs.iterrows():
    if c >= 0:
        l_g = l[l.obs[0].isin([row[0]])]
        counts = pd.DataFrame(scipy.sparse.csr_matrix.todense(l_g.X)).T
        counts.index = j.index
        counts = j.join(counts)
        for i in range(len(sample_names)):
            sample_name = sample_names[i]
            sample = counts[counts["samplename"].isin([sample_name])]
            df = pd.DataFrame(sample.groupby(by=group_by)[0].sum())
            df.columns = [row[0]]
            if indicator:
                dfs[i] = df
            else:
                dfs[i] = dfs[i].join(df)
        indicator = 0
        c += 1
for i in range(len(dfs)):
    dfff = dfs[i]
    sample_name = sample_names[i]
    dfff.to_csv("Sestan_DLPFC_Pseudobulk_" + group_by + "/Pseudobulked_Sestan_2022_DLPFC_" + group_by[0].upper() + group_by[1:] + "_" + d_caps[species] + "_" + sample_name + "_" + str(c) + ".txt", sep = "\t")