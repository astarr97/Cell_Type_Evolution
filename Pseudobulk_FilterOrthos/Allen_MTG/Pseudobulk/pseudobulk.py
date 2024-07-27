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

d_caps = {"chimp":"Chimp", "gorilla":"Gorilla", "human":"Human", "rhesus":"Rhesus", "marmoset":"Marmoset"}
l = sc.read(species + "_SCT_UMI_expression_matrix.h5ad")
l.var_names_make_unique()
v = pd.read_csv(d_caps[species] + "_Master_metadata.txt", sep = "\t")
v.index = v["sample_id"]
print(list(l.obs["sample_id"]) == list(v["sample_id"]))
l.obs = v
print(l)
l = l[l.obs["tech"].isin(["10Xv3"])].copy()
l = l.T
print(l)
#Takes approximately 8 hours to run
j = copy.deepcopy(l.var)

sample_names = list(set(l.var["donor"]))
print(sample_names)
dfs = list(np.repeat(0, len(sample_names)))
indicator = 1
c = 0
total = len(list(l.obs.index))
for index, row in l.obs.iterrows():
    if c >= 0:
        l_g = l[l.obs["features"].isin([row["features"]])]
        counts = pd.DataFrame(scipy.sparse.csr_matrix.todense(l_g.X)).T
        counts.index = j.index
        counts = j.join(counts)
        for i in range(len(sample_names)):
            sample_name = sample_names[i]
            sample = counts[counts["donor"].isin([sample_name])]
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
    dfff.to_csv("Pseudobulked_Allen_2023_MTG_" + d_caps[species] + "_" + sample_name + "_" + str(c) + ".txt", sep = "\t")