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
cell_level = sys.argv[2]
out_folder = sys.argv[3]

d_caps = {"chimp":"Chimp", "gorilla":"Gorilla", "human":"Human", "rhesus":"Rhesus", "marmoset":"Marmoset", "chimpanzee":"Chimpanzee"}
l = sc.read_mtx("snRNA-seq_" + d_caps[species] + "_counts.mtx.gz")
p = pd.read_csv("snRNA-seq_" + d_caps[species] + "_cell_meta.txt", sep = "\t")
g = pd.read_csv("snRNA-seq_" + d_caps[species] + "_genes.txt", sep = "\t", header = None)
l.var = p
l.obs = g
l.var_names_make_unique()
l.var.index = l.var["cell_name"]
l = l.T

ind = 1
for i in range(1, 101):
    print(i)
    for j in [25, 50, 100, 250, 500, 1000]:
        dfl = 0
        for k in list(set(list(l.obs[cell_level]))):
            l_cur = l[l.obs[cell_level].isin([k])]
            cells_l = list(l_cur.obs.index)
            np.random.seed(i + j)
            np.random.shuffle(cells_l)
            keep_l = cells_l[0:j]
            l_cur = l_cur[l_cur.obs.index.isin(keep_l)]
            df_lll = pd.DataFrame(l_cur.X.todense()).T
            df_lll.index = l_cur.var[0]
            df_lll["Gene Name"] = df_lll.index
            df_lll = df_lll.sort_values("Gene Name").drop("Gene Name", axis = 1)
            #Compute centroid in gene space
            cent_l = pd.DataFrame(np.sum(df_lll, axis = 1))
            if ind:
                dfl = cent_l
                dfl.columns = [k]
                ind = 0
            else:
                dfl = dfl.join(cent_l)
                dfl.columns = list(dfl.columns)[0:-1] + [k]
        ind = 1
        dfl.to_csv(out_folder + "/" + d_caps[species] + "_Sestan_DLPFC_SampSize_" + str(j) + "_" + "Round" + str(i) + ".txt", sep = "\t")