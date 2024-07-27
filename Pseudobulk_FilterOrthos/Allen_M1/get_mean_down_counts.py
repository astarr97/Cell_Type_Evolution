import pandas as pd
import numpy as np

exclude = ["Oligo_1", "Astro_2", "Astro_1", "Astro", "Endo", "Microglia/PVM", "Oligo_2", "exclude", "OPC", "VLMC", "SMC"]
excit = ["L5 IT_3", "L6 IT_1", "L5 ET_1", "L6 IT_2", "L6 CT_2", "L5 ET_2", "L5 IT_1", "L6 IT_3", "L6 CT_1", "L6b", "L5 IT_2", "L5/6 NP", "L2/3 IT"]
x = []
for i in range(1, 101):
    v = pd.read_csv("Allen_M1_Subsamples_cross_species_cluster_label/Human_Allen_M1_SampSize_50_Round" + str(i) + "_All121.txt", sep = "\t").set_index("features")
    v2 = pd.read_csv("Allen_M1_Subsamples_cross_species_cluster_label/Marmoset_Allen_M1_SampSize_50_Round" + str(i) + "_All121.txt", sep = "\t").set_index("features")
    v = v[np.setdiff1d(v.columns, exclude + excit)]
    v2 = v2[np.setdiff1d(v2.columns, exclude + excit)]
    print(len(list(v2.columns)))
    x.append(np.min([np.min(np.sum(v, axis = 0)), np.min(np.sum(v2, axis = 0))]))
print(np.mean(x))
