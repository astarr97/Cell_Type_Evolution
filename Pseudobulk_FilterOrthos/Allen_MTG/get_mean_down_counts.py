import pandas as pd
import numpy as np


inhib = ['Sst Chodl_1', 'Sst_8', 'Sst_7', 'Sst_9', 'Lamp5_Lhx6_1', 'Pax6_1', 'Vip_3', 'Vip_2', 'Sst_6', 'Sst_1', 'Chandelier_1', 'Sst_5', 'Vip_7', 'Sst_4', 'Sst_2', 'Sncg_1', 'Pvalb_2', 'Pax6_2', 'Vip_1', 'Vip_5', 'Lamp5_2', 'Sncg_2', 'Sst_3', 'Vip_6', 'Pvalb_3', 'Sncg_3', 'Pvalb_4', 'Pvalb_1', 'Vip_4', 'Lamp5_1', 'Vip_8']
print(len(inhib))
x = []
for i in range(1, 101):
    v = pd.read_csv("Allen_MTG_Subsamples_cross_species_cluster/Human_Allen_MTG_SampSize_50_Round" + str(i) + "_All121.txt", sep = "\t").set_index("features")
    v2 = pd.read_csv("Allen_MTG_Subsamples_cross_species_cluster/Marmoset_Allen_MTG_SampSize_50_Round" + str(i) + "_All121.txt", sep = "\t").set_index("features")
    v = v[inhib]
    v2 = v2[inhib]
    x.append(np.min([np.min(np.sum(v, axis = 0)), np.min(np.sum(v2, axis = 0))]))
print(np.mean(x))