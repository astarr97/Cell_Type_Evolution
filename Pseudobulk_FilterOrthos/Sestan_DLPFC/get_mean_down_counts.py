import pandas as pd
import numpy as np

inhib = ['InN ADARB2 CADPS2', 'InN ADARB2 CALCR', 'InN ADARB2 COL5A2', 'InN ADARB2 FAM19A1', 'InN ADARB2 PALMD', 'InN ADARB2 RGS10', 'InN ADARB2 SFRP2', 'InN LAMP5 ADAMTS20', 'InN LAMP5 BMP7', 'InN LAMP5 LHX6 PROX1', 'InN LAMP5 LHX6 TAC1', 'InN LAMP5 MEIS2', 'InN LAMP5 RELN', 'InN LHX6 HGF STON2', 'InN PVALB ANOS1', 'InN PVALB GRIN2C', 'InN PVALB HPSE', 'InN PVALB PDE3A', 'InN PVALB PIEZO2', 'InN PVALB PRKCH', 'InN PVALB SNTB1', 'InN PVALB SYT2', 'InN SST ADAMTSL1', 'InN SST ADRA1D', 'InN SST ANKRD33B', 'InN SST FREM1', 'InN SST HGF GABRQ', 'InN SST HS3ST5', 'InN SST HTR2C', 'InN SST KLHL14', 'InN SST NPY', 'InN SST PLPP4', 'InN SST STK32A', 'InN SST THSD7B', 'InN SST TRPC7', 'InN VIP ADAM12', 'InN VIP CLSTN2', 'InN VIP EGF', 'InN VIP EXPH5', 'InN VIP FREM2', 'InN VIP HS3ST3B1', 'InN VIP MDGA1', 'InN VIP MEGF11 LHFP', 'InN VIP PENK', 'InN VIP SCML4']
print(len(inhib))
x = []
for i in range(1, 101):
    v = pd.read_csv("Sestan_DLPFC_Subsamples_subtype/Human_Sestan_DLPFC_SampSize_50_Round" + str(i) + "_All121.txt", sep = "\t").set_index("0")
    v2 = pd.read_csv("Sestan_DLPFC_Subsamples_subtype/Marmoset_Sestan_DLPFC_SampSize_50_Round" + str(i) + "_All121.txt", sep = "\t").set_index("0")
    v = v[inhib]
    v2 = v2[inhib]
    x.append(np.min([np.min(np.sum(v, axis = 0)), np.min(np.sum(v2, axis = 0))]))
print(np.mean(x))
