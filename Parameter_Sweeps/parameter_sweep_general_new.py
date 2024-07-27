import sys
import numpy as np
import pandas as pd
import os
from collections import Counter
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics import pairwise_distances

#For SFARI, high and med are redundant and are both SFARI genes where as low is non-SFARI genes
#For GeneOrganizer, high is brain GeneOrganizer genes, med is genes in GeneOrganizer not associated with brain phenotypes, and low is all non-brain phenotype genes (includes genes in and not in GeneOrganizer)

#Parse arguments and get set up
resolution = sys.argv[1]
all_same = sys.argv[2]
down_counts = sys.argv[3]
data_source = sys.argv[4]
add_on = sys.argv[5]
print(add_on)
if add_on == "Pritchard" or add_on == "Tau" or add_on == "SFARI" or add_on == "GeneOrganizer" or add_on == "ExprLevel":
    try:
        control_bool = int(sys.argv[6])
    except:
        control_bool = 0
elif add_on == "ExprLevel":
    control_bool = 0

if add_on == "0":
    out_prefix = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_General"
    dist_out = open("Vanilla/" + out_prefix + "_distances.txt", 'w')
    trout = open("Vanilla/" + out_prefix + ".txt", 'w')
    cti = ["All", "Exc", "Inh"]
    cells = [500, 50, 100, 250]
    log = ["NoLog", "Log2"]
    filt_strats = [5, 1, 25, 50, 0, 10]
    filt_cpm = [1, 0, 5]
    dist_out.write("\t".join(["Comparison", "Iteration", "Number cells", "Raw filtering", "CPM filtering", "Log", "All genes same?", "Downsample counts?", "Cell type", "Distance spearman", "Distance pearson", "Distance euclidean", "Distance l1"]) + "\n")
    trout.write("\t".join(["Comparison", "Number cell types", "Iteration", "Number cells", "Raw filtering", "CPM filtering", "Log", "All genes same?", "Downsample counts?", "Cell type filtering", "Spearman correlation dist_met_spearman", "Spearman correlation p-value dist_met_spearman", "Spearman correlation dist_met_pearson", "Spearman correlation p-value dist_met_pearson", "Spearman correlation dist_met_euclidean", "Spearman correlation p-value dist_met_euclidean", "Spearman correlation dist_met_l1", "Spearman correlation p-value dist_met_l1", "Pearson correlation dist_met_Pearson", "Pearson correlation p-value dist_met_Pearson", "Pearson correlation dist_met_pearson", "Pearson correlation p-value dist_met_pearson", "Pearson correlation dist_met_euclidean", "Pearson correlation p-value dist_met_euclidean", "Pearson correlation dist_met_l1", "Pearson correlation p-value dist_met_l1"]) + "\n")
elif add_on == "SingleGene":
    cti = ["All", "Exc", "Inh"]
    if resolution in ["subtype", "cross_species_cluster", "cross_species_cluster_label"]:
        cells = [50]
    else:
        cells = [250]
    log = ["Log2"]
    filt_strats = [25]
    filt_cpm = [1]
elif add_on == "GeneSet":
    cti = ["All", "Exc", "Inh"]
    if resolution in ["subtype", "cross_species_cluster", "cross_species_cluster_label"]:
        cells = [50]
    else:
        cells = [250]
    log = ["Log2"]
    filt_strats = [25]
    filt_cpm = [1]
    out_prefix = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_General_" + add_on
    trout = open("GeneSet/" + out_prefix + ".txt", 'w')
    trout.write("\t".join(["Comparison", "Number cell types", "Iteration", "Number cells", "Raw filtering", "CPM filtering", "Log", "All genes same?", "Downsample counts?", "Cell type filtering", add_on, "Num Genes", "Spearman correlation dist_met_spearman", "Spearman correlation p-value dist_met_spearman", "Spearman correlation dist_met_pearson", "Spearman correlation p-value dist_met_pearson", "Spearman correlation dist_met_euclidean", "Spearman correlation p-value dist_met_euclidean", "Spearman correlation dist_met_l1", "Spearman correlation p-value dist_met_l1", "Pearson correlation dist_met_Pearson", "Pearson correlation p-value dist_met_Pearson", "Pearson correlation dist_met_pearson", "Pearson correlation p-value dist_met_pearson", "Pearson correlation dist_met_euclidean", "Pearson correlation p-value dist_met_euclidean", "Pearson correlation dist_met_l1", "Pearson correlation p-value dist_met_l1", "Spearman correlation dist_met_Spearman perms", "Spearman correlation dist_met_Pearson perms"]) + "\n")
else:
    if add_on != "Pritchard" and add_on != "Tau" and add_on != "SFARI" and add_on != "GeneOrganizer":
        if control_bool:
            out_prefix = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_General_ContTau_" + add_on
        else:
            out_prefix = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_General_" + add_on
    else:
        if control_bool:
            out_prefix = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_General_ContExp_" + add_on
        else:
            out_prefix = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_General_NotContExp_" + add_on
    dist_out = open(add_on + "/" + out_prefix + "_distances.txt", 'w')
    trout = open(add_on + "/" + out_prefix + ".txt", 'w')
    cti = ["All", "Exc", "Inh"]
    cells = [500, 50, 100, 250]
    log = ["Log2"]
    filt_strats = [10, 25]
    filt_cpm = [1, 5]
    dist_out.write("\t".join(["Comparison", "Iteration", "Number cells", "Raw filtering", "CPM filtering", "Log", "All genes same?", "Downsample counts?", "Cell type", add_on, "Distance spearman", "Distance pearson", "Distance euclidean", "Distance l1"]) + "\n")
    trout.write("\t".join(["Comparison", "Number cell types", "Iteration", "Number cells", "Raw filtering", "CPM filtering", "Log", "All genes same?", "Downsample counts?", "Cell type filtering", add_on, "Spearman correlation dist_met_spearman", "Spearman correlation p-value dist_met_spearman", "Spearman correlation dist_met_pearson", "Spearman correlation p-value dist_met_pearson", "Spearman correlation dist_met_euclidean", "Spearman correlation p-value dist_met_euclidean", "Spearman correlation dist_met_l1", "Spearman correlation p-value dist_met_l1", "Pearson correlation dist_met_Pearson", "Pearson correlation p-value dist_met_Pearson", "Pearson correlation dist_met_pearson", "Pearson correlation p-value dist_met_pearson", "Pearson correlation dist_met_euclidean", "Pearson correlation p-value dist_met_euclidean", "Pearson correlation dist_met_l1", "Pearson correlation p-value dist_met_l1"]) + "\n")


#Define needed functions

#Computes counts per million
def cpm(x):
    return x/np.sum(x)*1000000

#Computes tau, a measure of cell type-specificity    
def tau(x):
    max_x = max(x)
    n = len(x)
    return np.sum((1-x/max_x)/(n-1))

#Filters files based on cpm and raw counts
def filt(sp1_cpm, sp2_cpm, sp1, sp2, val_raw, val_cpm, i, t = "Or"):
    if t == "Or":
        filt_statement = ((sp1[i] > val_raw) | (sp2[i] > val_raw)) & ((sp1_cpm[i] > val_cpm) | (sp2_cpm[i] > val_cpm))
        sp1_to_dist = sp1_cpm[filt_statement][i]
        sp2_to_dist = sp2_cpm[filt_statement][i]
    elif t == "And":
        sp1_to_dist = sp1_cpm[(sp1[i] > val) & (sp2[i] > val)][i]
        sp2_to_dist = sp2_cpm[(sp1[i] > val) & (sp2[i] > val)][i]
    return sp1_to_dist, sp2_to_dist

#Downsample counts
def downsample_counts(x, min_counts, iteration):
    prob = x/np.sum(x)
    np.random.seed(iteration)
    return np.random.multinomial(n=min_counts, pvals = prob, size = 1)[0]

#Control for some variable, typically expression level
#Does so by taking all items in x and y with the same value and an equal number of items with 0 < log2(xi/yi) < 0.05 and 0 > log2(xi/yi) > -0.05 for xi in x and items yi in y
def control(x, y, iteration, cut=0.05):
    np.random.seed(iteration)
    #Compute log fold-change
    dif = np.log2(x[:, None] / y)
    
    #Get the positive, negative, and zero log fold-changes that are less than the cutoff
    pos = np.where((dif < cut) & (dif > 0))
    neg = np.where((dif > -cut) & (dif < 0))
    zero = np.where((dif == 0))
    
    #Randomly shuffle the dataframes then drop duplicates
    pos = pd.DataFrame(list(zip(pos[0], pos[1]))).sample(frac=1, replace = False).drop_duplicates(0).drop_duplicates(1)
    pos["PosNegZero"] = np.repeat("Pos", pos.shape[0])
    pos.index = range(0, pos.shape[0])
    neg = pd.DataFrame(list(zip(neg[0], neg[1]))).sample(frac=1, replace = False).drop_duplicates(0).drop_duplicates(1)
    neg["PosNegZero"] = np.repeat("Neg", neg.shape[0])
    neg.index = range(pos.shape[0], pos.shape[0] + neg.shape[0])
    zero = pd.DataFrame(list(zip(zero[0], zero[1]))).sample(frac=1, replace = False).drop_duplicates(0).drop_duplicates(1)
    zero["PosNegZero"] = np.repeat("Zero", zero.shape[0])
    zero.index = range(pos.shape[0] + neg.shape[0], pos.shape[0] + neg.shape[0] + zero.shape[0])
    
    #Concatenate everything together then make new dataframes for pos, neg, and zero
    df = pd.concat([pos, neg, zero]).sample(frac=1, replace = False).drop_duplicates(0).drop_duplicates(1)
    df_pos = df[df["PosNegZero"].isin(["Pos"])].copy()
    df_neg = df[df["PosNegZero"].isin(["Neg"])].copy()
    df_zero = df[df["PosNegZero"].isin(["Zero"])].copy()
    
    #Downsample so that there are the same number of pos and same number of neg genes in there.
    if df_pos.shape[0] > df_neg.shape[0]:
        df_pos = df_pos.sample(frac=df_neg.shape[0]/df_pos.shape[0], replace = False)
    elif df_neg.shape[0] > df_pos.shape[0]:
        df_neg = df_neg.sample(frac=df_pos.shape[0]/df_neg.shape[0], replace = False)
    df = pd.concat([df_pos, df_neg, df_zero])
    df[0] = df[0].astype(int)
    df[1] = df[1].astype(int)
    #Return the result
    return df

    
#Define functions to read in metadata
def read_meta_sestan_dlpfc(resolution=resolution):
    meta_hum = pd.read_csv("../Sestan_DLPFC_Processing/snRNA-seq_Human_cell_meta.txt", sep = "\t")
    meta_hum = meta_hum[~meta_hum[resolution].isin(exclude)]
    hum_ctp = Counter(meta_hum[resolution])
    
    meta_rhe = pd.read_csv("../Sestan_DLPFC_Processing/snRNA-seq_Rhesus_cell_meta.txt", sep = "\t")
    meta_rhe = meta_rhe[~meta_rhe[resolution].isin(exclude)]
    rhe_ctp = Counter(meta_rhe[resolution])
    
    meta_chp = pd.read_csv("../Sestan_DLPFC_Processing/snRNA-seq_Chimpanzee_cell_meta.txt", sep = "\t")
    meta_chp = meta_chp[~meta_chp[resolution].isin(exclude)]
    chp_ctp = Counter(meta_chp[resolution])
    
    meta_marm = pd.read_csv("../Sestan_DLPFC_Processing/snRNA-seq_Marmoset_cell_meta.txt", sep = "\t")
    meta_marm = meta_marm[~meta_marm[resolution].isin(exclude)]
    marm_ctp = Counter(meta_marm[resolution])

    #Make a dataframe of cell counts
    out = []
    for key in hum_ctp.keys():
        out.append([key, hum_ctp[key], rhe_ctp[key], chp_ctp[key], marm_ctp[key]])
    df = pd.DataFrame(out)
    df.columns = ["Cell type", "Count Human", "Count Rhesus", "Count Chimpanzee", "Count Marmoset"]
    #df.to_csv("Sestan_DLPFC_" + resolution + "_CellTypeNumbers.txt", sep = "\t")
    return df
    
def read_meta_allen_m1(resolution=resolution):
    meta_hum = pd.read_csv("../Allen_M1_Processing/Human_M1_10xV3_Metadata.txt", sep = "\t")
    meta_hum = meta_hum[~meta_hum[resolution].isin(exclude)]
    hum_ctp = Counter(meta_hum[resolution])
    
    meta_mouse = pd.read_csv("../Allen_M1_Processing/Mouse_M1_10xV3_Metadata.txt", sep = "\t")
    meta_mouse = meta_mouse[~meta_mouse[resolution].isin(exclude)]
    mouse_ctp = Counter(meta_mouse[resolution])
    
    meta_marm = pd.read_csv("../Allen_M1_Processing/Marmoset_M1_10xV3_Metadata.txt", sep = "\t")
    meta_marm = meta_marm[~meta_marm[resolution].isin(exclude)]
    marm_ctp = Counter(meta_marm[resolution])
    
    #Make a dataframe of cell counts
    out = []
    for key in hum_ctp.keys():
        out.append([key, hum_ctp[key], mouse_ctp[key], marm_ctp[key]])
    df = pd.DataFrame(out)
    df.columns = ["Cell type", "Count Human", "Count Mouse", "Count Marmoset"]
    #df.to_csv("Allen_M1_" + resolution + "_CellTypeNumbers.txt", sep = "\t")
    return df
    
def read_meta_allen_mtg(resolution=resolution):
    meta_hum = pd.read_csv("../Allen_MTG_Processing/Human_Master_metadata.txt", sep = "\t")
    meta_hum = meta_hum[meta_hum["tech"].isin(["10Xv3"])]
    meta_hum = meta_hum[~meta_hum[resolution].isin(exclude)]
    hum_ctp = Counter(meta_hum[resolution])

    meta_rhe = pd.read_csv("../Allen_MTG_Processing/Rhesus_Master_metadata.txt", sep = "\t")
    meta_rhe = meta_rhe[~meta_rhe[resolution].isin(exclude)]
    rhe_ctp = Counter(meta_rhe[resolution])

    meta_chp = pd.read_csv("../Allen_MTG_Processing/Chimp_Master_metadata.txt", sep = "\t")
    meta_chp = meta_chp[meta_chp["tech"].isin(["10Xv3"])]
    meta_chp = meta_chp[~meta_chp[resolution].isin(exclude)]
    chp_ctp = Counter(meta_chp[resolution])

    meta_marm = pd.read_csv("../Allen_MTG_Processing/Marmoset_Master_metadata.txt", sep = "\t")
    meta_marm = meta_marm[~meta_marm[resolution].isin(exclude)]
    marm_ctp = Counter(meta_marm[resolution])
    
    meta_gor = pd.read_csv("../Allen_MTG_Processing/Gorilla_Master_metadata.txt", sep = "\t")
    meta_gor = meta_gor[meta_gor["tech"].isin(["10Xv3"])]
    meta_gor = meta_gor[~meta_gor[resolution].isin(exclude)]
    gor_ctp = Counter(meta_gor[resolution])
    
    #Make a dataframe of cell counts
    out = []
    for key in hum_ctp.keys():
        out.append([key, hum_ctp[key], rhe_ctp[key], chp_ctp[key], marm_ctp[key], gor_ctp[key]])
    df = pd.DataFrame(out)
    df.columns = ["Cell type", "Count Human", "Count Rhesus", "Count Chimp", "Count Marmoset", "Count Gorilla"]
    #df.to_csv("Allen_MTG_" + resolution + "_CellTypeNumbers.txt", sep = "\t")
    return df

#Define species and cell types to keep or exclude
#Cell types were only excluded if they were non-neuronal or not found in every species
if data_source == "Sestan_DLPFC":
    species = ["Human", "Rhesus", "Chimpanzee", "Marmoset"]
    if resolution == "subclass":
        exclude = ["RB", "PC", "OPC", "Astro", "Micro", "Endo", "Immune", "Oligo", "SMC", "VLMC"]
        excit = ["L3-5 IT-1", "L2-3 IT", "L3-5 IT-3", "L6 CT", "L6 IT-1", "L5-6 NP", "L6 IT-2", "L3-5 IT-2", "L6B", "L5 ET"]
    elif resolution == "subtype":
        #Additionally exclude two neuronal subtypes not found in all species 
        #InN LAMP5 SYT10
        #L2-3 CUX2 ARHGAP18
        exclude = ['InN LAMP5 SYT10', 'L2-3 CUX2 ARHGAP18', 'Astro AQP4 OSMR', 'Astro AQP4 SLC1A2', 'Astro GFAP AQP1', 'Astro GFAP FABP7', 'B EBF1 IGKC', 'COP GPR17 SOX4', 'Endo CLDN5 DKK2', 'Endo CLDN5 IL1R1', 'Endo CLDN5 SLC7A5', 'Macro F13A1 COLEC12', 'Micro P2RY12 APBB1IP', 'Micro P2RY12 CCL3', 'Micro P2RY12 GLDN', 'Myeloid LSP1 LYZ', 'OPC PDGFRA PCDH15', 'Oligo MOG CDH7', 'Oligo MOG FRY', 'Oligo MOG GSN', 'Oligo MOG OPALIN', 'PC P2RY14 GRM8', 'RB HBA1 HBB', 'SMC ACTA2 CNN1', 'SMC ACTA2 CRISPLD2', 'SMC ACTA2 CYP1B1', 'T SKAP1 CD247', 'VLMC COL1A2 SLC13A3', 'VLMC COL1A2 VLMCA8', 'VLMC SLC47A1 SLC4A4']
        inhib = ['InN ADARB2 CADPS2', 'InN ADARB2 CALCR', 'InN ADARB2 COL5A2', 'InN ADARB2 FAM19A1', 'InN ADARB2 PALMD', 'InN ADARB2 RGS10', 'InN ADARB2 SFRP2', 'InN LAMP5 ADAMTS20', 'InN LAMP5 BMP7', 'InN LAMP5 LHX6 PROX1', 'InN LAMP5 LHX6 TAC1', 'InN LAMP5 MEIS2', 'InN LAMP5 RELN', 'InN LHX6 HGF STON2', 'InN PVALB ANOS1', 'InN PVALB GRIN2C', 'InN PVALB HPSE', 'InN PVALB PDE3A', 'InN PVALB PIEZO2', 'InN PVALB PRKCH', 'InN PVALB SNTB1', 'InN PVALB SYT2', 'InN SST ADAMTSL1', 'InN SST ADRA1D', 'InN SST ANKRD33B', 'InN SST FREM1', 'InN SST HGF GABRQ', 'InN SST HS3ST5', 'InN SST HTR2C', 'InN SST KLHL14', 'InN SST NPY', 'InN SST PLPP4', 'InN SST STK32A', 'InN SST THSD7B', 'InN SST TRPC7', 'InN VIP ADAM12', 'InN VIP CLSTN2', 'InN VIP EGF', 'InN VIP EXPH5', 'InN VIP FREM2', 'InN VIP HS3ST3B1', 'InN VIP MDGA1', 'InN VIP MEGF11 LHFP', 'InN VIP PENK', 'InN VIP SCML4']
        excit = ['L2-3 CUX2 ACVR1C INHBA', 'L2-3 CUX2 ACVR1C THSD7A', 'L2-3 CUX2 NTNG1 COL5A2', 'L2-3 CUX2 NTNG1 PALMD', 'L2-3 CUX2 NTNG1 PLCH1', 'L2-3 CUX2 PLCXD3', 'L2-4 CUX2 RORB CLMN', 'L3-5 RORB GABRG1 KCNH7', 'L3-5 RORB GABRG1 PLCH1', 'L3-5 RORB GABRG1 PTPRO', 'L3-5 RORB MKX DCC', 'L3-5 RORB MKX GALR1', 'L3-5 RORB MKX GRIN3A', 'L3-5 RORB MKX SCN5A', 'L3-5 RORB PCBP3 ARHGAP15', 'L3-5 RORB PCBP3 IL1RAPL2', 'L3-5 RORB PCBP3 LINGO2', 'L3-5 RORB TNNT2 PRRX1', 'L3-5 RORB TNNT2 TSHZ2', 'L5 FEZF2 BCL11B POU3F1', 'L5-6 FEZF2 NXPH2 CDH8', 'L5-6 FEZF2 NXPH2 GPC5', 'L5-6 FEZF2 NXPH2 GRIK1', 'L5-6 FEZF2 NXPH2 SORCS3', 'L6 FEZF2 AMOTL1 CDH6', 'L6 FEZF2 DCBLD1 MITF', 'L6 FEZF2 EYA1', 'L6 FEZF2 HCRTR2 THSD7B', 'L6 FEZF2 NPFFR2 KCNK2', 'L6 FEZF2 NPFFR2 TPD52L1', 'L6 FEZF2 SYT6 CDH9', 'L6 FEZF2 SYT6 ERBB4', 'L6 FEZF2 SYT6 PTPRT', 'L6 FEZF2 SYT6 SULF1', 'L6 OPRK1 SMYD1 ADAMTS17', 'L6 OPRK1 SMYD1 KCND2', 'L6 OPRK1 SMYD1 SNTB1', 'L6 OPRK1 SMYD1 TSHZ2', 'L6 OPRK1 THEMIS RGS6']
elif data_source == "Allen_M1":
    species = ["Human", "Mouse", "Marmoset"]
    if resolution == "subclass_label":
        exclude = ["Oligo", "Astro", "Endo", "Micro-PVM", "exclude", "OPC", "VLMC", "Peri", "SMC", "Meis2"]
        excit = ["L5 IT", "L6 CT", "L6 IT", "L6b", "L5/6 NP", "L2/3 IT", "L5 ET", "L6 IT Car3"]
    elif resolution == "cross_species_cluster_label":
        exclude = ["Oligo_1", "Astro_2", "Astro_1", "Astro", "Endo", "Microglia/PVM", "Oligo_2", "exclude", "OPC", "VLMC", "SMC"]
        excit = ["L5 IT_3", "L6 IT_1", "L5 ET_1", "L6 IT_2", "L6 CT_2", "L5 ET_2", "L5 IT_1", "L6 IT_3", "L6 CT_1", "L6b", "L5 IT_2", "L5/6 NP", "L2/3 IT"]
elif data_source == "Allen_MTG":
    species = ["Human", "Rhesus", "Gorilla", "Chimp", "Marmoset"]
    if resolution == "subclass":
        exclude = ["Astro", "VLMC", "Oligo", "Micro-PVM", "OPC", "Endo"]
        excit = ['L6 IT', 'L5 IT', 'L5/6 NP', 'L2/3 IT', 'L4 IT', 'L6 CT', 'L6b', 'L6 IT Car3', 'L5 ET']
        inhib = ['Lamp5', 'Chandelier', 'Sst', 'Pvalb', 'Sncg', 'Vip', 'Sst Chodl', 'Lamp5_Lhx6', 'Pax6']
    elif resolution == "cross_species_cluster":
        exclude = ['Astro_1', 'Endo_1', 'OPC_1', 'Oligo_1', 'OPC_2', 'Micro-PVM_1', 'VLMC_1']
        excit = ['L6 IT Car3_2', 'L6 IT Car3_1', 'L5 IT_1', 'L6 IT_1', 'L6b_1', 'L6b_3', 'L6 CT_1', 'L6 CT_2', 'L5 ET_2', 'L4 IT_1', 'L4 IT_2', 'L5 ET_1', 'L5/6 NP_2', 'L2/3 IT_2', 'L6b_2', 'L2/3 IT_1', 'L5 IT_2', 'L5/6 NP_1', 'L2/3 IT_3']
        inhib = ['Sst Chodl_1', 'Sst_8', 'Sst_7', 'Sst_9', 'Lamp5_Lhx6_1', 'Pax6_1', 'Vip_3', 'Vip_2', 'Sst_6', 'Sst_1', 'Chandelier_1', 'Sst_5', 'Vip_7', 'Sst_4', 'Sst_2', 'Sncg_1', 'Pvalb_2', 'Pax6_2', 'Vip_1', 'Vip_5', 'Lamp5_2', 'Sncg_2', 'Sst_3', 'Vip_6', 'Pvalb_3', 'Sncg_3', 'Pvalb_4', 'Pvalb_1', 'Vip_4', 'Lamp5_1', 'Vip_8']

#Depending on the add on to stratify by, read in the appropriate file
if add_on == "Pritchard":
    cut1 = 0.1
    cut2 = 0.01
    col_name = "post_mean"
    genes_to_filt = pd.read_csv("Pritchard_Constraint_Metric_GeneName.txt", sep = "\t")
    genes_high = genes_to_filt[genes_to_filt[col_name] >= cut1]["Gene name"]
    genes_med = genes_to_filt[(genes_to_filt[col_name] > cut2) & (genes_to_filt[col_name] < cut1)]["Gene name"]
    genes_low = genes_to_filt[genes_to_filt[col_name] <= cut2]["Gene name"]
elif add_on == "SFARI":
    genes_to_filt = pd.read_csv("SFARI-Gene_genes_03-28-2024release_05-09-2024export.csv", sep = ",")
    genes_aut = genes_to_filt["gene-symbol"]
elif add_on == "GeneOrganizer":
    genes_to_filt = pd.read_csv("GeneORGANizer-Confident-all-v12-BrainOnly.txt", sep = "\t")
    genes_brain = genes_to_filt[genes_to_filt["brain"] == 1]["Symbol"]
    genes_not_brain = genes_to_filt[genes_to_filt["brain"] == 0]["Symbol"]

#Function to compute distances between cell type-specific expression profiles for a single cell type from two species
def compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration, control_bool=0, genes_high=0, genes_med=0, genes_low=0, tau_df=0):
    if add_on == "0" or add_on == "GeneSet":
        x = filt(sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, i)
        sp1_to_dist = x[0]
        sp2_to_dist = x[1]
    
        if logg == "Log2":
            sp1_to_dist = np.log2(sp1_to_dist + 1)
            sp2_to_dist = np.log2(sp2_to_dist + 1)
    
        dist_sp = 1-spearmanr(sp1_to_dist, sp2_to_dist)[0]
        dist_pe = 1-pearsonr(sp1_to_dist, sp2_to_dist)[0]
        dist_euc = pairwise_distances(X=np.array(sp1_to_dist).reshape(1, -1), Y=np.array(sp2_to_dist).reshape(1, -1), metric="euclidean")[0][0]
        dist_l1 = pairwise_distances(X=np.array(sp1_to_dist).reshape(1, -1), Y=np.array(sp2_to_dist).reshape(1, -1), metric="l1")[0][0]
        if add_on == "0":
            return dist_sp, dist_pe, dist_euc, dist_l1
        else:
            return dist_sp, dist_pe, dist_euc, dist_l1, len(sp1_to_dist)
    else:
        x = filt(sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, i)
        sp1_to_dist = pd.DataFrame(x[0])
        sp2_to_dist = pd.DataFrame(x[1])
        
        sp1_to_dist_high = sp1_to_dist.loc[np.intersect1d(sp1_to_dist.index, genes_high)].copy()
        sp1_to_dist_med = sp1_to_dist.loc[np.intersect1d(sp1_to_dist.index, genes_med)].copy()
        sp1_to_dist_low = sp1_to_dist.loc[np.intersect1d(sp1_to_dist.index, genes_low)].copy()
        
        sp2_to_dist_high = sp2_to_dist.loc[np.intersect1d(sp2_to_dist.index, genes_high)].copy()
        sp2_to_dist_med = sp2_to_dist.loc[np.intersect1d(sp2_to_dist.index, genes_med)].copy()
        sp2_to_dist_low = sp2_to_dist.loc[np.intersect1d(sp2_to_dist.index, genes_low)].copy()
        keep_num = min([len(sp1_to_dist_high.index), len(sp1_to_dist_med.index), len(sp1_to_dist_low.index)])
        
        if not control_bool:
            #This version does not control for expression level
            #Randomly subsample to the correct number of genes
            keep_high = np.random.choice(sp1_to_dist_high.index, size=keep_num, replace = False)
            keep_low = np.random.choice(sp1_to_dist_low.index, size=keep_num, replace = False)
            
            #Create input arrays
            sp1_to_dist_low_high = sp1_to_dist_high.loc[keep_high].copy().sort_index()[i]
            sp1_to_dist_low_low = sp1_to_dist_low.loc[keep_low].copy().sort_index()[i]
            sp2_to_dist_low_high = sp2_to_dist_high.loc[keep_high].copy().sort_index()[i]
            sp2_to_dist_low_low = sp2_to_dist_low.loc[keep_low].copy().sort_index()[i]
            
            nf_keep = sp1_to_dist.sample(frac=keep_num/sp1_to_dist.shape[0], replace = False).index
            sp1_to_dist_nf_low = sp1_to_dist.loc[nf_keep].copy().sort_index()[i]
            sp2_to_dist_nf_low = sp2_to_dist.loc[nf_keep].copy().sort_index()[i]
            
            #Repeat for medium filtering level
            keep_med = np.random.choice(sp1_to_dist_med.index, size=keep_num, replace = False)
            sp1_to_dist_med_high = sp1_to_dist_high.loc[keep_high].copy().sort_index()[i]
            sp1_to_dist_med_med = sp1_to_dist_med.loc[keep_med].copy().sort_index()[i]
            sp2_to_dist_med_high = sp2_to_dist_high.loc[keep_high].copy().sort_index()[i]
            sp2_to_dist_med_med = sp2_to_dist_med.loc[keep_med].copy().sort_index()[i]
            
            sp1_to_dist_nf_med = sp1_to_dist.loc[nf_keep].copy().sort_index()[i]
            sp2_to_dist_nf_med = sp2_to_dist.loc[nf_keep].copy().sort_index()[i]
        elif control_bool and add_on == "ExprLevel":
            high_arr = np.array(tau_df.loc[sp1_to_dist_high.index]["Tau"])
            low_arr = np.array(tau_df.loc[sp1_to_dist_low.index]["Tau"])
            low_df = control(high_arr, low_arr, iteration, cut=0.01)
            #print(Counter(low_df["PosNegZero"]))
            keep_low_high = sp1_to_dist_high.index[list(low_df[0])]
            print(len(keep_low_high))
            keep_low_low = sp1_to_dist_low.index[list(low_df[1])]
                
            sp1_to_dist_low_high = sp1_to_dist_high.loc[keep_low_high].copy().sort_index()[i]
            sp1_to_dist_low_low = sp1_to_dist_low.loc[keep_low_low].copy().sort_index()[i]
            
            sp2_to_dist_low_high = sp2_to_dist_high.loc[keep_low_high].copy().sort_index()[i]
            sp2_to_dist_low_low = sp2_to_dist_low.loc[keep_low_low].copy().sort_index()[i]
            
            #Downsample to create a comparable version without stratifying
            nf_keep = sp1_to_dist.sample(frac=low_df.shape[0]/sp1_to_dist.shape[0], replace = False).index
            sp1_to_dist_nf_low = sp1_to_dist.loc[nf_keep].copy().sort_index()[i]
            sp2_to_dist_nf_low = sp2_to_dist.loc[nf_keep].copy().sort_index()[i]
            
            med_arr = np.array(tau_df.loc[sp1_to_dist_med.index]["Tau"])
            
            med_df = control(high_arr, med_arr, iteration, cut=0.01)
            keep_med_high = sp1_to_dist_high.index[list(med_df[0])]
            keep_med_med = sp1_to_dist_med.index[list(med_df[1])]
            sp1_to_dist_med_high = sp1_to_dist_high.loc[keep_med_high].copy().sort_index()[i]
            sp1_to_dist_med_med = sp1_to_dist_med.loc[keep_med_med].copy().sort_index()[i]
            sp2_to_dist_med_high = sp2_to_dist_high.loc[keep_med_high].copy().sort_index()[i]
            sp2_to_dist_med_med = sp2_to_dist_med.loc[keep_med_med].copy().sort_index()[i]
            
            nf_keep = sp1_to_dist.sample(frac=med_df.shape[0]/sp1_to_dist.shape[0], replace = False).index
            sp1_to_dist_nf_med = sp1_to_dist.loc[nf_keep].copy().sort_index()[i]
            sp2_to_dist_nf_med = sp2_to_dist.loc[nf_keep].copy().sort_index()[i]
        else:
            #Get summed CPM
            sp1_high = sp1_to_dist_high + sp2_to_dist_high
            sp1_low = sp1_to_dist_low + sp2_to_dist_low
            sp1_high = sp1_high.sample(frac=1, replace = False)
            sp1_low = sp1_low.sample(frac=1, replace = False)
                
            #Define two arrays
            #No need to add a pseudocount as we are filtering for genes that have non-zero CPM earlier
            high_arr = np.array(sp1_high[i])
            low_arr = np.array(sp1_low[i])
            
            #Call function to control for expression
            low_df = control(high_arr, low_arr, iteration, cut=0.05)
            
            keep_low_high = sp1_high.index[list(low_df[0])]
            keep_low_low = sp1_low.index[list(low_df[1])]
                
            sp1_to_dist_low_high = sp1_to_dist_high.loc[keep_low_high].copy().sort_index()[i]
            sp1_to_dist_low_low = sp1_to_dist_low.loc[keep_low_low].copy().sort_index()[i]
            
            sp2_to_dist_low_high = sp2_to_dist_high.loc[keep_low_high].copy().sort_index()[i]
            sp2_to_dist_low_low = sp2_to_dist_low.loc[keep_low_low].copy().sort_index()[i]
            
            #Downsample to create a comparable version without stratifying
            nf_keep = sp1_to_dist.sample(frac=low_df.shape[0]/sp1_to_dist.shape[0], replace = False).index
            sp1_to_dist_nf_low = sp1_to_dist.loc[nf_keep].copy().sort_index()[i]
            sp2_to_dist_nf_low = sp2_to_dist.loc[nf_keep].copy().sort_index()[i]
            
            sp1_med = sp1_to_dist_med + sp2_to_dist_med
            sp1_med = sp1_med.sample(frac=1, replace = False)
            med_arr = np.array(sp1_med[i])
            med_df = control(high_arr, med_arr, iteration, cut=0.05)
            keep_med_high = sp1_high.index[med_df[0]]
            keep_med_med = sp1_med.index[med_df[1]]
            sp1_to_dist_med_high = sp1_to_dist_high.loc[keep_med_high].copy().sort_index()[i]
            sp1_to_dist_med_med = sp1_to_dist_med.loc[keep_med_med].copy().sort_index()[i]
            sp2_to_dist_med_high = sp2_to_dist_high.loc[keep_med_high].copy().sort_index()[i]
            sp2_to_dist_med_med = sp2_to_dist_med.loc[keep_med_med].copy().sort_index()[i]
            
            nf_keep = sp1_to_dist.sample(frac=med_df.shape[0]/sp1_to_dist.shape[0], replace = False).index
            sp1_to_dist_nf_med = sp1_to_dist.loc[nf_keep].copy().sort_index()[i]
            sp2_to_dist_nf_med = sp2_to_dist.loc[nf_keep].copy().sort_index()[i]
        #print(i)
        #print(np.mean(sp1_to_dist_med_high), np.mean(sp1_to_dist_med_med), "Med")
        #print(np.mean(sp1_to_dist_low_high), np.mean(sp1_to_dist_low_low), "Low")
        #print(np.mean(sp2_to_dist_med_high), np.mean(sp2_to_dist_med_med), "Med")
        #print(np.mean(sp2_to_dist_low_high), np.mean(sp2_to_dist_low_low), "Low")
        if logg == "Log2":
            sp1_to_dist_nf_low = np.log2(sp1_to_dist_nf_low + 1)
            sp2_to_dist_nf_low = np.log2(sp2_to_dist_nf_low + 1)
            sp1_to_dist_nf_med = np.log2(sp1_to_dist_nf_med + 1)
            sp2_to_dist_nf_med = np.log2(sp2_to_dist_nf_med + 1)
            
            sp1_to_dist_low_high = np.log2(sp1_to_dist_low_high + 1)
            sp1_to_dist_low_low = np.log2(sp1_to_dist_low_low + 1)
            sp2_to_dist_low_high = np.log2(sp2_to_dist_low_high + 1)
            sp2_to_dist_low_low = np.log2(sp2_to_dist_low_low + 1)

            sp2_to_dist_med_high = np.log2(sp2_to_dist_med_high + 1)
            sp2_to_dist_med_med = np.log2(sp2_to_dist_med_med + 1)
            sp1_to_dist_med_high = np.log2(sp1_to_dist_med_high + 1)
            sp1_to_dist_med_med = np.log2(sp1_to_dist_med_med + 1)
        #print(len(sp1_to_dist_nf_low), len(sp1_to_dist_low_low))
        to_app = [i]
        types = ["Med High", "Med Med", "Low High", "Low Low", "NoFilt Low Size", "NoFilt Med Size"]
        tups = [(sp1_to_dist_med_high, sp2_to_dist_med_high), (sp1_to_dist_med_med, sp2_to_dist_med_med), (sp1_to_dist_low_high, sp2_to_dist_low_high), (sp1_to_dist_low_low, sp2_to_dist_low_low), (sp1_to_dist_nf_low, sp2_to_dist_nf_low), (sp1_to_dist_nf_med, sp2_to_dist_nf_med)]
        if control_bool == 0:
            types = ["Med High", "Med Med", "Low Low", "No Filt Low Size"]
            d_types = {"Med High":"High", "Med Med":"Med", "Low Low":"Low", "No Filt Low Size":"No Filt Matched Size"}
            tups = [(sp1_to_dist_med_high, sp2_to_dist_med_high), (sp1_to_dist_med_med, sp2_to_dist_med_med), (sp1_to_dist_low_low, sp2_to_dist_low_low), (sp1_to_dist_nf_low, sp2_to_dist_nf_low)]
        for ind in range(len(tups)):
            tup = tups[ind]
            df_tup1 = pd.DataFrame(tup[0].copy())
            df_tup1.columns = ["Tup1"]
            df_tup2 = pd.DataFrame(tup[1].copy())
            df_tup2.columns = ["Tup2"]
            df_tup = df_tup1.join(df_tup2)

            dist_sp = 1-spearmanr(tup[0], tup[1])[0]
            dist_pe = 1-pearsonr(tup[0], tup[1])[0]
            dist_euc = pairwise_distances(X=np.array(tup[0]).reshape(1, -1), Y=np.array(tup[1]).reshape(1, -1), metric="euclidean")[0][0]
            dist_l1 = pairwise_distances(X=np.array(tup[0]).reshape(1, -1), Y=np.array(tup[1]).reshape(1, -1), metric="l1")[0][0]
            to_app = to_app + [dist_sp, dist_pe, dist_euc, dist_l1]
            if control_bool == 0:
                dist_out.write("\t".join([str(elem) for elem in [species1 + "-" + species2, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, i, d_types[types[ind]]] + [dist_sp, dist_pe, dist_euc, dist_l1]]) + "\n")
                print("\t".join([str(elem) for elem in [species1 + "-" + species2, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, i, d_types[types[ind]]] + [dist_sp, dist_pe, dist_euc, dist_l1]]) + "\n")
            else:
                dist_out.write("\t".join([str(elem) for elem in [species1 + "-" + species2, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, i, types[ind]] + [dist_sp, dist_pe, dist_euc, dist_l1]]) + "\n")
                #print("\t".join([str(elem) for elem in [species1 + "-" + species2, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, i, types[ind]] + [dist_sp, dist_pe, dist_euc, dist_l1]]) + "\n")
        return to_app
    
#Function to compute correlations with cell type proportion        
def compute_corrs(dist_prop, dist1 = "Spearman", dist2 = "Pearson", dist3 = "Euclidean", dist4 = "L1"):
    corr_sp = spearmanr(np.log10(dist_prop["Proportion use"]), dist_prop[dist1])
    corr_pe = spearmanr(np.log10(dist_prop["Proportion use"]), dist_prop[dist2])
    corr_euc = spearmanr(np.log10(dist_prop["Proportion use"]), dist_prop[dist3])
    corr_l1 = spearmanr(np.log10(dist_prop["Proportion use"]), dist_prop[dist4])
    pcorr_sp = pearsonr(np.log10(dist_prop["Proportion use"]), dist_prop[dist1])
    pcorr_pe = pearsonr(np.log10(dist_prop["Proportion use"]), dist_prop[dist2])
    pcorr_euc = pearsonr(np.log10(dist_prop["Proportion use"]), dist_prop[dist3])
    pcorr_l1 = pearsonr(np.log10(dist_prop["Proportion use"]), dist_prop[dist4])
    
    return corr_sp, corr_pe, corr_euc, corr_l1, pcorr_sp, pcorr_pe, pcorr_euc, pcorr_l1

#Function to finish an dwrite out    
def finish(out, df_new, control_bool, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, trout):
    dists = pd.DataFrame(out)
    cols = ["Cell type"]
    types = ["Med High", "Med Med", "Low High", "Low Low", "NoFilt Low Size", "NoFilt Med Size"]
    #If control_bool is 0, then there is only one "High" estimate and the "Low" and "Med" sizes are identical
    if control_bool == 0:
        types = ["High", "Med", "Low", "No Filt Matched Size"]
    
    for t in types:
        cols = cols + [dist_metric + " " + t for dist_metric in ["Spearman", "Pearson", "Euclidean", "L1"]]
    dists.columns = cols
    dists = dists.set_index("Cell type")
    dist_prop = df_new.join(dists)
   
    for t in types:
        keep_cols = [dist_metric + " " + t for dist_metric in ["Spearman", "Pearson", "Euclidean", "L1"]]
        dist_prop_k = dist_prop[list(df_new.columns) + keep_cols].copy()
        corr_sp, corr_pe, corr_euc, corr_l1, pcorr_sp, pcorr_pe, pcorr_euc, pcorr_l1 = compute_corrs(dist_prop_k, keep_cols[0], keep_cols[1], keep_cols[2], keep_cols[3])
        #print("\t".join([str(elem) for elem in [species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, t, corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1]]]))
        trout.write("\t".join([str(elem) for elem in [species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, t, corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1]]]) + "\n")
        #print([species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, t, corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1]])

#For 100 iterations, go through all parameter combinations
for iteration in range(1, 101):
    print(iteration, resolution)
    for species1 in species:
        for species2 in species:
            if species1 != species2 and species.index(species1) < species.index(species2):
                if data_source == "Sestan_DLPFC":
                    df = read_meta_sestan_dlpfc()
                elif data_source == "Allen_M1":
                    df = read_meta_allen_m1()
                elif data_source == "Allen_MTG":
                    df = read_meta_allen_mtg()
                
                for num_cells in cells:
                    if data_source == "Sestan_DLPFC":
                        sp1 = pd.read_csv("../Sestan_DLPFC_Processing/Sestan_DLPFC_Subsamples_" + resolution + "/" + species1 + "_Sestan_DLPFC_SampSize_" + str(num_cells) + "_Round" + str(iteration) + "_All121.txt", sep= "\t").set_index("0")
                        sp2 = pd.read_csv("../Sestan_DLPFC_Processing/Sestan_DLPFC_Subsamples_" + resolution + "/" + species2 + "_Sestan_DLPFC_SampSize_" + str(num_cells) + "_Round" + str(iteration) + "_All121.txt", sep= "\t").set_index("0")
                    elif data_source == "Allen_M1":
                        sp1 = pd.read_csv("../Allen_M1_Processing/Allen_M1_Subsamples_" + resolution + "/" + species1 + "_Allen_M1_SampSize_" + str(num_cells) + "_Round" + str(iteration) + "_All121.txt", sep= "\t").set_index("features")
                        sp2 = pd.read_csv("../Allen_M1_Processing/Allen_M1_Subsamples_" + resolution + "/" + species2 + "_Allen_M1_SampSize_" + str(num_cells) + "_Round" + str(iteration) + "_All121.txt", sep= "\t").set_index("features")
                    elif data_source == "Allen_MTG":
                        sp1 = pd.read_csv("../Allen_MTG_Processing/Allen_MTG_Subsamples_" + resolution + "/" + species1 + "_Allen_MTG_SampSize_" + str(num_cells) + "_Round" + str(iteration) + "_All121.txt", sep= "\t").set_index("features")
                        sp2 = pd.read_csv("../Allen_MTG_Processing/Allen_MTG_Subsamples_" + resolution + "/" + species2 + "_Allen_MTG_SampSize_" + str(num_cells) + "_Round" + str(iteration) + "_All121.txt", sep= "\t").set_index("features")

                    sp1 = sp1[np.setdiff1d(sp1.columns, exclude)]
                    sp2 = sp2[np.setdiff1d(sp2.columns, exclude)]
                    for ct_keep in cti:
                        for logg in log:
                            for val_raw in filt_strats:
                                for val_cpm in filt_cpm:
                                    #Filter by num cells
                                    if data_source == "Sestan_DLPFC":
                                        df_new = df[(df["Count Human"] >= num_cells) & (df["Count Rhesus"] >= num_cells) & (df["Count Chimpanzee"] >= num_cells) & (df["Count Marmoset"] >= num_cells)]
                                    elif data_source == "Allen_M1":
                                        df_new = df[(df["Count Human"] >= num_cells) & (df["Count Mouse"] >= num_cells) & (df["Count Marmoset"] >= num_cells)]
                                    elif data_source == "Allen_MTG":
                                        df_new = df[(df["Count Human"] >= num_cells) & (df["Count Rhesus"] >= num_cells) & (df["Count Chimp"] >= num_cells) & (df["Count Marmoset"] >= num_cells) & (df["Count Gorilla"] >= num_cells)]
                                    #Filter by cell type
                                    df_new = df_new[~df_new["Cell type"].isin(exclude)].copy()
                                    if ct_keep == "Exc":
                                        sp1_new = sp1[excit].copy()
                                        sp2_new = sp2[excit].copy()
                                        df_new = df_new[df_new["Cell type"].isin(excit)].copy()
                                    elif ct_keep == "Inh":
                                        sp1_new = sp1[np.setdiff1d(sp1.columns, excit)].copy()
                                        sp2_new = sp2[np.setdiff1d(sp2.columns, excit)].copy()
                                        df_new = df_new[~df_new["Cell type"].isin(excit)].copy()
                                    else:
                                        sp1_new = sp1.copy()
                                        sp2_new = sp2.copy()

                                    sp1_new = sp1_new[df_new["Cell type"]].copy()
                                    sp2_new = sp2_new[df_new["Cell type"]].copy()
                                    df_new = df_new.set_index("Cell type")
                                    
                                    for spec in species:
                                        df_new["Proportion " + spec] = df_new["Count " + spec]/np.sum(df_new["Count " + spec])
                                    df_new["Proportion use"] = (df_new["Proportion " + species1] + df_new["Proportion " + species2])/2
                                    
                                    #Downsample counts if needed
                                    if down_counts == "DownCounts":
                                        lowest_counts = min(list(np.sum(sp1_new, axis = 0)) + list(np.sum(sp2_new, axis = 0)))
                                        sp1_new = sp1_new.apply(downsample_counts, axis = 0, min_counts = lowest_counts, iteration = iteration)
                                        sp2_new = sp2_new.apply(downsample_counts, axis = 0, min_counts = lowest_counts, iteration = iteration)

                                    sp1_cpm = sp1_new.apply(cpm, axis = 0)
                                    sp2_cpm = sp2_new.apply(cpm, axis = 0)

                                    assert(list(sp1_cpm.index) == list(sp2_cpm.index))
                                    assert(list(sp1_cpm.columns) == list(sp2_cpm.columns))
                                    
                                    #Filter to the same set of genes across all cell types if needed
                                    if all_same == "AllSame":
                                        sp1_truth = (sp1_new > val_raw).astype(int)
                                        sp2_truth = (sp2_new > val_raw).astype(int)
                                        sp1_cpm_truth = (sp1_cpm > val_cpm).astype(int)
                                        sp2_cpm_truth = (sp2_cpm > val_cpm).astype(int)
                                        truth = np.minimum(sp1_truth+sp2_truth, 1)
                                        truth_cpm = np.minimum(sp1_cpm_truth+sp2_cpm_truth, 1)
                                        truth = np.min(truth, axis = 1)
                                        truth_cpm = np.min(truth_cpm, axis = 1)
                                        keep = np.multiply(truth, truth_cpm)
                                        keep = list(keep[keep > 0].index)
                                        sp1_new = sp1_new.loc[keep]
                                        sp2_new = sp2_new.loc[keep]
                                        sp1_cpm = sp1_cpm.loc[keep]
                                        sp2_cpm = sp2_cpm.loc[keep]
                                    
                                    if add_on == "0":
                                        #Go through and compute distances
                                        out = []
                                        for i in list(sp1_cpm.columns):
                                            dist_sp, dist_pe, dist_euc, dist_l1 = compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration)
                                            out.append([i, dist_sp, dist_pe, dist_euc, dist_l1])
                                            if ct_keep == "All":
                                                dist_out.write("\t".join([str(elem) for elem in [species1 + "-" + species2, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, i, dist_sp, dist_pe, dist_euc, dist_l1]]) + "\n")
                                        
                                        #Finish by computing everything
                                        dists = pd.DataFrame(out)
                                        dists.columns = ["Cell type", "Spearman", "Pearson", "Euclidean", "L1"]
                                        dists = dists.set_index("Cell type")
                                        dist_prop = df_new.join(dists)
                                        corr_sp, corr_pe, corr_euc, corr_l1, pcorr_sp, pcorr_pe, pcorr_euc, pcorr_l1 = compute_corrs(dist_prop)
                                        #print("\t".join([str(elem) for elem in [species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1]]]))
                                        #trout.write("\t".join([str(elem) for elem in [species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1]]]) + "\n")
                                        #print([species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1]])
                                    #A routine if we are seeing how expression level affects things
                                    elif add_on == "ExprLevel":
                                        if control_bool:
                                            out = []
                                            df_mcpm = (sp1_cpm + sp2_cpm)/2
                                            #print(df_mcpm.columns)
                                            tau_df = pd.DataFrame(df_mcpm.apply(tau, axis = 1))
                                            tau_df.columns = ["Tau"]
                                            for i in list(sp1_cpm.columns):
                                                np.random.seed(iteration)
                                                
                                                #Partition into lowly, moderately, and highly expressed genes
                                                x = filt(sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, i)
                                                sp1_to_dist = pd.DataFrame(x[0]).copy()
                                                sp1_to_dist.columns = [i + " Sp1"]
                                                sp2_to_dist = pd.DataFrame(x[1]).copy()
                                                sp2_to_dist.columns = [i + " Sp2"]
                                                sp_dist = sp1_to_dist.join(sp2_to_dist)
                                                sp_dist["Mean CPM"] = (sp_dist[i + " Sp1"] + sp_dist[i + " Sp2"])/2
                                                sp_dist = sp_dist.sort_values("Mean CPM")
                                                step_size = len(sp_dist.index)//3
                                                genes_low = list(sp_dist.index)[0:step_size]
                                                genes_med = list(sp_dist.index)[step_size:step_size*2]
                                                genes_high = list(sp_dist.index)[step_size*2:step_size*3]
                                                
                                                
    
                                                #Compute distances and append
                                                to_app = compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration, 1, genes_high, genes_med, genes_low, tau_df = tau_df)
                                                out.append(to_app)
                                            
                                            #Finish the computation and write out to trout
                                            finish(out, df_new, control_bool, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, trout)
                                        else:
                                            out = []
                                            for i in list(sp1_cpm.columns):
                                                np.random.seed(iteration)
                                                
                                                #Partition into lowly, moderately, and highly expressed genes
                                                x = filt(sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, i)
                                                sp1_to_dist = pd.DataFrame(x[0]).copy()
                                                sp1_to_dist.columns = [i + " Sp1"]
                                                sp2_to_dist = pd.DataFrame(x[1]).copy()
                                                sp2_to_dist.columns = [i + " Sp2"]
                                                sp_dist = sp1_to_dist.join(sp2_to_dist)
                                                sp_dist["Mean CPM"] = (sp_dist[i + " Sp1"] + sp_dist[i + " Sp2"])/2
                                                sp_dist = sp_dist.sort_values("Mean CPM")
                                                step_size = len(sp_dist.index)//3
                                                genes_low = list(sp_dist.index)[0:step_size]
                                                genes_med = list(sp_dist.index)[step_size:step_size*2]
                                                genes_high = list(sp_dist.index)[step_size*2:step_size*3]
    
                                                #Compute distances and append
                                                to_app = compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration, 0, genes_high, genes_med, genes_low)
                                                out.append(to_app)
                                            
                                            #Finish the computation and write out to trout
                                            finish(out, df_new, control_bool, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, trout)
                                    #A routine to see how constraint on expression affects things
                                    elif add_on == "Pritchard":
                                        out = []
                                        for i in list(sp1_cpm.columns):
                                            #Set random seed
                                            np.random.seed(iteration)
                                            to_app = compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration, control_bool, genes_high, genes_med, genes_low)
                                            out.append(to_app)
                                    
                                        finish(out, df_new, control_bool, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, trout)
                                    #Adding a routine to work with the GeneOrganizer data
                                    elif add_on == "GeneOrganizer":
                                        out = []
                                        for i in list(sp1_cpm.columns):
                                            #Set random seed
                                            np.random.seed(iteration)
                                            to_app = compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration, control_bool, genes_brain, genes_not_brain, np.setdiff1d(sp1_cpm.index, genes_brain))
                                            out.append(to_app)
                                    
                                        finish(out, df_new, control_bool, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, trout)
                                    #Adding a routine to work with the SFARI data
                                    elif add_on == "SFARI":
                                        out = []
                                        for i in list(sp1_cpm.columns):
                                            np.random.seed(iteration)
                                            to_app = compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration, control_bool, genes_aut, genes_aut, np.setdiff1d(sp1_cpm.index, genes_aut))
                                            out.append(to_app)
                                        finish(out, df_new, control_bool, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, trout)
                                        
                                    #Now checking on tissue-specifity
                                    elif add_on == "Tau":
                                        assert(list(sp1_cpm.columns) == list(sp2_cpm.columns))
                                        
                                        df_mcpm = (sp1_cpm + sp2_cpm)/2
                                        #print(df_mcpm.columns)
                                        tau_df = pd.DataFrame(df_mcpm.apply(tau, axis = 1))
                                        tau_df.columns = ["Tau"]
                                        
                                        out = []
                                        for i in list(sp1_cpm.columns):
                                            #If we are in NotAllSame mode, genes can be considered to be in one tau category for one cell type but a different tau category for another cell type
                                            x = filt(sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, i)
                                            genes_tau = pd.DataFrame(x[0]).copy().index
                                            tau_df_i = tau_df.loc[genes_tau].copy().sort_values("Tau")
                                            out_prefix_tau = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_" + add_on + "_" + "_".join([ct_keep, str(num_cells), str(val_raw), str(val_cpm), str(iteration)])
                                            out_file_tau = add_on + "_TauValues" + "/" + out_prefix_tau + ".txt"
                                            tau_df_i.to_csv(out_file_tau, sep = "\t")
                                            step_size = len(genes_tau)//3
                                            
                                            genes_low = list(tau_df_i.index)[0:step_size]
                                            genes_med = list(tau_df_i.index)[step_size:step_size*2]
                                            genes_high = list(tau_df_i.index)[step_size*2:step_size*3]
                                            
                                            #print(np.mean(tau_df_i.loc[genes_low]["Tau"]), np.mean(tau_df_i.loc[genes_med]["Tau"]), np.mean(tau_df_i.loc[genes_high]["Tau"]))
                                            #Set random seed
                                            np.random.seed(iteration)
                                            to_app = compute_dist(i, sp1_cpm, sp2_cpm, sp1_new, sp2_new, val_raw, val_cpm, add_on, iteration, control_bool, genes_high, genes_med, genes_low)
                                            out.append(to_app)
                                    
                                        finish(out, df_new, control_bool, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, trout)
                                    #Routine to compute single gene correlations
                                    elif add_on == "SingleGene":
                                        out_prefix = data_source + "_" + resolution + "_" + all_same + "_" + down_counts + "1-100_" + add_on + "_" + "_".join([ct_keep, str(num_cells), str(val_raw), str(val_cpm), str(iteration)])
                                        out_file = add_on + "/" + out_prefix + ".txt"
                                        assert(all_same == "AllSame")
                                        #Compute the log2 fold-change
                                        sp_lfc = np.abs(np.log2((sp1_cpm + 1)/(sp2_cpm + 1))).T.sort_index()
                                        sp_lfc_df = sp_lfc.join(df_new).copy().sort_index()
                                        
                                        ### NEED TO CHANGE THIS TO LOG10 AND RERUN IF WE USE PEARSON ###
                                        spear = sp_lfc.apply(spearmanr, b=np.log2(sp_lfc_df["Proportion use"]))
                                        pear = sp_lfc.apply(pearsonr, y=np.log2(sp_lfc_df["Proportion use"]))
                                        spear.index = ["Spearman rho", "Spearman p-value"]
                                        pear.index = ["Pearson rho", "Pearson p-value"]
                                        df_out = pd.concat([spear, pear]).T.dropna()
                                        df_out.to_csv(out_file, sep = "\t")
                                    #Routine to do computations for gene sets
                                    elif add_on == "GeneSet":
                                        gene_sets = ["GO_Biological_Process_2023", "Human_Phenotype_Ontology"]
                                        for gene_set in gene_sets:
                                            map_gene_sets = pd.read_csv(gene_set + "_GSEAPY_Efficient_GEQ100.txt", sep = "\t", header = None)
                                            map_gene_sets.columns = ["Category", "Genes"]
                                            map_gene_sets["Genes"] = [x.split(",") for x in list(map_gene_sets["Genes"])]
                                            c = 0
                                            
                                            for index, row in map_gene_sets.iterrows():
                                                genes_keep = row["Genes"]
                                                sp1_cpm_gk = sp1_cpm.loc[np.intersect1d(sp1_cpm.index, genes_keep)].copy().sort_index()
                                                
                                                if len(sp1_cpm_gk) >= 100:
                                                    try:
                                                        c += 1
                                                        sp2_cpm_gk = sp2_cpm.loc[np.intersect1d(sp2_cpm.index, genes_keep)].copy().sort_index()
                                                        sp1_new_gk = sp1_new.loc[np.intersect1d(sp1_new.index, genes_keep)].copy().sort_index()
                                                        sp2_new_gk = sp2_new.loc[np.intersect1d(sp2_new.index, genes_keep)].copy().sort_index()
                                                        out = []
                                                        sizes = []
                                                        for i in list(sp1_cpm_gk.columns):
                                                            dist_sp, dist_pe, dist_euc, dist_l1, size = compute_dist(i, sp1_cpm_gk, sp2_cpm_gk, sp1_new_gk, sp2_new_gk, val_raw, val_cpm, add_on, iteration)
                                                            out.append([i, dist_sp, dist_pe, dist_euc, dist_l1])
                                                            sizes.append(size)
                                                        
                                                        #Finish by computing everything
                                                        dists = pd.DataFrame(out)
                                                        dists.columns = ["Cell type", "Spearman", "Pearson", "Euclidean", "L1"]
                                                        dists = dists.set_index("Cell type")
                                                        dist_prop = df_new.join(dists)
                                                        corr_sp, corr_pe, corr_euc, corr_l1, pcorr_sp, pcorr_pe, pcorr_euc, pcorr_l1 = compute_corrs(dist_prop)
                                                        
                                                        #Do 100 matched permutations to see if we could get the correlation for that gene set by chance
                                                        perm_sp = []
                                                        perm_pear = []
                                                        size = len(sp1_cpm_gk)
                                                        for _ in range(100):
                                                            out = []
                                                            #Pick a random subset of genes
                                                            keep_genes_perm = np.random.choice(list(sp1_cpm.index), size, replace=False)
                                                            sp1_cpm_perm = sp1_cpm.loc[keep_genes_perm].copy().sort_index()
                                                            sp2_cpm_perm = sp2_cpm.loc[keep_genes_perm].copy().sort_index()
                                                            sp1_new_perm = sp1_new.loc[keep_genes_perm].copy().sort_index()
                                                            sp2_new_perm = sp2_new.loc[keep_genes_perm].copy().sort_index()
                                                            
                                                            #Go through and compute dists for these genes
                                                            for i in list(sp1_cpm.columns):
                                                                dist_sp, dist_pe, dist_euc, dist_l1, size_not_used = compute_dist(i, sp1_cpm_perm, sp2_cpm_perm, sp1_new_perm, sp2_new_perm, val_raw, val_cpm, add_on, iteration)
                                                                out.append([i, dist_sp, dist_pe, dist_euc, dist_l1])
                                                                dists = pd.DataFrame(out)
                                                            dists.columns = ["Cell type", "Spearman", "Pearson", "Euclidean", "L1"]
                                                            dists = dists.set_index("Cell type")
                                                            dist_prop = df_new.join(dists)
                                                            corr_sp, corr_pe, corr_euc, corr_l1, pcorr_sp, pcorr_pe, pcorr_euc, pcorr_l1 = compute_corrs(dist_prop)
                                                            perm_sp.append(str(corr_sp[0]))
                                                            perm_pear.append(str(corr_pe[0]))
                                                        perm_sp = ",".join(perm_sp)
                                                        perm_pear = ",".join(perm_pear)
                                                        trout.write("\t".join([str(elem) for elem in [species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, row["Category"], len(sp1_cpm_gk.index), corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1], perm_sp, perm_pear]]) + "\n")
                                                        #print([str(elem) for elem in [species1 + "-" + species2, dist_prop.shape, iteration, num_cells, val_raw, val_cpm, logg, all_same, down_counts, ct_keep, row["Category"], len(sp1_cpm_gk.index), corr_sp[0], corr_sp[1], corr_pe[0], corr_pe[1], corr_euc[0], corr_euc[1], corr_l1[0], corr_l1[1], pcorr_sp[0], pcorr_sp[1], pcorr_pe[0], pcorr_pe[1], pcorr_euc[0], pcorr_euc[1], pcorr_l1[0], pcorr_l1[1], perm_sp, perm_pear]])
                                                    except:
                                                        pass
if add_on != "SingleGene":
    trout.close()
    dist_out.close()
