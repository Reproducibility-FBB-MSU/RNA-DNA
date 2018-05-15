
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import collections
import seaborn
import matplotlib.pyplot as plt


GRID_file = "K562_chr18_19.tab.txt"
ench_file = "ench"
chr_number = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17","chr18", "chr19", "chr20", "chr21", "chr22", "X", "Y", "M"]

window = 10000
bin_size = 100


def chr_extract(elem):
    chr_num = elem[3:]
    if chr_num == "X":
        elem = 23
    elif chr_num == "Y":
        elem = 24
    elif chr_num == "M":
        elem = 25
    else:
        elem = int(chr_num)
    return elem



def count_info(df_win_start, df_win_stop, df_chr_start, df_chr_score, bin_size):
    
    count_final = collections.Counter({0: 0})
    
    for line in range(len(df_win_start)):
        c = collections.Counter()
        win_start = df_win_start[line]
        win_end = df_win_stop[line]
        ind = np.logical_and((df_chr_start > win_start), (df_chr_start < win_end))
        
        read_db = df_chr_start[ind]
        score_db = df_chr_score[ind]
        read_db = read_db - win_start
        
        for i in range (len(read_db)):
            a = read_db[i] // bin_size
            c[a] += score_db[i]
            count_final += c
    return count_final


f, ax = plt.subplots()
ax.set_title("Enchancer distribution", loc='center')

df_chr_all = pd.read_csv(GRID_file, sep="\t", header=None)
df_ench_all = pd.read_csv(ench_file, sep="\t", header=None, skiprows=1)

for chrn in chr_number:
    
    
    df_chr = df_chr_all[df_chr_all[4] == chrn].copy()
    chr_extract = np.vectorize(chr_extract)
    df_chr[4] = df_chr[4].apply(chr_extract)
    
    del df_chr[0]
    del df_chr[1]
    del df_chr[2]
    del df_chr[3]
    del df_chr[6]
    
    df_chr_start = df_chr[5].values
    df_chr_score = df_chr[7].values
    
    df_ench = df_ench_all[df_ench_all[1] == chrn].copy()
    df_ench[1] = df_ench[1].apply(chr_extract)
    del df_ench[0]
    del df_ench[3]
    del df_ench[4]
    
    #df_ench = df_ench[2].values
    df_ench[4] = df_ench[2] - (window/2)
    df_ench[5] = df_ench[2] + (window/2)
    
    df_ench_analyse = df_ench.values
    df_win_start = df_ench[4].values
    df_win_stop = df_ench[5].values
    
    
    #building dataframe for experimental data
    count_final_exp = count_info(df_win_start, df_win_stop, df_chr_start, df_chr_score, bin_size)
    
    bins = list(count_final_exp.keys())
    reads = list(count_final_exp.values())
    
    df_exp = pd.DataFrame({'bin_number': bins, 'reads_number': reads}).sort_values('bin_number').reset_index(drop=True)
    df_exp['reads_number'] = df_exp['reads_number'] / len(reads)
    
    #plt.title("enchancers", loc='center')
    df_exp_plot = df_exp.plot(x='bin_number', y='reads_number' ,figsize=(12,8), grid=True, label="Experimental data {}".format(chrn), ax = ax)
    #fig = df_exp_plot.get_figure()
#%matplotlib notebook
fig.savefig("all_genome.png")

#df_exp_plot = df_exp.plot(x='bin_number', y='reads_number' ,figsize=(12,8), grid=True, label="Experimental data", color="red", ax=ax)
#df_ran_plot = df_exp.plot(x='bin_number', y='ran_reads_number' ,figsize=(12,8), grid=True, label="Randomized data", color="blue", ax=ax)
#fig1 = df_exp_plot.get_figure()
#fig2 = df_ran_plot.get_figure()



