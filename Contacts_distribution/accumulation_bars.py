import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import collections
import math
from matplotlib import rc
import argparse

parser = argparse.ArgumentParser(
    description='Creating cis-trans-local contact distribution')
parser.add_argument('-i', action='store', dest='input_filename',
                    help='input filename with RNA-DNA contacts')
parser.add_argument('-o', action='store', dest='output_filename',
                    help='output filename')    
res = parser.parse_args()
input_filename = res.input_filename
output_filename = res.output_filename
df = pd.read_csv(input_filename, sep="\t", header=None)


def chromosome_extract(elem):
    chr_num = elem[3:]
    if chr_num == 'X':
        elem = 23
    elif chr_num == 'Y':
        elem = 24
    elif chr_num == 'M':
        elem = 25       
    else:
        elem = int(chr_num)
    return elem

chromosome_extract = np.vectorize(chromosome_extract)
df[0] = df[0].apply(chromosome_extract)
df[4] = df[4].apply(chromosome_extract)
del df[3]
df_for_analyse = df.values


# counting cis/trans/local
def get_cis_trans_loc(array):
    dict_contact = collections.Counter()
    for i in range(len(array)):
        if array[i, 0] != array[i, 3]:
            dict_contact['trans'] += array[i, 6]
        else:
            if math.fabs(array[i, 1] - array[i, 4]) <= 10000:
                dict_contact['local'] += array[i, 6]
            else:
                dict_contact['cis'] += array[i, 6]
    return dict_contact

chromosomes = set(df_for_analyse[:, 0])

chr_contacts = pd.DataFrame(columns=['chr', 'cis', 'trans', 'local'])
for elem in chromosomes:
    local_array = df_for_analyse[df_for_analyse[:, 0] == elem]
    chr_dic = get_cis_trans_loc(local_array)
    new_line = ['chr' + str(elem), chr_dic['cis'], chr_dic['trans'], chr_dic['local']]
    chr_contacts.loc[len(chr_contacts)] = new_line
chr_contacts = chr_contacts.T
chr_contacts.to_csv(output_filename+'_bar_matrix', encoding='utf-8', index=False)       
# y-axis in bold
rc('font', weight='bold')
 
# Values of each group

bars1 = chr_contacts.iloc[1].values
bars2 = chr_contacts.iloc[2].values
bars3 = chr_contacts.iloc[3].values
 
# Heights of bars1 + bars2 (TO DO better)
bars = []
for x in range(len(bars1)):
    bars.append(bars1[x]+bars2[x])
    

# The position of the bars on the x-axis
r = [x for x in range(len(bars1))]
 
# Names of group and bar width
names = chr_contacts.iloc[0]
barWidth = 0.6

width, height = 20, 10
mpl.rcParams['figure.figsize'] = [width, height]
# Create green bars
plt.bar(r, bars1, color='g', edgecolor='black', width=barWidth, label='Cis')
# Create blue bars (middle), on top of the firs ones
plt.bar(r, bars2, bottom=bars1, color='b', edgecolor='black', width=barWidth, label='Trans')
# Create yellow bars (top)
plt.bar(r, bars3, bottom=bars, color='yellow', edgecolor='black', width=barWidth, label='Local')
 
# Custom X axis
plt.xticks(r, names, fontweight='bold', fontsize = 12)
plt.legend(fontsize=20)
plt.ylabel("Number of contacts", size=25) 
# Save graphic
plt.savefig(output_filename)
