import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import time
start_time = time.time()
# decoration of program
parser = argparse.ArgumentParser(
    description='Creating rna-dna contacts distribution')
parser.add_argument('--all', action='store_true', help='analyze all genome')
parser.add_argument('--start_chr', action='store', dest='start_chr',
                    help='input start chromosome number')
parser.add_argument('--start_position', action='store', default=1,
                    dest='start_point', help='input start position', type=int)
parser.add_argument('--end_chr', action='store', dest='end_chr',
                    help='input end chromosome number')
parser.add_argument('--end_position', action='store', dest='end_point',
                    default='end', help='input end position')
parser.add_argument('-i', action='store', dest='input_filename',
                    help='input filename with RNA-DNA contacts')
parser.add_argument('-o', action='store', dest='output_filename',
                    help='output filename')
parser.add_argument('-b', action='store', dest='my_bin', required=True,
                    help='type bin size in nt', type=int)
res = parser.parse_args()
input_filename = res.input_filename
output_filename = res.output_filename
start_chr = res.start_chr
end_chr = res.end_chr
start_point = res.start_point
end_point = res.end_point
my_bin = res.my_bin

df = pd.read_csv(input_filename, sep="\t", header=None)
# read file with dna-rna contacts
lengths_df = pd.read_csv("length_chrom", sep="\t", header=None)
# read file containing lengths of chromosome

if res.all:
    start_chr = 'chr1'
    end_chr = 'chr25'
    start_point = 1
    end_point = 16569

elif res.end_point == 'end':
    if end_chr == 'chrX':
        ending = 23
    else:
        ending = int(end_chr[3:])
    end_point = int(lengths_df[1][ending - 1])

print('creating bins')

end_point = int(end_point)
chr_set = []
chr_set.append(start_chr[3:])
chr_set.append(end_chr[3:])
for i in range(2):
    if chr_set[i] == 'X':
        chr_set[i] = 23
    elif chr_set[i] == 'Y':
        chr_set[i] = 24
    elif chr_set[i] == 'M':
        chr_set[i] = 25
    else:
        chr_set[i] = int(chr_set[i])

heat_df = pd.DataFrame(columns=['chr', 'position', 'bin_number'])
bin_dict = {}
if chr_set[0] == chr_set[1]:
    for i in range(start_point, end_point, my_bin):
        heat_df.loc[i] = ['chr'+str(chr_set[0]), int(i), len(heat_df)]
    bin_dict[chr_set[0]] = len(heat_df)-1
    
else:
        # first chromosome input:
    for i in range(start_point, int(lengths_df[1][chr_set[0] - 1]), my_bin):
        heat_df.loc[len(heat_df)] = ['chr'+str(chr_set[0]), int(i), len(heat_df)]

    bin_dict[chr_set[0]] = len(heat_df)-1

    # core        
    for chr_num in range(chr_set[0] + 1, chr_set[1]):
        for i in range(1, int(lengths_df[1][chr_num-1]), my_bin):
            heat_df.loc[len(heat_df)] = ['chr'+str(chr_num), int(i), len(heat_df)]
        bin_dict[chr_num] = len(heat_df) - 1
    # last chromosome
    for i in range(1, end_point, my_bin):
            heat_df.loc[len(heat_df)] = ['chr'+str(chr_set[1]), int(i), len(heat_df)]
    bin_dict[chr_set[1]] = len(heat_df) - 1
           
    # insert previous chromosome (required for algoritm)    
    bin_dict[chr_set[0]-1] = -1

heat_df.to_csv(output_filename+'bins', encoding='utf-8', index=False)
print('done')
print('preparing df')


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

df[4] = df[4].apply(chromosome_extract)
df = df.sort_values(by=[4, 5])
df_for_analyse = df.iloc[:, 4:]
del df_for_analyse[6]

df_for_analyse = df_for_analyse.values

# determing contacts
final_df_for_heatmap = np.empty((0, 2), int)
for i in range(len(df_for_analyse)):
    key_for_bin = df_for_analyse[i, 0] - 1
    dna_bin = bin_dict[key_for_bin] + df_for_analyse[i, 1] // my_bin + 1 
    futher = np.array([[dna_bin, df_for_analyse[i, 2]]])
    final_df_for_heatmap = np.append(final_df_for_heatmap, futher, axis=0)
index_list = [i in range(len(final_df_for_heatmap))]
final_df_for_heatmap = pd.DataFrame(final_df_for_heatmap) 
final_df_for_gisto = final_df_for_heatmap.groupby([0])[1].sum()
final_df_for_gisto.to_csv(output_filename+'_contact_matrix', encoding='utf-8', index=False)  

# creating plot
bins_df = pd.Series(bin_dict, name="bins")
bins_df.index.name = 'chr'
bins_df.reset_index()
bins_list = list(bins_df.values)
bins_list = bins_list[1:]
label_list = bins_df.index.values.tolist()
label_list = label_list[1:]

a = sns.barplot(x= final_df_for_gisto.index, y=final_df_for_gisto.values, color="blue")
plt.xticks(bins_list[1:], label_list, size=20)
plt.xlabel('')
plt.ylabel('RNA contacts')
fig = a.get_figure()
fig.savefig(output_filename)
print("--- %s seconds ---" % (time.time() - start_time))
