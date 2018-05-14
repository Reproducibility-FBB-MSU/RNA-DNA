import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

df = pd.read_csv("dna_reads/dna_read_exon", sep="\t", header=None)
dna_exon = sum(df[7])
df = pd.read_csv("dna_reads/dna_read_intron", sep="\t", header=None)
dna_intron = sum(df[7])
df = pd.read_csv("dna_reads/dna_read_downstream", sep="\t", header=None)
dna_down = sum(df[7])
df = pd.read_csv("dna_reads/dna_read_upstream", sep="\t", header=None)
dna_up = sum(df[7])
df = pd.read_csv("dna_reads/dna_read_nongenic", sep="\t", header=None)
dna_non = sum(df[7])

# Make data: I have 2 groups and 7 subgroups
group_names = ['Genic', 'Intergenic']
group_size = [dna_exon + dna_intron, dna_non]
subgroup_names = ['Exon', 'Intron', 'Upstream', 'Downstream', 'Distal intergenic']
subgroup_size = [dna_exon, dna_intron, dna_up, dna_down, dna_non - dna_up - dna_down]
 
# Create colors
a, b, c = [plt.cm.Blues, plt.cm.YlOrRd, plt.cm.Greens]
out_circle = [c(0.8), a(0.4),'purple', 'grey', b(0.5)]

# First Ring (outside)
fig, ax = plt.subplots()
ax.axis('equal')
mypie2, _ = ax.pie(subgroup_size, radius=1.3, colors=out_circle)
plt.setp(mypie2, width=0.4, edgecolor='white')
plt.margins(0, 0)

# Second Ring (Inside)
mypie, _ = ax.pie(group_size, radius=1.3 - 0.3, colors=[a(0.4), b(0.3)], labels=group_names, labeldistance=0.4)
plt.setp(mypie, edgecolor='white')
 
# show it
plt.legend(mypie2, subgroup_names, bbox_to_anchor=(1.2,0.5), loc="center right", fontsize=10, 
           bbox_transform=plt.gcf().transFigure)

plt.title("RNA reads", size=18)
plt.show()
plt.savefig("RNA reads")   

df = pd.read_csv("rna_reads/rna_read_exon", sep="\t", header=None)
rna_exon = sum(df[7])
df = pd.read_csv("rna_reads/rna_read_intron", sep="\t", header=None)
rna_intron = sum(df[7])
df = pd.read_csv("rna_reads/rna_read_downstream", sep="\t", header=None)
rna_down = sum(df[7])
df = pd.read_csv("rna_reads/rna_read_upstream", sep="\t", header=None)
rna_up = sum(df[7])
df = pd.read_csv("rna_reads/rna_read_nongenic", sep="\t", header=None)
rna_non = sum(df[7])
# Make data: I have 2 groups and 7 subgroups
group_names = ['Genic', 'Intergenic']
group_size = [rna_exon + rna_intron, rna_non]
subgroup_names = ['Exon', 'Intron', 'Upstream', 'Downstream', 'Distal intergenic']
subgroup_size = [rna_exon, rna_intron, rna_up, rna_down, rna_non - rna_up - rna_down]
 
# Create colors
a, b, c =[plt.cm.Blues, plt.cm.YlOrRd, plt.cm.Greens]
out_circle =[c(0.6), a(0.8),'purple', 'grey', b(0.5)]

# First Ring (outside)
fig, ax = plt.subplots()
ax.axis('equal')
mypie2, _ = ax.pie(subgroup_size, radius=1.3,  colors=out_circle)
plt.setp(mypie2, width=0.4, edgecolor='white')
plt.margins(0, 0)

# Second Ring (Inside)
mypie, _ = ax.pie(group_size, radius=1.3 - 0.3, colors=[a(0.4), b(0.3)], labels=group_names, labeldistance=0.4)
plt.setp(mypie, edgecolor='white')
 

 
# show it
plt.legend(mypie2, subgroup_names, bbox_to_anchor=(1.2, 0.5), loc="center right", fontsize=10, 
           bbox_transform=plt.gcf().transFigure)

plt.title("DNA reads", size=18)
plt.savefig("DNA reads") 
