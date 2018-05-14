import pandas as pd
input_file_name = str(input())
df = pd.read_csv(input_file_name, sep="\t", header=None)
my_df = df.iloc[:, 4:7]
my_df[7] = df.iloc[:, 3]
my_df[8] = df.iloc[:, 0]
my_df[9] = df.iloc[:, 1]
my_df[10] = df.iloc[:, 2]
my_df[11] = df.iloc[:, 7]
my_df.to_csv('rna_dna_' + input_file_name, encoding='utf-8', sep="\t", index=False)
