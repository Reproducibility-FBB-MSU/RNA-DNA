import numpy as np
import pandas as pd
import math

df = pd.read_csv("rna_intersected", sep="\t", header=None)
# in this file first - rna, second - dna part of contact
list_contact = []
for i in range(len(df)):
    if df.iloc[i, 0] != df.iloc[i, 4]:
        list_contact.append('trans')
    else:
        if math.fabs(int(df.iloc[i, 1]) - int(df.iloc[i, 5])) <= 10000:
            list_contact.append('local')
        else:
            list_contact.append('cis')
list_contact = pd.Series(list_contact)
df[13] = list_contact  # append type of contact

new = df.groupby([11, 13, 12])[7].sum().to_frame(name='count').reset_index()
# group by name and type of contact
new['sum'] = (new.apply(lambda x: x[11] + ', ' + x[12], 1))
del new[11]
del new[12]

df_final = new.pivot(index='sum', columns=13, values='count')
df_final['Total'] = df_final.sum(axis=1)

df_final['cis'] = round(df_final['cis']/df_final['Total']*100, 2)
df_final['local'] = round(df_final['local']/df_final['Total']*100, 2)
df_final['trans'] = round(df_final['trans']/df_final['Total']*100, 2)

del df_final['Total']
df_final.columns = ['cis', 'local', 'trans']
a = df_final.index
new_ind = list(range(len(df_final)))
df_final.index = new_ind
a = pd.Series(a)
df_final['sum'] = a
df_final['gene_id'], df_final['type'] = df_final['sum'].str.split(',').str.get(0), df_final['sum'].str.split(',').str.get(1)
del df_final['sum']
df_final.to_csv('df_triangle', encoding='utf-8', sep="\t", index=False)
