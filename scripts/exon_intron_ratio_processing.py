#!/usr/bin/env python

import pandas as pd
import sys

def single_gene(gene):
    if gene.intron.all():
        return None
    gene = gene.groupby('intron').sum()
    retval = pd.Series([gene.loc[False, 'TPM'], gene.loc[False, 'TPM'] / (1e-3 + gene.loc[True, 'TPM'])], index=['total', 'ratio'] )
    return retval 

infile = sys.argv[1]

df = pd.read_csv(infile, sep='\t', header=0, index_col=0)
df['gene'] = df.index.str.split(':').str.get(0)
df['intron'] = df.index.str.endswith('intron')

results = df.groupby('gene').apply(single_gene)
results = results.sort_values('total')
ratios = pd.np.log2(results.iloc[results.shape[0] - 1100:results.shape[0] - 100].ratio)
mean = 2 ** ratios.mean()
median = 2 ** ratios.median()

open(infile.replace('quant.sf', 'mean_intron_ratio'), 'w').write('%.2f' % mean)
open(infile.replace('quant.sf', 'median_intron_ratio'), 'w').write('%.2f' % median)

#import IPython ; IPython.embed()

