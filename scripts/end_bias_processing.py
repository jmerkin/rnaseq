#!/usr/bin/env python

import pandas as pd
import sys

def single_gene(gene):
    gene = gene.groupby('end').sum()
    retval = pd.Series([gene.TPM.sum(), gene.loc[True, 'TPM'] / (1e-3 + gene.loc[False, 'TPM'])], index=['total', 'ratio'] )
    return retval 

infile = sys.argv[1]

df = pd.read_csv(infile, sep='\t', header=0, index_col=0)
df['gene'] = df.index.str.split(':').str.get(0)
df['end'] = df.index.str.endswith('3p')

results = df.groupby('gene').apply(single_gene)
results = results.sort_values('total').dropna()
ratios = pd.np.log2(results.iloc[results.shape[0] - 1100:results.shape[0] - 100].ratio)
ratios = ratios.loc[pd.np.isfinite(ratios)]
mean = 2 ** ratios.mean()
median = 2 ** ratios.median()

#import IPython ; IPython.embed()
open(infile.replace('quant.sf', 'mean_end_ratio'), 'w').write('%.2f' % mean)
open(infile.replace('quant.sf', 'median_end_ratio'), 'w').write('%.2f' % median)


