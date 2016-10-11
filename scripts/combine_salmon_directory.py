#!/usr/bin/python

import pandas as pd
import os
import sys

test = 1
test = 0

indir = sys.argv[1]

potential_samples = os.listdir(indir)
samples = [ii for ii in os.listdir(indir) if os.path.exists('%s/%s/quant.sf' % (indir, ii))]

if test:
    samples = samples[:5]

print "combining %s runs" % len(samples)
print "excluding: %s" % ','.join([ii for ii in potential_samples if ii not in set(samples)])
files = ['%s/%s/quant.sf' % (indir, ii) for ii in samples] 
dfs = [pd.read_csv(ii, index_col=0, sep='\t') for ii in files]

genes = dfs[0].index.str.split('|').str.get(1)
txs = dfs[0].index.str.split('|').str.get(0)
df = pd.concat([result.loc[:, ['TPM']].rename(columns={'TPM': sample}) for result, sample in zip(dfs, samples)], axis=1)
df.to_csv('%s/tpm.csv' % indir)
df.index = txs
df.to_csv('%s/tpm.tsv' % indir, sep='\t')
df['gene'] = genes
df = df.groupby('gene').sum()
df.to_csv('%s/gene_tpm.csv' % indir)
df.to_csv('%s/gene_tpm.tsv' % indir, sep='\t')

df = pd.concat([result.loc[:, ['NumReads']].rename(columns={'NumReads': sample}) for result, sample in zip(dfs, samples)], axis=1)
df.to_csv('%s/reads.csv' % indir)
df['gene'] = genes
df = df.groupby('gene').sum()
df.to_csv('%s/gene_reads.csv' % indir)

df = pd.concat([result.loc[:, ['EffectiveLength']].rename(columns={'EffectiveLength': sample}) for result, sample in zip(dfs, samples)], axis=1)
df.to_csv('%s/efflength.csv' % indir)
df['gene'] = genes
df = df.groupby('gene').mean()
df.to_csv('%s/gene_efflength.csv' % indir)
