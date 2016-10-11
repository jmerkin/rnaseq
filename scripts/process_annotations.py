#!/usr/bin/python

import subprocess
import generate_regions as GR
import sys
import os
from optparse import OptionParser

script_directory = os.path.dirname(os.path.realpath(__file__))
source_directory = script_directory.rstrip('scripts')

parser = OptionParser()
parser.add_option('-s', '--seqlen', dest='seqlen', default='650')
parser.add_option('-e', '--extra_files', dest='extra_files') 
parser.add_option('-i', '--infile', dest='infile', 
    default='%s/annotations/gencode.v24.annotation.gff3.gz' % source_directory)
parser.add_option('-o', '--outdir', dest='outdir', 
    default= '%s/processed_annotations/' % source_directory)

op, args = parser.parse_args()

infile = op.infile
seqlen = int(op.seqlen)
extra_files = [] if op.extra_files is None else op.extra_files.split(',')

outdir = os.path.join(op.outdir, infile.split('/')[-1])
if not os.path.exists(outdir):
    os.makedirs(outdir)
subprocess.Popen('ln -s %s %s' % (infile, outdir), shell=1).wait()

file_names = {
        'tx_fasta': os.path.join(outdir, 'transcript.fa'), 
        'intron_ratio_fa': os.path.join(outdir, 'exon_intron_ratio.fa'), 
        'end_bias_ratio_fa': os.path.join(outdir, 'end_bias_ratio.fa'), 
        }
if extra_files:
    file_names['full_fasta'] = os.path.join(outdir, 'complete.fa')

tx_file = open(os.path.join(outdir, 'transcripts.tsv'), 'w')
ce_file = open(os.path.join(outdir, 'const_exons.bed'), 'w')

files = {key: open(value, 'w') for key, value in file_names.iteritems()}

for name, gene in GR.generate_gene(infile):
    txes = [line for tx in gene for line in tx if line[2] == 'transcript']
    txnames = [txline[-1].split(';')[0].split('=')[1] for txline in txes]
    for txname in txnames:
        tx_file.write('%s\t%s\n' % (txname, name))
    gene = [[exon for exon in tx if exon[2] == 'exon'] for tx in gene]
    reverse = gene[0][0][6] == '-'
    chromosome = gene[0][0][0] 

    for exon in GR.generate_const_exons(gene):
        ce_file.write('%s\t%s\t%s\n' % exon)

    introns = GR.generate_introns(gene) 
    introns = list(introns)
    for seq, txname in zip(GR.generate_sequence(gene, gff=1, reverse=reverse, nested=1), txnames):
        #files['tx_fasta'].write('>%s:%s\n%s\n' % (name, txname, seq))
        #files['full_fasta'].write('>%s:%s\n%s\n' % (name, txname, seq))
        files['tx_fasta'].write('>%s\n%s\n' % (txname, seq))
        if extra_files:
            files['full_fasta'].write('>%s\n%s\n' % (txname, seq))

    long_txs = [sum(int(exon[4]) - int(exon[3]) + 1 for exon in tx) > (3 * seqlen) 
            for tx in gene]
    moderate_introns = [1000 < (int(intron[2]) - int(intron[1])) < 10000 for intron in introns]

    if any(long_txs): 
        long_seqs = GR.generate_sequence(
            [tx for tx, is_long in zip(gene, long_txs) if is_long], 
            gff=1, reverse=reverse, nested=1)
        for snum, seq in enumerate(long_seqs):
            files['end_bias_ratio_fa'].write('>%s:%s:5p\n%s\n' % (name, snum, seq[:seqlen]))
            files['end_bias_ratio_fa'].write('>%s:%s:3p\n%s\n' % (name, snum, seq[-seqlen:]))

    if any(moderate_introns):
        mod_seqs = [GR.generate_sequence([intron], gff=0, reverse=reverse, nested=0).next() 
                for intron, is_mod in zip(introns, moderate_introns) if is_mod]
        for snum, seq in enumerate(mod_seqs):
            files['intron_ratio_fa'].write('>%s:%s:intron\n%s\n' % (name, snum, seq))
        for txnum, seq in enumerate(GR.generate_sequence(gene, gff=1, reverse=reverse, nested=1)):
            files['intron_ratio_fa'].write('>%s:%s:tx\n%s\n' % (name, txnum, seq))

map(lambda xx: xx.close(), files.values())
tx_file.close()
ce_file.close()

if extra_files:
    for ff in extra_files:
        cmd = "cat %s | tr '[:lower:]' '[:upper:]' >> %s" % (ff, os.path.join(outdir, 'complete.fa'))
        subprocess.Popen(cmd, shell=1).wait()

for fname in file_names.values():
    print 'indexing', fname
    cmd = 'salmon index -i %s.index -t %s' % (fname, fname)
    subprocess.Popen(cmd, shell=1,).wait()

