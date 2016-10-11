#!/usr/bin/python

import subprocess
import gzip
import sys
from string import maketrans
trans = maketrans('ATGCNatcgn', 'TACGNtagcn')

def process_attributes(unproc_attr):
    unproc_attr = [ii.split('=') for ii in unproc_attr.split(';')]
    attr = {ii[0]: ii[1] for ii in unproc_attr}
    return attr

def overlap(ex1, ex2):
    if (ex2[0] <= ex1[0] <= ex2[1]) or (ex2[0] <= ex1[1] <= ex2[1]):
        return True
    return False

def exon_in_tx(exon, tx):
    return any(overlap(exon, ex2) for ex2 in tx)

def generate_const_exons(gene):
    chrm = gene[0][0][0]
    gene = [[ (int(exon[3]) - 1, int(exon[4])) for exon in tx] for tx in gene]
    exons = set([exon for tx in gene for exon in tx])
    exons = [(chrm, str(exon[0]), str(exon[1])) for exon in exons if all([exon_in_tx(exon, tx) for tx in gene])]
    exons = ['\t'.join(exon) for exon in exons]
    cmd = 'bedtools sort -i <( echo -e "%s" ) | bedtools merge -i stdin ' % '\n'.join(exons)
    pp = subprocess.Popen(cmd, shell=1, executable='/bin/bash', stdout=subprocess.PIPE)
    exons = map(tuple, [ii.strip('\n').split('\t') for ii in pp.stdout])
    for exon in exons:
        yield exon

def generate_gene(infile):
    of = gzip.open if infile.endswith('gz') else open
    fo = of(infile)
    for line in fo:
        if not line.startswith('#'):
            break
    #lines = [line.strip('\n').split('\t')]

    line = line.strip('\n').split('\t')
    gene, tx = [], []
    prev_attr = process_attributes(line[8])
    gene_name = prev_attr['gene_name'] if 'gene_name' in prev_attr else None

    for line in fo:
        if line.startswith('#'):
            continue
        line = line.strip('\n').split('\t')
        attr = process_attributes(line[8])
        if prev_attr['gene_id'] == attr['gene_id']:
            if ('transcript_id' not in prev_attr):
                tx.append(line)
            elif prev_attr['transcript_id'] == attr['transcript_id']:
                tx.append(line)
            else:
                gene.append(tx)
                tx = [line]
        else:
            gene.append(tx)
            yield gene_name, gene
            gene, tx = [], []
            gene_name = attr['gene_name'] if 'gene_name' in attr else None
        prev_attr = attr

    gene.append(tx)
    yield gene_name, gene

def generate_introns(gene, gff=True):
    offset = -1 if gff else 0
    reverse = gene[0][0][6] == '-'
    #if reverse:
        #gene = [tx[::-1] for tx in gene]

    if reverse:
        introns = [(exon1[0], exon2[3], str(int(exon1[4]) + offset)) 
            for tx in gene for exon1, exon2 in zip(tx[:-1], tx[1:]) if len(tx) > 2]
    else:
        introns = [(exon1[0], str(int(exon1[4]) + offset), exon2[3]) 
            for tx in gene for exon1, exon2 in zip(tx[:-1], tx[1:]) if len(tx) > 2]

    introns = ['\\t'.join(intron) for intron in introns]
    if introns:
        cmd = 'bedtools sort -i <( echo -e "%s" ) | bedtools merge -i stdin ' % '\n'.join(introns)
        pp = subprocess.Popen(cmd, shell=1, executable='/bin/bash', stdout=subprocess.PIPE)
        introns = map(tuple, [ii.strip('\n').split('\t') for ii in pp.stdout])
        for intron in introns:
            yield intron

def generate_sequence(gene, gff=True, reverse=None, nested=True,
    #fasta='/home/jmerkin/genomes/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa',):
    fasta='/home/jmerkin/genomes/hg38_ucsc/hg38.fa',):
    offset, beg, end = (-1, 3, 4) if gff else (0, 1, 2)
    if not nested:
        gene = [gene,]
    if reverse is None:
        reverse = gene[0][0][6] == '-'
    for tx in gene:
        regions = [(region[0], str(int(region[beg]) + offset), region[end]) for region in tx]
        regions  = ['\t'.join(region) for region in regions]
        cmd = 'bedtools getfasta -fi %s -bed <( echo -e "%s" ) ' % (fasta, '\n'.join(regions))
        pp = subprocess.Popen(cmd, shell=1, executable='/bin/bash', stdout=subprocess.PIPE)
        fa = [ii.strip('\n') for ii in pp.stdout if not ii.startswith('>')]
        if reverse:
            fa = [seq.translate(trans)[::-1] for seq in fa]
        seq = ''.join(fa).upper()
        yield seq



