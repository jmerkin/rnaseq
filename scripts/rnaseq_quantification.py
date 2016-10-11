#!/usr/bin/env python

import os
import sys
import subprocess
from optparse import OptionParser

def get_print_cmd(read_input, fifos):
    if read_input is None:
        return 'echo no read 2'
    read_input = read_input.split(',')
    cmds = ['(', ]
    for fin in read_input:
        if fin.startswith('http'):
            #cmd = 'wget -0- -q %s' % fin if fin.endswith('gz') else 'cat %s' % fin 
            cmd = 'scp %s /dev/stdout' % fin 
            if fin.endswith('gz'):
                cmd = cmd + ' | zcat'
        else:
            cmd = 'zcat %s' % fin if fin.endswith('gz') else 'cat %s' % fin 
        cmd = cmd + ' ;'
        cmds.append(cmd)
    cmd = '%s ) ' % ' '.join(cmds)

    if len(fifos) == 1:
        cmd = '%s > %s' % (cmd, fifos[0])
    else:
        cmd = '%s | tee %s > %s' % (cmd, ' '.join(fifos[:-1]), fifos[-1])

    return cmd

def quant_main(sample_name, r1, r2, outdir, gene_file, index_dir, script_directory,
        runs={'quant': 'transcript.fa.index', 
            'full_quant': 'complete.fa.index', 
            'end_bias': 'end_bias_ratio.fa.index', 
            'intron_ratio': 'exon_intron_ratio.fa.index'}):

    # remove complete index in case it wasn't generated (i.e. no extra files given)
    runs = {rtype: index for rtype, index in runs.iteritems() 
            if os.path.exists(os.path.join(index_dir, index))}

    fifos1 = ['/tmp/%s.r1.%s.fifo' % (sample_name, num) for num in xrange(len(runs))]
    if r2 is not None:
        fifos2 = ['/tmp/%s.r2.%s.fifo' % (sample_name, num) for num in xrange(len(runs))]
    else:
        fifos2 = [None for _ in fifos1]

    full_cmd = ['mkfifo %s ; ' % fifo for fifo in (fifos1 + fifos2) if fifo is not None]

    r1_print = get_print_cmd(r1, fifos1)
    r2_print = get_print_cmd(r2, fifos2)
    full_cmd.append('%s & %s &' % (r1_print, r2_print) )

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for run_num, (run_type, index) in enumerate(runs.iteritems()):
        boots = 50 if 'quant' in run_type else 0
        fifo1 = fifos1[run_num]
        fifo2 = fifos2[run_num]
        run_output = os.path.join(outdir, run_type)
        icmd = ['salmon quant',
                '--index', os.path.join(index_dir, index),
                '--output', run_output,
                '--libType', 'A',
                '--seqBias',
                '--gcBias',
                '--numBootstraps', str(boots),
                '--dumpEq',
                '-w 200',
                '-p', '3',
                ]
        if gene_file is not None:
            icmd.extend(['--geneMap', gene_file])
        if fifo2 is not None:
            icmd.append('-1 %s -2 %s' % (fifo1, fifo2))
        else:
            icmd.append('-r %s' % fifo1)
        icmd.append(' & ')
        icmd = ' '.join(icmd)
        full_cmd.append(icmd)

    full_cmd.append(' wait ; ')
    full_cmd.append('python %s/end_bias_processing.py %s/end_bias/quant.sf ; ' % \
            (script_directory, run_output) )
    full_cmd.append('python %s/exon_intron_ratio_processing.py %s/end_bias/quant.sf ; ' % \
            (script_directory, run_output) )

    full_cmd += ['rm -f %s ; ' % fifo for fifo in (fifos1 + fifos2)]
    full_cmd = ' '.join(full_cmd)
    #print full_cmd
    subprocess.Popen(full_cmd, shell=1, 
            stdout=open('%s/log_out.txt' % outdir, 'w'),
            stderr=open('%s/log_err.txt' % outdir, 'w'),
            executable='/bin/bash',
            ).wait()

    return

if __name__ == '__main__':

    script_directory = os.path.dirname(os.path.realpath(__file__)).rstrip('scripts/')

    parser = OptionParser()
    parser.add_option('-s', '--sample', dest='sample')
    parser.add_option('-1', '--readone', dest='r1')
    parser.add_option('-2', '--readtwo', dest='r2')
    parser.add_option('-o', '--outdir', dest='output_directory')
    parser.add_option('-g', '--gene_map', dest='gene_file',
        default='%s/processed_annotations/gencode.v24.annotation.gff3.gz/transcripts.tsv' % script_directory)
    parser.add_option('-i', '--index_dir', dest='index_dir', 
        default='%s/processed_annotations/gencode.v24.annotation.gff3.gz/' % script_directory)

    op, args = parser.parse_args()

    quant_main(op.sample, op.r1, op.r2, op.output_directory, op.gene_file, op.index_dir, script_directory)

