#!/usr/bin/env python

__author__ = 'Francesco Asnicar (f.asnicar@unitn.it)'
__version__ = '0.2.13'
__date__ = '16 November 2023'


import os
import sys
import bz2
import glob
import gzip
import pandas as pd
import argparse
import subprocess as sb
import multiprocessing as mp
import time
import datetime
import shutil
import numpy as np
from Bio import SeqIO
import math
import traceback
from concurrent import futures

if sys.version_info[0] < 3:
    raise Exception("Preprocessing requires Python 3, your current Python version is {}.{}.{}"
                    .format(sys.version_info[0], sys.version_info[1], sys.version_info[2]))


def info(s, init_new_line=False, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write('{}'.format(s))
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)


def error(s, init_new_line=False, exit=False, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write('[e] {}\n'.format(s))
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)


def read_params():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('-i', '--input_dir', required=True, type=str, help="Path to input directory")
    p.add_argument('-o', '--output_dir', required=True, type=str, help="Path to output directory")
    p.add_argument('-e', '--extension', required=False, default=".fastq.gz",
                   choices=[".fastq.gz", ".fq.gz", ".fastq.bz2", ".fq.bz2"], help="The extension of the raw input files")
    p.add_argument('-m', '--mapper', required=False, type=str, help="The mapper (aligner) being used for removing human reads. The options are bowtie2, kraken or stringent (both mappers)",
                    choices=["bowtie2", "kraken", "stringent"], default="bowtie2")
    p.add_argument('-q', '--quality_control', required=False, type=str, help="The tool used for quality trimming. The options are trimgalore (default) or fastp (paired).", choices=['trimgalore', 'fastp'], default='trimgalore')
    p.add_argument('-s', '--samplename', required=False, default="", help="Specify the sample name")

    p.add_argument('-f', '--forward', required=False, default="R1",
                   help="Identifier to distinguish forward reads in the input folder")
    p.add_argument('-r', '--reverse', required=False, default="R2",
                   help="Identifier to distinguish reverse reads in the input folder")

    procs = p.add_argument_group("Params for the number of processors to use")
    procs.add_argument('-n', '--nproc', required=False, default=2, type=int, help="Number of threads to use")
    procs.add_argument('-t', '--nproc_mapper', required=False, default=2, type=int, help="Number of mapper processors")

    run = p.add_argument_group("Params for the Kraken2 run")
    run.add_argument('-c', '--run_confidence', required=False, default=0.05, type=float, help="Confidence score threshold; must be in [0, 1]")
    run.add_argument('-g', '--run_minimum_hit_groups', required=False, default=3, type=int, help="Minimum number of hit groups (overlapping k-mers")

    rm = p.add_argument_group('Available host genomes that can be removed')
    rm.add_argument('--rm_hsap', required=False, default=False, action='store_true', 
                    help="Remove H. sapiens genome")
    rm.add_argument('--rm_mmus', required=False, default=False, action='store_true', 
                    help="Remove Mus musculus C57BL/6J (black 6) genome")
    rm.add_argument('--rm_rrna', required=False, default=False, action='store_true', 
                    help="Remove rRNA (for mRNA datasets)")
    rm.add_argument('--rm_pcin', required=False, default=False, action='store_true', 
                    help="Remove Phascolarctos cinereus GCA_900166895 (Koala) genome")
    rm.add_argument('--rm_pcoq', required=False, default=False, action='store_true', 
                    help="Remove Propithecus coquereli (lemur) genome")
    rm.add_argument('--rm_mmur', required=False, default=False, action='store_true', 
                    help="Remove Microcebus murinus (grey mouse lemur) genome")
    rm.add_argument('--rm_mmul', required=False, default=False, action='store_true', 
                    help="Remove Macaca mulatta GCF_003339765 (rhesus monkey) genome")
    rm.add_argument('--rm_ptro', required=False, default=False, action='store_true', 
                    help="Remove Pan troglodytes GCF_028858775 (chimpanzee) genome")
    rm.add_argument('--rm_sbol', required=False, default=False, action='store_true', 
                    help="Remove Saimiri boliviensis GCF_016699345 (squirrel monkey) genome")
    rm.add_argument('--rm_vvar', required=False, default=False, action='store_true', 
                    help="Remove Varecia variegata GCA_028533085 (black-and-white ruffed lemur) genome")
    rm.add_argument('--rm_clup', required=False, default=False, action='store_true', 
                    help="Remove Canis lupus familiaris GCF_000002285 (dog) genome")
    rm.add_argument('--rm_sscr', required=False, default=False, action='store_true',
                    help="Remove Sus scrofa GCF_000003025 (pig) genome")
    rm.add_argument('--rm_eubi', required=False, default=False, action='store_true',
                    help="Remove both Canis lupus familiaris GCF_000002285 (dog) genome and Felis catus GCF_018350175 (cat) genome")
    rm.add_argument('--rm_mdom', required=False, default=False, action='store_true',
                    help="Remove Malus domestica GDDH13_1-1 (domestic apple) genome")
    rm.add_argument('--rm_btau', required=False, default=False, action='store_true',
                    help="Remove Bos taurus GCA_002263795.4 (domestic cattle) genome")
    rm.add_argument('--rm_fana', required=False, default=False, action='store_true',
                    help="Remove Fragaria x ananassa (cultivated octoploid strawberry) genome")
    p.add_argument('-p', '--paired_end', required=False, default=False, action='store_true',
                   help="Specify this when providing paired-end sequencing reads as input")
    p.add_argument('-k', '--keep_intermediate', required=False, default=False, action='store_true',
                   help="If specified the script won't remove intermediate files")
    p.add_argument('-x', '--bowtie2_indexes', required=False, default='/shares/CIBIO-Storage/CM/scratch/databases/bowtie2_indexes',
                   type=str, help="Folder containing the bowtie2 indexes of the genomes to be removed from the samples")
    p.add_argument('-d', '--kraken2_database', required=False, default='/shares/CIBIO-Storage/CM/scratch/databases/kraken2_databases',
                   type=str, help="Folder containing the bowtie2 indexes of the genomes to be removed from the samples")
    p.add_argument('--dry_run', required=False, default=False, action='store_true', help="Print commands do not execute them")
    p.add_argument('--verbose', required=False, default=False, action='store_true', help="Makes preprocessing verbose")
    p.add_argument('-v', '--version', action='version', version='Preprocessing version {} ({})'.format(__version__, __date__),
                   help="Prints the current Preprocessing version and exit")
    return p.parse_args()


def init_terminating(terminating_):
    global terminating
    terminating = terminating_

def parallel_execution(arguments):
    function, args = arguments
    if not terminating.is_set():
        try:
            return function(args)
        except Exception as e:
            terminating.set()
            traceback.print_exc()
            print("Error executing {}".format(arguments))
            print("Parallel execution fails: "+str(e))
    else:
        terminating.set()

def make_gen_fastq(reader):
        """Read fastq file in chunks
        
        Args:
            reader (function): the reader function
        """
        b = reader(1024 * 1024) 
        while (b):
            yield b
            b = reader(1024 * 1024)
            
def rawpycount_bases_mp(x):
    if not terminating.is_set():
        try:
            filename, samplename, output_dir, tot, dry_run, verbose, keep_headers = x
            
            if dry_run:
                return
            
            f = bz2.BZ2File(filename, 'rb') if filename.endswith(".bz2") else gzip.open(filename, 'rb') if filename.endswith('.gz') else open(filename, 'rb')
            f_gen = make_gen_fastq(f.read)

            intype = 'fasta' if 'fasta' in filename else 'fastq'

            # Determine which lines to count based on the input type
            i = 4 if intype == 'fastq' else 2
            j = 2 if intype == 'fastq' else 0

            nbases = 0 if tot else list()
            read_headers = list()
            read_num = -1
            line_num = 0
            buffer = b''

            for chunk in f_gen:
                lines = (buffer + chunk).split(b'\n')
                buffer = lines[-1]  # Save the last incomplete line to be processed in the next chunk   
                for line in lines[:-1]:  # Process all but the last line
                    line_num += 1
                    if (line_num+1) % i == j and keep_headers:
                        read_headers.append(line.decode('ascii').split(" ")[0][1:])
                    if line_num % i == j:
                        if tot:
                            nbases += len(line)
                        else:
                            read_num += 1
                            nbases.append(len(line))

            # Check if there's remaining content in the buffer
            if buffer:
                line_num += 1
                if (line_num+1) % i == j and keep_headers:
                        read_headers.append(line.decode('ascii').split(" ")[0][1:])
                if line_num % i == j:
                    if tot:
                        nbases += len(line)
                    else:
                        read_num += 1
                        nbases.append(len(line))
            f.close()
            
            samplename = "/".join(filename.split("/")[:-1]) + "/" + samplename + "." + filename.split("/")[-1][len(samplename)+1:]
            samplename = ".".join(samplename.split(".")[:-2]) if ".gz" in samplename or ".bz2" in samplename else ".".join(samplename.split(".")[:-1])
            samplename = "{}/{}".format(output_dir, samplename.split("/")[-1])
            print(samplename + ".stats")
            
            with open(samplename + ".stats", 'w') as outf:
                outf.write('\t'.join(["#samplename", "n_of_bases", "n_of_reads", "min_read_len", "median_read_len", "mean_read_len", "max_read_len"]) + '\n')
                outf.write('\t'.join([str(a) for a in [samplename.split("/")[-1], np.sum(nbases), len(nbases), np.min(nbases), np.median(nbases), np.mean(nbases), np.max(nbases)]]) + '\n')
        
            if keep_headers: return read_headers
        
        except Exception as e:
            terminating.set()
            
            error('rawpycount_bases_mp()\n    x: {}\n    e: {}'.format(x, str(e)), init_new_line=True)
            raise
    else:
        terminating.set()

def pool_rawpycount_bases(arg):
    terminating = mp.Event()

    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, arg, chunksize=2)]
        except Exception as e:
            error('rawpycount_bases()\ntasks: {}\n    e: {}'.format(arg, e), init_new_line=True, exit=True)        

def rawpycount_bases(inputs_r1s_r2s, samplename, output_dir, tot, dry_run=False, verbose=False, keep_headers=False):
    if dry_run or verbose:
        info('rawpycount_bases()\n', init_new_line=True)
    
    tasks = [(R[0], samplename, output_dir, tot, dry_run, verbose, keep_headers) for R in inputs_r1s_r2s]

    return pool_rawpycount_bases(((rawpycount_bases_mp, args) for args in tasks))            
            
def check_params(args):
    if not os.path.isdir(args.input_dir):
        error('input folder "{}" does not exists'.format(args.input_dir), exit=True)

    if args.input_dir.endswith('/'):
        args.input_dir = args.input_dir[:-1]


def initt(terminating_):
    # This places terminating in the global namespace of the worker subprocesses.
    # This allows the worker function to access `terminating` even though it is
    # not passed as an argument to the function.
    global terminating
    terminating = terminating_


def preflight_check(dry_run=False, verbose=False):
    if dry_run or verbose:
        info('preflight_check()\n', init_new_line=True)

    cmds = ['fna_len.py -h', 'trim_galore --help', 'bowtie2 -h',
            'kraken2 -h', 'split_and_sort.new.update.test2.py -h',
            'cat_stats.py -h']

    for cmd in cmds:
        if dry_run or verbose:
            info('{}\n'.format(cmd))

        if dry_run:
            continue

        try:
            with open(os.devnull, 'w') as devnull:
                sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
        except Exception as e:
            error('preflight_check()\n{}\n{}'.format(cmd, e), exit=True)


def get_inputs(input_dir, fwd, rev, sn, ext, verbose=False):
    if verbose:
        info('get_inputs()\n', init_new_line=True)

    R1, R2 = [], []

    if (fwd not in sn) and (rev not in sn):
        R1 = [os.path.join(input_dir, i) for i in os.listdir(input_dir) if (fwd in i) and i.endswith(ext)]
        R2 = [os.path.join(input_dir, i) for i in os.listdir(input_dir) if (rev in i) and i.endswith(ext)]
    else:
        count_fwd = sn.count(fwd) + 1 if fwd in sn else 1
        count_rev = sn.count(rev) + 1 if rev in sn else 1

        for i in os.listdir(input_dir):
            if not i.endswith(ext):
                continue

            if i.count(fwd) == count_fwd:
                R1.append(os.path.join(input_dir, i))
            elif i.count(rev) == count_rev:
                R2.append(os.path.join(input_dir, i))

    return (sorted(R1), sorted(R2))

def write_gz(final_file, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('write_gz()\n', init_new_line=True)
    
    tasks = [(final_file, "1", dry_run, verbose),
             (final_file, "2", dry_run, verbose)]
    return pool_write_gz(((write_gz_mp, args) for args in tasks))

def pool_write_gz(arg):
    terminating = mp.Event()

    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, arg, chunksize=2)]
        except Exception as e:
            error('write_gz()\ntasks: {}\n    e: {}'.format(arg, e), init_new_line=True, exit=True)
    
def write_gz_mp(x):
    if not terminating.is_set():
        try:
            final_file, R, dry_run, verbose = x
            start = datetime.datetime.now()
            
            cmd = ('pigz -p {} {}').format(args.nproc, final_file.replace("#",R))
            
            if dry_run or verbose:
                info('{}\n'.format(cmd))

            if not dry_run:
                with open(os.devnull, 'w') as devnull:
                    sb.check_call(cmd.split(' '), stdout=devnull)
                    os.rename("{}.gz".format(final_file.replace("#",R)), "{}.fastq.gz".format(".".join(final_file.replace("#",R).split(".")[:-1])))
#                 with gzip.open("{}.fastq.gz".format(".".join(final_file.replace("#",R).split(".")[:-1])), 'wt') as outf:
#                     for record in SeqIO.parse(final_file.replace("#",R), format='fastq'):
#                         SeqIO.write(record, outf, 'fastq')
#             if not args.keep_intermediate:
#                 os.remove(final_file.replace("#",R))
#             else:
#                 print("{}.fastq.gz".format(".".join(final_file.replace("#",R).split(".")[:-1])))
            print("{}.fastq.gz: ".format(".".join(final_file.replace("#",R).split(".")[:-1])), datetime.datetime.now() - start)
        except Exception as e:
            terminating.set()

            error('write_gz_mp()\n    x: {}\n    e: {}'.format(x, str(e)), init_new_line=True)
            raise
    else:
        terminating.set()

def concatenate_reads(input_dir, inputs_r1s_r2s, nproc=1, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('concatenate_reads()\n', init_new_line=True)

    out_prefix = os.path.basename(input_dir)
    R1s, R2s = inputs_r1s_r2s
    outR1fastq = '{}.R1.fastq'.format(os.path.join(input_dir, out_prefix))
    outR2fastq = '{}.R2.fastq'.format(os.path.join(input_dir, out_prefix))
    outR1stats = '{}.R1.stats'.format(os.path.join(input_dir, out_prefix))
    outR2stats = '{}.R2.stats'.format(os.path.join(input_dir, out_prefix))
    tasks = [(R1s, outR1fastq, outR1stats, dry_run, verbose),
             (R2s, outR2fastq, outR2stats, dry_run, verbose)]

    return pool_concatenate_reads(((concatenate_reads_mp, args) for args in tasks))

def pool_concatenate_reads(arg):
    terminating = mp.Event()

    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, arg, chunksize=2)]
        except Exception as e:
            error('concatenate_reads()\ntasks: {}\n    e: {}'.format(arg, e), init_new_line=True, exit=True)

def concatenate_reads_mp(x):
    if not terminating.is_set():
        try:
            inps, out_fastq, out_stats, dry_run, verbose = x

            # decompress
            if not os.path.isfile(out_fastq):
                if dry_run or verbose:
                    info('cat {} > {}\n'.format(' '.join(inps), out_fastq))

                if not dry_run:
                    g = open(out_fastq, 'wb')

                    for inpR in inps:
                        # decompress input file
                        if inpR.endswith('.bz2'):
                            with bz2.open(inpR, 'rb') as f:
                                shutil.copyfileobj(f, g)
                        elif inpR.endswith('.gz'):
                            with gzip.open(inpR, 'rb') as f:
                                shutil.copyfileobj(f, g)

                    g.close()

            # stats
            if not os.path.isfile(out_stats):
                cmd = 'fna_len.py -q --stat {} {}'.format(out_fastq, out_stats)

                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    sb.check_call(cmd.split(' '))

            return os.path.basename(out_fastq)
        except Exception as e:
            terminating.set()

            for i in [out_fastq, out_stats]:
                if os.path.exists(i):
                    os.remove(i)

            error('concatenate_reads_mp()\n    x: {}\n    e: {}'.format(x, str(e)), init_new_line=True)
            raise
    else:
        terminating.set()


def quality_control(input_dir, merged_r1_r2, output_dir, keep_intermediate, sn, nproc=1, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('quality_control()\n', init_new_line=True)

    if not args.paired_end:
        tasks = zip(merged_r1_r2, [input_dir] * len(merged_r1_r2),
                    [output_dir] * len(merged_r1_r2),
                [keep_intermediate] * len(merged_r1_r2),
                [dry_run] * len(merged_r1_r2),
                [verbose] * len(merged_r1_r2))
        r1, r2 = pool_quality_control(((quality_control_mp, args) for args in tasks))

        qc = [r1[0],r2[0]]
        reads_headers = [r1[1],r2[1]]

        samplename = "_".join(".".join(merged_r1_r2[0][0].split(".")[:-2]).split("_")[:-1]).split("/")[-1]
        unpaired_file = samplename + '_unpaired.txt.bz2'

        with bz2.open(os.path.join(output_dir, unpaired_file), 'wt') as f:
            sym_difference = set(reads_headers[0]).symmetric_difference(set(reads_headers[1]))
            f.write("\n".join( (i for i in sym_difference) ) + "\n")

        fwd, rev = 'R1', 'R2'
        r1, r2 = [], []
        count_fwd = sn.count(fwd) + 1 if fwd in sn else 1
        count_rev = sn.count(rev) + 1 if rev in sn else 1

        for i in qc:
            if i.count(fwd) == count_fwd:
                r1.append(i)
            elif i.count(rev) == count_rev:
                r2.append(i)

        if len(r1) > 1:
            error('quality_control(): more than one R1 detected: [{}]'.format(', '.join(r1)), exit=True)

        if len(r2) > 1:
            error('quality_control(): more than one R2 detected: [{}]'.format(', '.join(r2)), exit=True)

        return tuple([r1[0], r2[0]]), unpaired_file
    else:
        if args.quality_control == "trimgalore":
            tasks = (merged_r1_r2, input_dir, output_dir,
                    keep_intermediate,
                    dry_run,
                    verbose)
            r1, r2 = quality_control_paired(tasks)

            qc = [r1,r2]

            fwd, rev = '_R1', '_R2'
            r1, r2 = [], []
            count_fwd = sn.count(fwd) + 1 if fwd in sn else 1
            count_rev = sn.count(rev) + 1 if rev in sn else 1

            for i in qc:
                if i.split("/")[-1].count(fwd) == count_fwd:
                    r1.append(i)
                elif i.split("/")[-1].count(rev) == count_rev:
                    r2.append(i)

            if len(r1) > 1:
                error('quality_control(): more than one R1 detected: [{}]'.format(', '.join(r1)), exit=True)

            if len(r2) > 1:
                error('quality_control(): more than one R2 detected: [{}]'.format(', '.join(r2)), exit=True)

            return tuple([r1[0], r2[0]]), None
        if args.quality_control == "fastp":
            tasks = (merged_r1_r2, input_dir, output_dir, keep_intermediate, dry_run, verbose)
            r1, r2 = quality_control_fastp(tasks)

            qc = [r1,r2]

            fwd, rev = '_R1', '_R2'
            r1, r2 = [], []
            count_fwd = sn.count(fwd) + 1 if fwd in sn else 1
            count_rev = sn.count(rev) + 1 if rev in sn else 1

            for i in qc:
                if i.split("/")[-1].count(fwd) == count_fwd:
                    r1.append(i)
                elif i.split("/")[-1].count(rev) == count_rev:
                    r2.append(i)

            if len(r1) > 1:
                error('quality_control(): more than one R1 detected: [{}]'.format(', '.join(r1)), exit=True)

            if len(r2) > 1:
                error('quality_control(): more than one R2 detected: [{}]'.format(', '.join(r2)), exit=True)

            return tuple([r1[0], r2[0]]), None

def pool_quality_control(arg):
    terminating = mp.Event()

    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, arg, chunksize=2)]
        except Exception as e:
            error('quality_control()\ntasks: {}\n    e: {}'.format(arg, e), init_new_line=True, exit=True)

def quality_control_fastp(x):
    R, input_dir, output_dir, keep_intermediate, dry_run, verbose = x
    R1 = R[0][0]
    R2 = R[1][0]
            
    oR1 = ".".join(R1.split(".")[:-2]).split("/")[-1]
    oR2 = ".".join(R2.split(".")[:-2]).split("/")[-1]
    
    if not (os.path.isfile('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR1))) and os.path.isfile('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR2)))):
        cmd = ('fastp -i {} -o {}_trimmed.fq.gz -I {} -O {}_trimmed.fq.gz -z 4 -e 20 -l 75 -w {} -q 20 --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --trim_poly_g --trim_poly_x -y -n 2 -h {}.html -j {}.json').format(os.path.join(input_dir, R1), os.path.join(output_dir, oR1), os.path.join(input_dir, R2), os.path.join(output_dir, oR2), args.nproc, os.path.join(output_dir, oR1[:-3]), os.path.join(output_dir, oR1[:-3]))
        
        if dry_run or verbose:
                info('{}\n'.format(cmd))

        if not dry_run:
            with open(os.devnull, 'w') as devnull:
                sb.check_call(cmd.split(' '), stdout=devnull)

    rawpycount_bases([["{}_trimmed.fq.gz".format(os.path.join(output_dir, oR1))], ["{}_trimmed.fq.gz".format(os.path.join(output_dir, oR2))]], args.samplename, output_dir, False, dry_run, verbose, False)

    return ('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR1)), '{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR2)))
            
def quality_control_mp(x):
    trimmed_reads = []
    if not terminating.is_set():
        try:
            R, input_dir, output_dir, keep_intermediate, dry_run, verbose = x
            R = R[0]
            
            oR = ".".join(R.split(".")[:-2]).split("/")[-1]

            if not os.path.isfile('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR))):
                # # METATRANSCRIPTOMES (RNA) as they have usually shorter fragments - cut length 50
                # cmd = ('trim_galore --nextera --stringency 5 --length 50 --2colour 20 --max_n 2 --trim-n -j 1 --dont_gzip '
                #        '--no_report_file --suppress_warn --output_dir {} {}').format(input_dir, os.path.join(input_dir, R))

                # STANDARD DNA (meta)genomes
                cmd = ('trim_galore --nextera --stringency 5 --length 75 --2colour 20 --max_n 2 --trim-n -j 1 '
                       '--no_report_file --suppress_warn --output_dir {} {}').format(output_dir, os.path.join(input_dir, R))

                # # command for Moreno, no --nextera
                # cmd = ('trim_galore --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip '
                #        '--no_report_file --suppress_warn --output_dir {} {}').format(input_dir, os.path.join(input_dir, R))

                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    with open(os.devnull, 'w') as devnull:
                        sb.check_call(cmd.split(' '), stdout=devnull)

            reads_headers = rawpycount_bases_mp(["{}_trimmed.fq.gz".format(os.path.join(output_dir, oR)), args.samplename, output_dir, False, dry_run, verbose, True])
            
            return ('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR)), reads_headers)

        except Exception as e:
            terminating.set()

            error('quality_control_mp()\n    x: {}\n    e: {}'.format(x, e), init_new_line=True)
            raise

    else:
        terminating.set()

def quality_control_paired(x):
    R, input_dir, output_dir, keep_intermediate, dry_run, verbose = x
    R1 = R[0][0]
    R2 = R[1][0]

    oR1 = ".".join(R1.split(".")[:-2]).split("/")[-1]
    oR2 = ".".join(R2.split(".")[:-2]).split("/")[-1]

    if not (os.path.isfile('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR1))) or os.path.isfile('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR2)))):
        # # METATRANSCRIPTOMES (RNA) as they have usually shorter fragments - cut length 50
        # cmd = ('trim_galore --nextera --stringency 5 --length 50 --2colour 20 --max_n 2 --trim-n -j 1 --dont_gzip '
        #        '--no_report_file --suppress_warn --output_dir {} {}').format(input_dir, os.path.join(input_dir, R))

        # STANDARD DNA (meta)genomes
        cmd = ('trim_galore --nextera --stringency 5 --length 75 --2colour 20 --max_n 2 --trim-n -j 1 --no_report_file --suppress_warn --paired --output_dir {} {} {}').format(output_dir, os.path.join(input_dir, R1), os.path.join(input_dir, R2))

        # # command for Moreno, no --nextera
        # cmd = ('trim_galore --stringency 5 --length 75 --quality 20 --max_n 2 --trim-n --dont_gzip '
        #        '--no_report_file --suppress_warn --output_dir {} {}').format(input_dir, os.path.join(input_dir, R))

        if dry_run or verbose:
            info('{}\n'.format(cmd))

        if not dry_run:
            with open(os.devnull, 'w') as devnull:
                sb.check_call(cmd.split(' '), stdout=devnull)

        if not dry_run:
            os.rename("{}_val_1.fq.gz".format(os.path.join(output_dir, oR1)), "{}_trimmed.fq.gz".format(os.path.join(output_dir, oR1)))
            os.rename("{}_val_2.fq.gz".format(os.path.join(output_dir, oR2)), "{}_trimmed.fq.gz".format(os.path.join(output_dir, oR2)))
    
    rawpycount_bases([["{}_trimmed.fq.gz".format(os.path.join(output_dir, oR1))], ["{}_trimmed.fq.gz".format(os.path.join(output_dir, oR2))]], args.samplename, output_dir, False, dry_run, verbose, False)

    return ('{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR1)), '{}_trimmed.fq.gz'.format(os.path.join(output_dir, oR2)))

def screen_contaminating_dnas(input_dir, qced_r1_r2, bowtie2_indexes, keep_intermediate, rm_hsap, rm_rrna, rm_mmus, rm_pcin, rm_pcoq,
                              rm_mmur, rm_mmul, rm_ptro, rm_sbol, rm_vvar, rm_clup, rm_sscr, rm_eubi, rm_btau, rm_mdom, rm_fana, kraken2_database, confidence_score, minimum_hit_groups,
                              mapper, nprocs_mapper=1, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('screen_contaminating_dnas()\n', init_new_line=True)
    if not args.paired_end:
        tasks = zip(qced_r1_r2, [input_dir] * len(qced_r1_r2),
                    [bowtie2_indexes] * len(qced_r1_r2),
                    [keep_intermediate] * len(qced_r1_r2),
                    [rm_hsap] * len(qced_r1_r2),
                    [rm_rrna] * len(qced_r1_r2),
                    [rm_mmus] * len(qced_r1_r2),
                    [rm_pcin] * len(qced_r1_r2),
                    [rm_pcoq] * len(qced_r1_r2),
                    [rm_mmur] * len(qced_r1_r2),
                    [rm_mmul] * len(qced_r1_r2),
                    [rm_ptro] * len(qced_r1_r2),
                    [rm_sbol] * len(qced_r1_r2),
                    [rm_vvar] * len(qced_r1_r2),
                    [rm_clup] * len(qced_r1_r2),
                    [rm_sscr] * len(qced_r1_r2),
                    [rm_eubi] * len(qced_r1_r2),
                    [rm_btau] * len(qced_r1_r2),
                    [rm_mdom] * len(qced_r1_r2),
                    [rm_fana] * len(qced_r1_r2),
                    [kraken2_database] * len(qced_r1_r2),
                    [confidence_score] * len(qced_r1_r2),
                    [minimum_hit_groups] * len(qced_r1_r2),
                    [mapper] * len(qced_r1_r2),
                    [nprocs_mapper] * len(qced_r1_r2),
                    [dry_run] * len(qced_r1_r2),
                    [verbose] * len(qced_r1_r2))

        return pool_screen_contaminating_dnas(((screen_contaminating_dnas_mp, args) for args in tasks))
    else:
        tasks = (qced_r1_r2, input_dir,
                    bowtie2_indexes,
                    keep_intermediate,
                    rm_hsap,
                    rm_rrna,
                    rm_mmus,
                    rm_pcin,
                    rm_pcoq,
                    rm_mmur,
                    rm_mmul,
                    rm_ptro,
                    rm_sbol,
                    rm_vvar,
                    rm_clup,
                    rm_sscr,
                    rm_eubi,
                    rm_btau,
                    rm_fana,
                    rm_mdom,
                    kraken2_database,
                    confidence_score,
                    minimum_hit_groups,
                    mapper,
                    nprocs_mapper,
                    dry_run,
                    verbose)
        return screen_contaminating_dnas_paired(tasks)

def pool_screen_contaminating_dnas(arg):
    terminating = mp.Event()

    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, arg, chunksize=2)]
        except Exception as e:
            error('screen_contaminating_dnas()\ntasks: {}\n    e: {}'.format(arg, e), init_new_line=True, exit=True)

def screen_contaminant_bowtie2(x):
    R, input_dir, bowtie2_indexes, keep_intermediate, rm_hsap, rm_rrna, rm_mmus, rm_pcin, rm_pcoq, rm_mmur, rm_mmul, rm_ptro, rm_sbol, rm_vvar, rm_clup, rm_sscr, rm_eubi, rm_btau, rm_fana, rm_mdom, kraken2_database, confidence_score, minimum_hit_groups, mapper, nprocs_mapper, dry_run, verbose = x
#   screened = []
    to_removes = []
#   cont_dnas = ['phiX174']
    cont_dnas = []

    if rm_hsap:
        cont_dnas += ['hg38']

    if rm_mmus:
        cont_dnas += ['mmusculus_black6_GCA_000001635_8']

    if rm_rrna:
        cont_dnas += ['SILVA_138.1_LSURef_NR99_tax_silva_DNA', 'SILVA_138.1_SSURef_NR99_tax_silva_DNA']

    if rm_pcin:
        cont_dnas += ['Phascolarctos_cinereus__GCA_900166895.1__tgac_v2.0']

    if rm_pcoq:
        cont_dnas += ['Propithecus_coquereli_GCF_000956105.1']

    if rm_mmur:
        cont_dnas += ['Microcebus_murinus_GCF_000165445.2']

    if rm_mmul:
        cont_dnas += ['mmulatta_GCF_003339765.1']

    if rm_ptro:
        cont_dnas += ['ptroglodytes_GCF_028858775.1']

    if rm_sbol:
        cont_dnas += ['sboliviensis_GCF_016699345.2']

    if rm_vvar:
        cont_dnas += ['vvariegata_GCA_028533085.1']

    if rm_clup:
        cont_dnas += ['CanFam3.1']

    if rm_sscr:
        cont_dnas += ['Sscrofa11.1']
    
    if rm_eubi:
        cont_dnas += ['Eubiome']

    if rm_btau:
        cont_dnas += ['bosTau9-ARS-UCD2.0']
    
    if rm_mdom:
        cont_dnas += ['Mdomestica_GDDH13_1-1']
    
    if rm_fana:
        cont_dnas += ['F_ana_Camarosa_6_28_17']

    outf = R[:R.rfind('.')]
    Rext = ".".join(R.split(".")[-2:])
    final = None

    for cont_dna in cont_dnas:
        iR = outf
        
        suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
        outf += suffix

        if not os.path.isfile('{}_bt2.fastq.gz'.format(os.path.join(input_dir, outf))):
                    # -x is the argument for indices
                    # --un <path>        write unpaired reads that didn't align to <path>
                    # -p is the cores number
                    # -U input fastq files
            cmd = ('bowtie2 -x {} -U {} -p {} --sensitive-local --un-gz {}_bt2.fastq.gz'
                    .format(os.path.join(bowtie2_indexes, cont_dna), os.path.join(input_dir, iR + Rext), nprocs_mapper,
                    os.path.join(input_dir, outf)))


            if dry_run or verbose:
                info('{}\n'.format(cmd))

            if not dry_run:
                try:
                    with open(os.devnull, 'w') as devnull:
                        sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
                except Exception as e:
                    if os.path.exists('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf))):
                        os.remove('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf)))

                    error('screen_contaminating_dnas()\n{}\n{}'.format(cmd, e), exit=True)

        rawpycount_bases_mp(["{}_bt2.fastq.gz".format(os.path.join(input_dir, outf)), args.samplename, output_dir, False, dry_run, verbose, False])
                    
        if not keep_intermediate:
            to_removes.append(iR + Rext)

        Rext = '_bt2.fastq.gz'
        final = outf + '_bt2.fastq.gz'

    remove(to_removes, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)

    return final

def screen_contaminant_bowtie2_paired(x):
    R, input_dir, bowtie2_indexes, keep_intermediate, rm_hsap, rm_rrna, rm_mmus, rm_pcin, rm_pcoq, rm_mmur, rm_mmul, rm_ptro, rm_sbol, rm_vvar, rm_clup, rm_sscr, rm_eubi, rm_btau, rm_fana, rm_mdom, kraken2_database, confidence_score, minimum_hit_groups, mapper, nprocs_mapper, dry_run, verbose = x
#   screened = []
    to_removes = []
#   cont_dnas = ['phiX174']
    cont_dnas = []

    if rm_hsap:
        cont_dnas += ['hg38']

    if rm_mmus:
        cont_dnas += ['mmusculus_black6_GCA_000001635_8']

    if rm_rrna:
        cont_dnas += ['SILVA_138.1_LSURef_NR99_tax_silva_DNA', 'SILVA_138.1_SSURef_NR99_tax_silva_DNA']

    if rm_pcin:
        cont_dnas += ['Phascolarctos_cinereus__GCA_900166895.1__tgac_v2.0']

    if rm_pcoq:
        cont_dnas += ['Propithecus_coquereli_GCF_000956105.1']

    if rm_mmur:
        cont_dnas += ['Microcebus_murinus_GCF_000165445.2']

    if rm_mmul:
        cont_dnas += ['mmulatta_GCF_003339765.1']

    if rm_ptro:
        cont_dnas += ['ptroglodytes_GCF_028858775.1']

    if rm_sbol:
        cont_dnas += ['sboliviensis_GCF_016699345.2']

    if rm_vvar:
        cont_dnas += ['vvariegata_GCA_028533085.1']

    if rm_clup:
        cont_dnas += ['CanFam3.1']

    if rm_sscr:
        cont_dnas += ['Sscrofa11.1']
    
    if rm_eubi:
        cont_dnas += ['Eubiome']

    if rm_btau:
        cont_dnas += ['bosTau9-ARS-UCD2.0']
    
    if rm_mdom:
        cont_dnas += ['Mdomestica_GDDH13_1-1']
        
    if rm_fana:
        cont_dnas += ['F_ana_Camarosa_6_28_17']

    outf1 = ".".join(R[0].split(".")[:-2]) if 'gz' in R[0] or 'bz2' in R[0] else ".".join(R[0].split(".")[:-2])
    outf2 = ".".join(R[1].split(".")[:-2]) if 'gz' in R[1] or 'bz2' in R[1] else ".".join(R[1].split(".")[:-2])
    Rext = "." + ".".join(R[0].split(".")[-2:]) if 'gz' in R[0] or 'bz2' in R[0] else "." + ".".join(R[0].split(".")[-1:])
    
    final = None

    for cont_dna in cont_dnas:
        iR1 = outf1
        iR2 = outf2
        
        suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
        outf1 += suffix
        outf2 += suffix

        if not (os.path.isfile('{}_bt2.fastq.bz2'.format(os.path.join(input_dir, outf1))) or os.path.isfile('{}_bt2.fastq.bz2'.format(os.path.join(input_dir, outf2)))):
                    # -x is the argument for indices
                    # --un <path>        write unpaired reads that didn't align to <path>
                    # -p is the cores number
                    # -U input fastq files
            cmd = ('bowtie2 -x {} -1 {} -2 {} -p {} --sensitive-local | SplitUnmappedSAM.py -1 {}_bt2.fastq.gz -2 {}_bt2.fastq.gz'
                    .format(os.path.join(bowtie2_indexes, cont_dna), os.path.join(input_dir, iR1 + Rext), os.path.join(input_dir, iR2 + Rext), nprocs_mapper, os.path.join(input_dir, outf1), os.path.join(input_dir, outf2)))


            if dry_run or verbose:
                info('{}\n'.format(cmd))

            if not dry_run:
                try:
                    with open(os.devnull, 'w') as devnull:
                         sb.check_output(cmd, shell=True, stderr=devnull)
#                         sb.run(cmd.split('|')[1].split(" "), stdin=sb.run(cmd.split('|')[0].split(" "), stdout=sb.PIPE, stderr=devnull))
#                         sb.Popen(cmd, shell=True,
#                          stdin=sb.PIPE,
#                          stdout=sb.PIPE,
#                          stderr=sb.PIPE)
                except Exception as e:
                    if os.path.exists('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf1))):
                        os.remove('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf1)))
                    if os.path.exists('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf2))):
                        os.remove('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf2)))

                    error('screen_contaminating_dnas()\n{}\n{}'.format(cmd, e), exit=True)

        rawpycount_bases([["{}_bt2.fastq.gz".format(os.path.join(input_dir, outf1))], ["{}_bt2.fastq.gz".format(os.path.join(input_dir, outf2))]], args.samplename, output_dir, False, dry_run, verbose, False)
                    
        if not keep_intermediate:
            to_removes.append(iR1 + Rext)
            to_removes.append(iR2 + Rext)

        Rext = "_bt2.fastq.gz"
#     final = outf + '_bt2.fastq.bz2'
    
    if not (os.path.isfile("{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R1"))) or os.path.isfile("{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R2")))):
        os.rename("{}_bt2.fastq.gz".format(os.path.join(input_dir, outf1)), "{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R1")))
        os.rename("{}_bt2.fastq.gz".format(os.path.join(input_dir, outf2)), "{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R2")))
        
    if not (os.path.isfile("{}.stats".format(os.path.join(input_dir, args.samplename + "_R1"))) or os.path.isfile("{}.stats".format(os.path.join(input_dir, args.samplename + "_R2")))):
        stats_R1 = pd.read_csv("{}_bt2.stats".format(os.path.join(input_dir, "_".join(outf1.split("_")[:-(len(cont_dnas)+2)]) + ".R1_" + "_".join(outf1.split("_")[-(len(cont_dnas)+1):]))), sep='\t')
        stats_R1['#samplename'] = [".".join(i.split(".")[:-1]) + "_" + i.split(".")[-1][:2] + ".fastq" for i in stats_R1['#samplename']]
        stats_R1.to_csv("{}.stats".format(os.path.join(input_dir, args.samplename + "_R1")), sep='\t', index=None)
        stats_R2 = pd.read_csv("{}_bt2.stats".format(os.path.join(input_dir, "_".join(outf2.split("_")[:-(len(cont_dnas)+2)]) + ".R2_" + "_".join(outf2.split("_")[-(len(cont_dnas)+1):]))), sep='\t')
        stats_R2['#samplename'] = [".".join(i.split(".")[:-1]) + "_" + i.split(".")[-1][:2] + ".fastq" for i in stats_R2['#samplename']]
        stats_R2.to_csv("{}.stats".format(os.path.join(input_dir, args.samplename + "_R2")), sep='\t', index=None)

    remove(to_removes, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)

    return (outf1 + Rext, outf2 + Rext)

def screen_contaminant_kraken(x):
    R, input_dir, bowtie2_indexes, keep_intermediate, rm_hsap, rm_rrna, rm_mmus, rm_pcin, rm_pcoq, rm_mmur, rm_mmul, rm_ptro, rm_sbol, rm_vvar, rm_clup, rm_sscr, rm_eubi, rm_btau, rm_mdom, rm_fana, kraken2_database, confidence_score, minimum_hit_groups, mapper, nprocs_mapper, dry_run, verbose = x
    
    screened = []
    to_removes = []
#   cont_dnas = ['phiX174']
    cont_dnas = []
#    intermediate_files = []
#    kraken_reads, trimmed_reads = [],[]

    if rm_hsap:
        cont_dnas += ['hg38']
    
    if rm_rrna:
        cont_dnas += ['SILVA_138.1_LSURef_NR99_tax_silva_DNA', 'SILVA_138.1_SSURef_NR99_tax_silva_DNA']

    if rm_mmus:
        cont_dnas += ['mmusculus_black6_GCA_000001635_8']
    
    if rm_clup:
        cont_dnas += ['CanFam3.1']
    
    if rm_eubi:
        cont_dnas += ['Eubiome']

    if rm_btau:
        cont_dnas += ['bosTau9-ARS-UCD2.0']
    
    if rm_mdom:
        cont_dnas += ['Mdomestica_GDDH13_1-1']
    
    if rm_fana:
        cont_dnas += ['F_ana_Camarosa_6_28_17']

    outf = ".".join(R.split(".")[:-2])
    
    Rext = ".".join(R.split(".")[-2:])
    final = None

    for cont_dna in cont_dnas:
        iR = outf
        suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
        outf += suffix

        if not os.path.isfile('{}.fastq'.format(os.path.join(input_dir, outf))):
            cmd = ('kraken2  --db {}  --confidence {}  --minimum-hit-groups {}  --threads  {} --gzip-compressed --unclassified-out {}_krk2.fastq -- {}'
                    .format(os.path.join(kraken2_database, cont_dna), confidence_score, minimum_hit_groups, nprocs_mapper,
                    os.path.join(input_dir, outf), os.path.join(input_dir, iR + "." + Rext)))

            if dry_run or verbose:
                info('{}\n'.format(cmd))

            if not dry_run:
                try:
                    with open(os.devnull, 'w') as devnull:
                        sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
                except Exception as e:
                    if os.path.exists('{0}_krk2.fastq'.format(os.path.join(input_dir, outf))):
                        os.remove('{0}_krk2.fastq'.format(os.path.join(input_dir, outf)))

                    error('screen_contaminating_dnas()\n{}\n{}'.format(cmd, e), exit=True)

        rawpycount_bases_mp(["{}_krk2.fastq".format(os.path.join(input_dir, outf)), args.samplename, output_dir, False, dry_run, verbose, False])

        if not keep_intermediate:
            to_removes.append(iR + "." + Rext)
            
        Rext = '_krk2.fastq'
#        intermediate_files.append(outf + '_krk2.fastq')
        final = outf + '_krk2.fastq'
    
    file = open('{0}.stats'.format(os.path.join(input_dir, args.samplename + "." + iR.split("/")[-1][len(args.samplename)+1:])), "r")
    file.readline()
    trimmed_reads = float(file.readline().split("\t")[2])
    file.close()

    file = open('{0}_krk2.stats'.format(os.path.join(input_dir, args.samplename + "." + outf.split("/")[-1][len(args.samplename)+1:])), "r")
    file.readline()
    kraken_reads = float(file.readline().split("\t")[2])
    file.close()

    host = (trimmed_reads - kraken_reads) / trimmed_reads
    print("\nPercentage of host reads: {:.2f}%".format(host*100))

    if host < 0.75 and args.mapper!="stringet":
        remove(to_removes, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)
        return final
    else:
        cont_dnas = []
        if rm_hsap:
            cont_dnas += ['hg38']
            
        if rm_rrna:
            cont_dnas += ['SILVA_138.1_LSURef_NR99_tax_silva_DNA', 'SILVA_138.1_SSURef_NR99_tax_silva_DNA']
        
        if rm_mmus:
            cont_dnas += ['mmusculus_black6_GCA_000001635_8']
        
        if rm_clup:
            cont_dnas += ['CanFam3.1']

        if rm_eubi:
            cont_dnas += ['Eubiome']

        if rm_btau:
            cont_dnas += ['bosTau9-ARS-UCD2.0']
            
        if rm_mdom:
            cont_dnas += ['Mdomestica_GDDH13_1-1']
        
        if rm_fana:
            cont_dnas += ['F_ana_Camarosa_6_28_17']

        outf = final[:final.rfind('.')]
        Rext = final[final.rfind('.'):]
        final = None

        for cont_dna in cont_dnas:
            iR = outf
            suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
            outf += suffix

            if not os.path.isfile('{}_bt2.fastq.gz'.format(os.path.join(input_dir, outf))):
                        # -x is the argument for indices
                        # --un <path>        write unpaired reads that didn't align to <path>
                        # -p is the cores number
                        # -U input fastq files
                cmd = ('bowtie2 -x {} -U {} -p {} --sensitive-local --un-gz {}_bt2.fastq.gz'
                        .format(os.path.join(bowtie2_indexes, cont_dna), os.path.join(input_dir, iR + Rext), nprocs_mapper,
                        os.path.join(input_dir, outf)))


                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    try:
                        with open(os.devnull, 'w') as devnull:
                            sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
                    except Exception as e:
                        if os.path.exists('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf))):
                            os.remove('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf)))

                        error('screen_contaminating_dnas()\n{}\n{}'.format(cmd, e), exit=True)

            rawpycount_bases_mp(["{}_bt2.fastq.gz".format(os.path.join(input_dir, outf)), args.samplename, output_dir, False, dry_run, verbose, False])

            if not keep_intermediate:
                to_removes.append(iR + Rext)
                
            Rext = '_bt2.fastq.gz'
            final = outf + '_bt2.fastq.gz'
        
        print(to_removes)

        remove(to_removes, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)

        return final

def screen_contaminant_kraken_paired(x):
    R, input_dir, bowtie2_indexes, keep_intermediate, rm_hsap, rm_rrna, rm_mmus, rm_pcin, rm_pcoq, rm_mmur, rm_mmul, rm_ptro, rm_sbol, rm_vvar, rm_clup, rm_sscr, rm_eubi, rm_btau, rm_fana, rm_mdom, kraken2_database, confidence_score, minimum_hit_groups, mapper, nprocs_mapper, dry_run, verbose = x
    
    screened = []
    to_removes = []
#   cont_dnas = ['phiX174']
    cont_dnas = []
#    intermediate_files = []
#    kraken_reads, trimmed_reads = [],[]

    if rm_hsap:
        cont_dnas += ['hg38']
        
    if rm_rrna:
        cont_dnas += ['SILVA_138.1_LSURef_NR99_tax_silva_DNA', 'SILVA_138.1_SSURef_NR99_tax_silva_DNA']

    if rm_mmus:
        cont_dnas += ['mmusculus_black6_GCA_000001635_8']
    
    if rm_clup:
        cont_dnas += ['CanFam3.1']
    
    if rm_eubi:
        cont_dnas += ['Eubiome']

    if rm_btau:
        cont_dnas += ['bosTau9-ARS-UCD2.0']
        
    if rm_mdom:
        cont_dnas += ['Mdomestica_GDDH13_1-1']
    
    if rm_fana:
        cont_dnas += ['F_ana_Camarosa_6_28_17']

    outf = ".".join(R[0].split(".")[:-2]) if 'gz' in R[0] or 'bz2' in R[0] else ".".join(R[0].split(".")[:-1])
    Rext = "." + ".".join(R[0].split(".")[-2:]) if 'gz' in R[0] or 'bz2' in R[0] else "." + ".".join(R[0].split(".")[-1:])
    final = None

    for cont_dna in cont_dnas:
        iR1 = "/".join(outf.split("/")[:-1]) + "/" + outf.split("/")[-1][:len(args.samplename)+1] + "R1" + outf.split("/")[-1][len(args.samplename)+3:]
        iR2 = "/".join(outf.split("/")[:-1]) + "/" + outf.split("/")[-1][:len(args.samplename)+1] + "R2" + outf.split("/")[-1][len(args.samplename)+3:]
        suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
        outf += suffix
        outf = "/".join(outf.split("/")[:-1]) + "/" + outf.split("/")[-1][:len(args.samplename)+1] + "R#" + outf.split("/")[-1][len(args.samplename)+3:]
        
        if not (os.path.isfile('{}.fastq'.format(os.path.join(input_dir, outf.replace("#","_1")))) or os.path.isfile('{}.fastq'.format(os.path.join(input_dir, outf.replace("#","_2"))))):
            if "gz" in Rext:
                cmd = ('kraken2 --db {} --confidence {} --minimum-hit-groups {} --threads {} --gzip-compressed --paired --unclassified-out {}_krk2.fastq {} {}'.format(os.path.join(kraken2_database, cont_dna), confidence_score, minimum_hit_groups, nprocs_mapper,
                    os.path.join(input_dir, outf), os.path.join(input_dir, iR1 + Rext), os.path.join(input_dir, iR2 + Rext)))
            else:
                cmd = ('kraken2 --db {} --confidence {} --minimum-hit-groups {} --threads {} --paired --unclassified-out {}_krk2.fastq {} {}'.format(os.path.join(kraken2_database, cont_dna), confidence_score, minimum_hit_groups, nprocs_mapper,
                    os.path.join(input_dir, outf), os.path.join(input_dir, iR1 + Rext), os.path.join(input_dir, iR2 + Rext)))

            if dry_run or verbose:
                info('{}\n'.format(cmd))

            if not dry_run:
                try:
                    with open(os.devnull, 'w') as devnull:
                        sb.check_call(cmd.split(' '), stdout=devnull, stderr=devnull)
                    os.rename('{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","_1"))), '{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","1"))))
                    os.rename('{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","_2"))), '{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","2"))))
                except Exception as e:
                    if os.path.exists('{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","_1")))):
                        os.remove('{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","_1"))))
                    if os.path.exists('{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","_2")))):
                        os.remove('{0}_krk2.fastq'.format(os.path.join(input_dir, outf.replace("#","_2"))))
                

                    error('screen_contaminating_dnas()\n{}\n{}'.format(cmd, e), exit=True)
        
        final = outf + '_krk2.fastq'
        
        rawpycount_bases([["{}_krk2.fastq".format(os.path.join(input_dir, outf.replace("#","1")))],["{}_krk2.fastq".format(os.path.join(input_dir, outf.replace("#","2")))]], args.samplename, output_dir, False, dry_run, verbose, False)

        if not keep_intermediate:
            to_removes.append(iR1 + Rext)
            to_removes.append(iR2 + Rext)
            
        Rext = '_krk2.fastq'
        final = outf + '_krk2.fastq'
#        intermediate_files.append(outf + '_krk2.fastq')
    
    if not dry_run:
        file = open('{0}.stats'.format(os.path.join(input_dir, args.samplename + ".R1_trimmed")), "r")
        file.readline()
        trimmed_reads = float(file.readline().split("\t")[2])
        file.close()

        file = open('{0}_krk2.stats'.format(os.path.join(input_dir, args.samplename + "." + outf.replace("#","1").split("/")[-1][len(args.samplename)+1:])), "r")
        file.readline()
        kraken_reads = float(file.readline().split("\t")[2])
        file.close()

        host = (trimmed_reads - kraken_reads) / trimmed_reads
        print("\nPercentage of host reads: {:.2f}%".format(host*100))

    if host < 0.75 and args.mapper!="stringet":
        if not (os.path.isfile("{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R1"))) or os.path.isfile("{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R2")))):
            write_gz(os.path.join(input_dir, final), dry_run, verbose)
            os.rename("{}_krk2.fastq.gz".format(os.path.join(input_dir, outf.replace("#","1"))), "{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R1")))
            os.rename("{}_krk2.fastq.gz".format(os.path.join(input_dir, outf.replace("#","2"))), "{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R2")))

        if not (os.path.isfile("{}.stats".format(os.path.join(input_dir, args.samplename + "_R1"))) or os.path.isfile("{}.stats".format(os.path.join(input_dir, args.samplename + "_R2")))):
            stats_R1 = pd.read_csv("{}_krk2.stats".format(os.path.join(input_dir, "_".join(outf.split("_")[:-(len(cont_dnas)+2)]) + ".R1_" + "_".join(outf.split("_")[-(len(cont_dnas)+1):]))), sep='\t')
            stats_R1['#samplename'] = [".".join(i.split(".")[:-1]) + "_" + i.split(".")[-1][:2] + ".fastq" for i in stats_R1['#samplename']]
            stats_R1.to_csv("{}.stats".format(os.path.join(input_dir, args.samplename + "_R1")), sep='\t', index=None)
            stats_R2 = pd.read_csv("{}_krk2.stats".format(os.path.join(input_dir, "_".join(outf.split("_")[:-(len(cont_dnas)+2)]) + ".R2_" + "_".join(outf.split("_")[-(len(cont_dnas)+1):]))), sep='\t')
            stats_R2['#samplename'] = [".".join(i.split(".")[:-1]) + "_" + i.split(".")[-1][:2] + ".fastq" for i in stats_R2['#samplename']]
            stats_R2.to_csv("{}.stats".format(os.path.join(input_dir, args.samplename + "_R2")), sep='\t', index=None)
        remove(to_removes, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)
        return (final.replace("#","1"), final.replace("#","2"))
    else:
        cont_dnas = []
        if rm_hsap:
            cont_dnas += ['hg38']
            
        if rm_rrna:
            cont_dnas += ['SILVA_138.1_LSURef_NR99_tax_silva_DNA', 'SILVA_138.1_SSURef_NR99_tax_silva_DNA']
        
        if rm_mmus:
            cont_dnas += ['mmusculus_black6_GCA_000001635_8']
        
        if rm_clup:
            cont_dnas += ['CanFam3.1']

        if rm_eubi:
            cont_dnas += ['Eubiome']

        if rm_btau:
            cont_dnas += ['bosTau9-ARS-UCD2.0']
            
        if rm_mdom:
            cont_dnas += ['Mdomestica_GDDH13_1-1']
        
        if rm_fana:
            cont_dnas += ['F_ana_Camarosa_6_28_17']

        outf1 = ".".join(final.replace("#","1").split(".")[:-2]) if 'gz' in final.replace("#","1") or 'bz2' in final.replace("#","1") else ".".join(final.replace("#","1").split(".")[:-1])
        outf2 = ".".join(final.replace("#","2").split(".")[:-2]) if 'gz' in final.replace("#","2") or 'bz2' in final.replace("#","2") else ".".join(final.replace("#","2").split(".")[:-1])
        Rext = "." + ".".join(final.replace("#","1").split(".")[-2:]) if 'gz' in final.replace("#","1") or 'bz2' in final.replace("#","1") else "." + ".".join(final.replace("#","1").split(".")[-1:])
#         outf1 = final.replace("#","1")[:final.replace("#","1").rfind('.')]
#         outf2 = final.replace("#","2")[:final.replace("#","2").rfind('.')]
#         Rext = final[final.rfind('.'):]
        final = None

        for cont_dna in cont_dnas:
            iR1 = outf1
            iR2 = outf2

            suffix = '_{}'.format(cont_dna.replace('_', '-').replace('.', '-'))
            outf1 += suffix
            outf2 += suffix

            if not (os.path.isfile('{}_bt2.fastq.gz'.format(os.path.join(input_dir, outf1))) or os.path.isfile('{}_bt2.fastq.gz'.format(os.path.join(input_dir, outf2)))):
                        # -x is the argument for indices
                        # --un <path>        write unpaired reads that didn't align to <path>
                        # -p is the cores number
                        # -U input fastq files
                cmd = ('bowtie2 -x {} -1 {} -2 {} -p {} --sensitive-local | SplitUnmappedSAM.py -1 {}_bt2.fastq.gz -2 {}_bt2.fastq.gz'
                        .format(os.path.join(bowtie2_indexes, cont_dna), os.path.join(input_dir, iR1 + Rext), os.path.join(input_dir, iR2 + Rext), nprocs_mapper, os.path.join(input_dir, outf1), os.path.join(input_dir, outf2)))


                if dry_run or verbose:
                    info('{}\n'.format(cmd))

                if not dry_run:
                    try:
                        with open(os.devnull, 'w') as devnull:
                             sb.check_output(cmd, shell=True, stderr=devnull)
    #                         sb.run(cmd.split('|')[1].split(" "), stdin=sb.run(cmd.split('|')[0].split(" "), stdout=sb.PIPE, stderr=devnull))
    #                         sb.Popen(cmd, shell=True,
    #                          stdin=sb.PIPE,
    #                          stdout=sb.PIPE,
    #                          stderr=sb.PIPE)
                    except Exception as e:
                        if os.path.exists('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf1))):
                            os.remove('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf1)))
                        if os.path.exists('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf2))):
                            os.remove('{0}_bt2.fastq.gz'.format(os.path.join(input_dir, outf2)))

                        error('screen_contaminating_dnas()\n{}\n{}'.format(cmd, e), exit=True)

            rawpycount_bases([["{}_bt2.fastq.gz".format(os.path.join(input_dir, outf1))], ["{}_bt2.fastq.gz".format(os.path.join(input_dir, outf2))]], args.samplename, output_dir, False, dry_run, verbose, False)

            if not keep_intermediate:
                to_removes.append(iR1 + Rext)
                to_removes.append(iR2 + Rext)

            Rext = "_bt2.fastq.gz"
    #     final = outf + '_bt2.fastq.bz2'

        if not (os.path.isfile("{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R1"))) or os.path.isfile("{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R2")))):
            os.rename("{}_bt2.fastq.gz".format(os.path.join(input_dir, outf1)), "{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R1")))
            os.rename("{}_bt2.fastq.gz".format(os.path.join(input_dir, outf2)), "{}.fastq.gz".format(os.path.join(input_dir, args.samplename + "_R2")))

        if not (os.path.isfile("{}.stats".format(os.path.join(input_dir, args.samplename + "_R1"))) or os.path.isfile("{}.stats".format(os.path.join(input_dir, args.samplename + "_R2")))):
            stats_R1 = pd.read_csv("{}_bt2.stats".format(os.path.join(input_dir, "_".join(outf1.split("_")[:outf1.split("_").index("trimmed")-1]) + ".R1_" + "_".join(outf1.split("_")[outf1.split("_").index("trimmed"):]))), sep='\t')
            stats_R1['#samplename'] = [".".join(i.split(".")[:-1]) + "_" + i.split(".")[-1][:2] + ".fastq" for i in stats_R1['#samplename']]
            stats_R1.to_csv("{}.stats".format(os.path.join(input_dir, args.samplename + "_R1")), sep='\t', index=None)
            stats_R2 = pd.read_csv("{}_bt2.stats".format(os.path.join(input_dir, "_".join(outf2.split("_")[:outf2.split("_").index("trimmed")-1]) + ".R2_" + "_".join(outf2.split("_")[outf2.split("_").index("trimmed"):]))), sep='\t')
            stats_R2['#samplename'] = [".".join(i.split(".")[:-1]) + "_" + i.split(".")[-1][:2] + ".fastq" for i in stats_R2['#samplename']]
            stats_R2.to_csv("{}.stats".format(os.path.join(input_dir, args.samplename + "_R2")), sep='\t', index=None)

        remove(to_removes, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)

        return (outf1 + Rext, outf2 + Rext)

def screen_contaminating_dnas_mp(x):

    if args.mapper=='bowtie2' :
        return screen_contaminant_bowtie2(x)

    else:
        return screen_contaminant_kraken(x)
                                                                                                                          
def screen_contaminating_dnas_paired(x):

    if args.mapper=='bowtie2' :
        return screen_contaminant_bowtie2_paired(x)

    else:
        return screen_contaminant_kraken_paired(x)

def get_unpaired(input_dir, r1_r2, samplename):
    tasks = zip(r1_r2, [input_dir] * len(r1_r2))

    r1_index, r2_index = pool_get_unpaired(((get_unpaired_mp, args) for args in tasks))

    unpaired_file = samplename + '_unpaired.txt.bz2'

    with bz2.open(os.path.join(input_dir, unpaired_file), 'wt') as f:
        sym_difference = set(r1_index).symmetric_difference(set(r2_index))
        f.write('\n'.join( (i for i in sym_difference) ) + '\n')
#        f.write('\n'.join( (i for i in r1_index if i not in r2_index) ) + '\n' +
#                '\n'.join( (i for i in r2_index if i not in r1_index) ) + '\n')

    return unpaired_file

def get_unpaired_mp(x):
    R, input_dir = x
    r_index = SeqIO.index(os.path.join(input_dir, R), "fastq")
    return r_index

def pool_get_unpaired(arg):
    terminating = mp.Event()

    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, arg, chunksize=2)]
        except Exception as e:
            error('get_unpaired()\ntasks: {}\n    e: {}'.format(arg, e), init_new_line=True, exit=True)

def split_and_sort(input_dir, screened_r1_r2, samplename, keep_intermediate, unpaired_file, 
                   nproc=1, dry_run=False, verbose=False):
    if dry_run or verbose:
        info('split_and_sort()\n', init_new_line=True)

    R1, R2 = screened_r1_r2

    if not (os.path.isfile(samplename + '_R1.fastq.bz2') and
            os.path.isfile(samplename + '_R2.fastq.bz2') and
            os.path.isfile(samplename + '_UN.fastq.bz2')):
        ## Changed on 23/04/24
        cmd = 'split_and_sort.new.update.test2.py --R1 {} --R2 {} --prefix {}'.format(os.path.join(input_dir, R1),
                                                                         os.path.join(input_dir, R2),
                                                                         os.path.join(input_dir, samplename))

        if unpaired_file and os.path.isfile(os.path.join(input_dir, unpaired_file)):
            cmd += ' --unpaired {}'.format(os.path.join(input_dir, unpaired_file))

        if dry_run or verbose:
            info('{}\n'.format(cmd))

        if not dry_run:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                for i in [os.path.isfile(samplename + '_R1.fastq.bz2'),
                          os.path.isfile(samplename + '_R2.fastq.bz2'),
                          os.path.isfile(samplename + '_UN.fastq.bz2')]:
                    if os.path.exists(i):
                        os.remove(i)

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    if not os.path.isfile(input_dir + "/" + samplename + '_summary.stats'):
        print(input_dir + samplename + '_summary.stats')
        cmd = 'cat_stats.py -i {} -o {}'.format(input_dir, os.path.join(input_dir, samplename + '_summary.stats'))
        
        if dry_run or verbose:
            info('{}\n'.format(cmd))

        if not dry_run:
            try:
                sb.check_call(cmd.split(' '))
            except Exception as e:
                if os.path.exists(samplename + '_summary.stats'):
                    os.remove(samplename + '_summary.stats')

                error('split_and_sort()\n{}\n{}'.format(cmd, e), exit=True)

    remove(screened_r1_r2, keep_intermediate, folder=input_dir, dry_run=dry_run, verbose=verbose)
    return (samplename + '_R1.fastq.bz2', samplename + '_R2.fastq.bz2', samplename + '_UN.fastq.bz2')

def pool_cat_stats(arg):
    terminating = mp.Event()

    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, arg, chunksize=2)]
        except Exception as e:
            error('cat_stats()\ntasks: {}\n    e: {}'.format(arg, e), init_new_line=True, exit=True)

def cat_stats_mp(x):
    if not terminating.is_set():
        try:
            output = x

        except Exception as e:
            terminating.set()

            if os.path.exists(output):
                os.remove(output)

            error('cat_stats_mp()\n    x: {}\n    e: {}'.format(x, e), init_new_line=True)
            raise
    else:
        terminating.set()


def remove(to_remove, keep_intermediate, folder=None, dry_run=False, verbose=False):
    if verbose:
        info('remove()\n', init_new_line=True)

    if not keep_intermediate:
        for r in to_remove:
            rf = os.path.join(folder, r) if folder else r

            if os.path.isfile(rf):
                if verbose:
                    info('rm {}\n'.format(rf))

                if not dry_run:
                    os.remove(rf)


if __name__ == "__main__":
    t0 = time.time()
    args = read_params()

    if args.verbose:
        info('Preprocessing version {} ({})\n'.format(__version__, __date__))
        info('Command line: {}\n'.format(' '.join(sys.argv)), init_new_line=True)

    check_params(args)
    #preflight_check(dry_run=args.dry_run, verbose=args.verbose)
    
    output_dir = args.output_dir
    
    os.makedirs(output_dir, exist_ok=True)
    
    inputs_r1s_r2s = get_inputs(args.input_dir, args.forward, args.reverse, args.samplename, args.extension, verbose=args.verbose)

#    if args.mapper=='kraken' and (not args.rm_hsap and not args.rm_mmus):
#        error('Using kraken2, the host should be human or mouse, Please set one of --rm_mmus or --rm_hsap as True', exit=True)

    if (len(inputs_r1s_r2s[0]) == 0) or (len(inputs_r1s_r2s[1]) == 0):
        error('No input files detected!\nR1s: {}\nR2s: {}'.format(inputs_r1s_r2s[0], inputs_r1s_r2s[1]), exit=True)

    if args.dry_run or args.verbose:
        info('inputs_r1s: {}\n'.format('\n            '.join(inputs_r1s_r2s[0])), init_new_line=True)
        info('inputs_r2s: {}\n'.format('\n            '.join(inputs_r1s_r2s[1])))

#    merged_r1_r2 = concatenate_reads(args.input_dir, inputs_r1s_r2s, nproc=args.nproc, dry_run=args.dry_run, verbose=args.verbose)
    rawpycount_bases(inputs_r1s_r2s, args.samplename, output_dir, False, args.dry_run, args.verbose, False)

    qced_r1_r2, unpaired_file = quality_control(args.input_dir, inputs_r1s_r2s, output_dir, args.keep_intermediate, args.samplename,
                                 nproc=args.nproc, dry_run=args.dry_run, verbose=args.verbose)
#    remove(merged_r1_r2, args.keep_intermediate, folder=args.input_dir, dry_run=args.dry_run, verbose=args.verbose)

    if args.dry_run or args.verbose:
        info('qced_r1: {}\n'.format(qced_r1_r2[0]), init_new_line=True)
        info('qced_r2: {}\n'.format(qced_r1_r2[1]))

#    unpaired_file = get_unpaired(args.input_dir, qced_r1_r2, args.samplename) if args.paired_end else None

    if any([args.rm_hsap,args.rm_rrna,args.rm_mmus,args.rm_pcin,args.rm_pcoq,args.rm_mmur,args.rm_mmul,args.rm_ptro,args.rm_sbol,args.rm_vvar,args.rm_clup,args.rm_sscr,args.rm_eubi,args.rm_btau,args.rm_mdom,args.rm_fana]):
        screened_r1_r2 = screen_contaminating_dnas(output_dir, qced_r1_r2, args.bowtie2_indexes, args.keep_intermediate,
                                                args.rm_hsap, args.rm_rrna, args.rm_mmus, args.rm_pcin, args.rm_pcoq, args.rm_mmur,
                                                args.rm_mmul, args.rm_ptro, args.rm_sbol, args.rm_vvar, args.rm_clup, args.rm_sscr,
                                                args.rm_eubi, args.rm_btau, args.rm_mdom, args.rm_fana, args.kraken2_database, args.run_confidence, args.run_minimum_hit_groups,
                                                args.mapper, nprocs_mapper=args.nproc_mapper if args.nproc_mapper > args.nproc else args.nproc, dry_run=args.dry_run, verbose=args.verbose)
        remove(qced_r1_r2, args.keep_intermediate, folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)

        if args.dry_run or args.verbose:
            info('screened_r1: {}\n'.format(screened_r1_r2[0]), init_new_line=True)
            info('screened_r2: {}\n'.format(screened_r1_r2[1]))

        if not args.paired_end:
            splitted_and_sorted = split_and_sort(output_dir, screened_r1_r2, args.samplename, args.keep_intermediate, unpaired_file,
                                             nproc=args.nproc, dry_run=args.dry_run, verbose=args.verbose)
            remove(screened_r1_r2, args.keep_intermediate, output_dir, dry_run=args.dry_run, verbose=args.verbose)
            remove([f for f in os.listdir(output_dir) if f.endswith('.stats') and not f.endswith('_summary.stats')], args.keep_intermediate, 
               folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)
            if args.dry_run or args.verbose:
                info('splitted_and_sorted: {}\n'.format(splitted_and_sorted[0]), init_new_line=True)
                info('                     {}\n'.format(splitted_and_sorted[1]))
                info('                     {}\n'.format(splitted_and_sorted[2]))
        else:
            remove(screened_r1_r2, args.keep_intermediate, folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)
            if not os.path.isfile(output_dir + "/" + args.samplename + '_summary.stats'):
                cmd = 'cat_stats_v2.py -i {} -o {}'.format(output_dir + "/", output_dir + "/" + args.samplename + '_summary.stats')

                if args.dry_run or args.verbose:
                    info('{}\n'.format(cmd))

                if not args.dry_run:
                    try:
                        sb.check_call(cmd.split(' '))
                    except Exception as e:
                        if os.path.exists(output_dir + "/" + output_dir + '_summary.stats'):
                            os.remove(output_dir + "/" + output_dir + '_summary.stats')
            remove([f for f in os.listdir(output_dir) if f.endswith('.stats') and not f.endswith('_summary.stats')], args.keep_intermediate, 
               folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)
    else:
        if not args.paired_end:
            splitted_and_sorted = split_and_sort(output_dir, qced_r1_r2, args.samplename, args.keep_intermediate, unpaired_file,
                                         nproc=args.nproc, dry_run=args.dry_run, verbose=args.verbose)
            remove(qced_r1_r2, args.keep_intermediate, folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)
            remove([f for f in os.listdir(output_dir) if f.endswith('.stats') and not f.endswith('_summary.stats')], args.keep_intermediate,
                folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)
            if args.dry_run or args.verbose:
                info('splitted_and_sorted: {}\n'.format(splitted_and_sorted[0]), init_new_line=True)
                info('                     {}\n'.format(splitted_and_sorted[1]))
                info('                     {}\n'.format(splitted_and_sorted[2]))
        else:
            if not (os.path.isfile("{}.fastq.gz".format(os.path.join(output_dir, args.samplename + "_R1"))) or os.path.isfile("{}.fastq.gz".format(os.path.join(output_dir, args.samplename + "_R2")))):
                for R in qced_r1_r2:
                    os.rename("{}".format(os.path.join(output_dir, R)), "{}.fastq.gz".format(os.path.join(output_dir, R.split("_trimmed")[0])))

            if not (os.path.isfile("{}.stats".format(os.path.join(output_dir, args.samplename + "_R1"))) or os.path.isfile("{}.stats".format(os.path.join(output_dir, args.samplename + "_R2")))):
                stats_R1 = pd.read_csv("{}.R1_trimmed.stats".format(os.path.join(output_dir, args.samplename)), sep='\t')
                stats_R1['#samplename'] = args.samplename + "_R1.fastq"
                stats_R1.to_csv("{}.stats".format(os.path.join(output_dir, args.samplename + "_R1")), sep='\t', index=None)
                stats_R2 = pd.read_csv("{}.R2_trimmed.stats".format(os.path.join(output_dir, args.samplename)), sep='\t')
                stats_R2['#samplename'] = args.samplename + "_R2.fastq"
                stats_R2.to_csv("{}.stats".format(os.path.join(output_dir, args.samplename + "_R2")), sep='\t', index=None)
            remove(qced_r1_r2, args.keep_intermediate, folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)
            if not os.path.isfile(output_dir + "/" + args.samplename + '_summary.stats'):
                cmd = 'cat_stats_v2.py -i {} -o {}'.format(output_dir + "/", output_dir + "/" + args.samplename + '_summary.stats')

                if args.dry_run or args.verbose:
                    info('{}\n'.format(cmd))

                if not args.dry_run:
                    try:
                        sb.check_call(cmd.split(' '))
                    except Exception as e:
                        if os.path.exists(output_dir + "/" + output_dir + '_summary.stats'):
                            os.remove(output_dir + "/" + output_dir + '_summary.stats')
            remove([f for f in os.listdir(output_dir) if f.endswith('.stats') and not f.endswith('_summary.stats')], args.keep_intermediate,
               folder=output_dir, dry_run=args.dry_run, verbose=args.verbose)


    if args.verbose:
        info('time elapsed: {} s\n'.format(int(time.time() - t0)), init_new_line=True)

    sys.exit(0)
