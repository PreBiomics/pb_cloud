#!/mnt/Prebiomics_Data/tools/anaconda3/envs/mpa4.1/bin/python

from Bio import SeqIO
import numpy as np
from datetime import datetime
from multiprocessing import Event
import traceback
from concurrent import futures
import argparse
import bz2
import gzip

def init_terminating(terminating_):
    global terminating
    terminating = terminating_

def parallel_execution(arguments):
    function, *args = arguments
    if not terminating.is_set():
        try:
            return function(*args)
        except Exception as e:
            terminating.set()
            traceback.print_exc()
            print("Error executing {}".format(arguments))
            print("Parallel execution fails: "+str(e))
    else:
        terminating.set()

def execute_pool(args, nprocs):
    terminating = Event()
    with futures.ThreadPoolExecutor(initializer=init_terminating, initargs=(terminating,), max_workers=2) as pool:
        try:
            return [_ for _ in pool.map(parallel_execution, args, chunksize=2)]
        except Exception as e:
            print('Parallel execution fails: '+str(e))

def parallel_read_unpaired(input_file):
    if not 'gz' in input_file:
        return SeqIO.index(input_file, format='fastq')
    else:
        with gzip.open(input_file,'rt') as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, format="fastq"))

def parallel_process_unpaired_sorted(to_check, where_to_check, reads_index, out, samplename):
    unpair_reads, lenn = [], []
    for x in to_check:
        if x in where_to_check:
            SeqIO.write(reads_index[x], out, 'fastq')
            lenn.append(len(reads_index[x].seq))
        elif x in unpaired_reads:
            unpair_reads.append(reads_index[x])
    with open(samplename+".stats", 'w') as outf:
        outf.write('\t'.join(["#samplename", "n_of_bases", "n_of_reads", "min_read_len", "median_read_len", "mean_read_len", "max_read_len"]) + '\n')
        outf.write('\t'.join([str(a) for a in [samplename.split("/")[-1]+".fastq", np.sum(lenn), len(lenn), np.min(lenn), np.median(lenn), np.mean(lenn), np.max(lenn)]]) + '\n')
    return unpair_reads

init = datetime.now()

parser = argparse.ArgumentParser()
parser.add_argument("--R1", help="R1 in", required=True)
parser.add_argument("--R2", help="R2 in", required=True)
parser.add_argument("-p", '--prefix', help="Destination prefix", default='out')
parser.add_argument("-u", '--unpaired', help="File with unpaired read IDs")
args = parser.parse_args()

filR1 = bz2.open(args.prefix + '_R1.fastq.bz2', 'wt')
filR2 = bz2.open(args.prefix + '_R2.fastq.bz2', 'wt')
filUN = bz2.open(args.prefix + '_UN.fastq.bz2', 'wt')

global unpaired_reads
unpaired_reads = set([i.strip() for i in bz2.open(args.unpaired, 'rt')])

R1_sortedIndex, R2_sortedIndex = execute_pool(((parallel_read_unpaired, readfile) for readfile in [args.R1, args.R2]), 2)

i1 = sorted(R1_sortedIndex.keys()).copy()
c1 = set((R2_sortedIndex.keys()))

i2 = sorted(R2_sortedIndex.keys()).copy()
c2 = set(R1_sortedIndex.keys())

unpaired = execute_pool(((parallel_process_unpaired_sorted, *args) for args in [(i1, c1, R1_sortedIndex, filR1, args.prefix + '_R1'), (i2, c2, R2_sortedIndex , filR2, args.prefix + '_R2')]), 2)

if len(unpaired[0])>0 and len(unpaired[1])>0:
    lenn = []
    for pair in unpaired:
        for read in pair:
            lenn.append(len(read.seq))
            SeqIO.write(read, filUN, 'fastq')

    samplename = args.prefix + '_UN'
    with open(samplename+".stats", 'w') as outf:
            outf.write('\t'.join(["#samplename", "n_of_bases", "n_of_reads", "min_read_len", "median_read_len", "mean_read_len", "max_read_len"]) + '\n')
            outf.write('\t'.join([str(a) for a in [samplename.split("/")[-1]+".fastq", np.sum(lenn), len(lenn), np.min(lenn), np.median(lenn), np.mean(lenn), np.max(lenn)]]) + '\n')
    filUN.close()
        
filR1.close()
filR2.close()

print('Computation time: ', datetime.now()-init)
