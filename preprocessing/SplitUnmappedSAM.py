#!/usr/bin/env python

import sys, bz2, gzip, time
import argparse as ap


def read_params():
    p = ap.ArgumentParser(description="")                   
    p.add_argument('-1', '--R1', type=str, default=None,
                   help="The output R1")                   
    p.add_argument('-2', '--R2', type=str, default=None,
                   help="The output R2")   
    return p.parse_args()

def check_params(args):
    if not args.R1:
        print('ERROR: -1 (or --R1) must be specified')
        exit(1)
    if not args.R2:
        print('ERROR: -2 (or --R2) must be specified')
        exit(1)

def split_reads(file_r1, file_r2):
    r1 = None
    r2 = None

    w1 =  bz2.open(file_r1, 'wt') if file_r1.endswith('.bz2') else gzip.open(file_r1, 'wt') if file_r1.endswith('.gz') else open(file_r1, 'w')
    w2 =  bz2.open(file_r2, 'wt') if file_r2.endswith('.bz2') else gzip.open(file_r2, 'wt') if file_r2.endswith('.gz') else open(file_r2, 'w')

    for line in sys.stdin:
        if not line.startswith('@'):
            if r1 is None:
                r1 = line.strip().split('\t')
            elif r2 is None:
                r2 = line.strip().split('\t')

            if r1 is not None and r2 is not None:
                flag1 = "{0:b}".format(int(r1[1]))
                flag2 = "{0:b}".format(int(r2[1]))
                
                if r1[0] != r2[0]:
                    print('ERROR: Something went wrong, headers in the SAM are not paired')
                    w1.close()
                    w2.close()
                    exit(1)

                if flag1[-3] == flag2[-3] == '1':
                    w1.write('@' + r1[0] + '\n')
                    w1.write(r1[9]+'\n')
                    w1.write('+\n')
                    w1.write(r1[10]+'\n')
                    w2.write('@' + r2[0]+'\n')
                    w2.write(r2[9]+'\n')
                    w2.write('+\n')
                    w2.write(r2[10]+'\n')
                r1 = None
                r2 = None

    w1.close()
    w2.close()


if __name__ == "__main__":
    print("Start execution")
    t0 = time.time()
    args = read_params()
    check_params(args)
    split_reads(args.R1, args.R2)
    exec_time = time.time() - t0
    print("Finish execution ({}) seconds)".format(round(exec_time, 2)))
