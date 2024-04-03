#!/usr/bin/env python3

import argparse
from Bio import SeqIO, pairwise2
import sys
import multiprocessing

parser = argparse.ArgumentParser(
        description='''orientate - use the beginning of a ref FASTA to orient 
            outcomes FASTA records to start in same approx place'''
        )

parser.add_argument('-r','--ref',help='what to orient input to')
parser.add_argument('-n','--number-bases',help='how many bases to use in aln',
        default=50)
parser.add_argument('-c','--cpus',help='how many parallel jobs to launch',
        type=int,default=1)
parser.add_argument('-i','--input',default="-",help='input to orient, default stdin')
parser.add_argument('-o','--output',default="-",help='outputs, default stdout')
args = parser.parse_args()

ref = SeqIO.parse(args.ref,'fasta').__next__()
ref_header = ref[0:args.number_bases]


if args.input == '-':
    inputs = sys.stdin
else:
    inputs = open(args.input,"r")


if args.output == '-':
    outputs_go = sys.stdout
else:
    outputs_go = open(args.output,"wa")

def align2ref(query,ref):

    forward_aln = pairwise2.align.localxs(
            query.seq,
            ref.seq,
            -1, #open   I think?
            -1, #extend I think?
            one_alignment_only=True
            )

    reverse_aln = pairwise2.align.localxs(
            query.reverse_complement().seq,
            ref.seq,
            -1, #open   I think?
            -1, #extend I think?
            one_alignment_only=True
            )

    if forward_aln[0].score >= reverse_aln[0].score:
        outseq = query
        the_output = ">"+outseq.id+"\n"+ \
                str( (outseq[forward_aln[0].start:len(outseq)] + 
                    outseq[0:forward_aln[0].start]).seq) +"\n" 
    else:
        outseq = query
        outseq.seq = query.reverse_complement().seq
        the_output = ">"+outseq.id+"\n"+ \
                str( (outseq[reverse_aln[0].start:len(outseq)] + 
                    outseq[0:reverse_aln[0].start]).seq) +"\n" 

    return the_output

with multiprocessing.Pool(args.cpus) as mp:

    for i in mp.starmap(align2ref,[ (k,ref_header) for k in SeqIO.parse(inputs,'fasta') ]):
        outputs_go.write(i)

inputs.close()

