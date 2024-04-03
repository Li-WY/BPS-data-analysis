#!/bin/env/python3

from Bio import SeqIO as si
from Bio import pairwise2
import Bio.Seq
import Bio.SeqRecord
import re
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--pipe-quads',action='store_true',default=False)
parser.add_argument('--query',type=str)
parser.add_argument('--ref',type=str)
parser.add_argument('--query-string',type=str)
parser.add_argument('--ref-string',type=str)
parser.add_argument('--score-only',action='store_true')
parser.add_argument('--output',default=sys.stdout,type=str)
parser.add_argument('--verbose','-v',action='count',default=0)
args = parser.parse_args()

def read_seq_file(path,verbosity=args.verbose):
    if re.search("\.fasta$",path) is not None:
        the_format = "fasta"
    elif re.search("\.fastq$",path) is not None:
        the_format = "fastq"

    if verbosity >= 1:
        print("Parsing",path,"as format",the_format)

    return( si.parse(path,the_format) )


def read_seq_string(string,separator=',',verbosity=args.verbose):
    parts = string.split(separator)
    if verbosity >= 1:
        print("Putting string",string,"into FASTA format")
    return( Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(parts[1],'fasta'),id=parts[0]))

if args.pipe_quads:

    for i in sys.stdin:
        (qid, qseq, rid, rseq) = i.split("\t")
        query = read_seq_string(qid+","+qseq)
        reference = read_seq_string(rid+","+rseq)

        if args.score_only:

            this_alignment = pairwise2.align.localxx(
                    query.seq, reference.seq, 
                    #1,-1,-1,-1,
                    one_alignment_only=True,
                    score_only=True,
                    penalize_end_gaps=False)

            if args.verbose >= 1:
                sys.stderr.write(query.id+"\t"+reference.id)

            sys.stdout.write(
                    query.id+"\t"+
                    reference.id+"\t"+
                    str(this_alignment)+"\t"+
                    str(len(query.seq))+"\n"
                    )

        else:

            this_alignment = pairwise2.align.localxx(
                    query.seq, reference.seq, 
                    #1,0,0,0,
                    one_alignment_only=True,
                    penalize_end_gaps=False)

            if args.verbose >= 2:
                print( pairwise2.format_alignment(*this_alignment[0]))
            #print(this_alignment[0])

            if args.verbose >= 1:
                sys.stderr.write(query.id+"\t"+reference.id)

            sys.stdout.write(
                    query.id+"\t"+
                    reference.id+"\t"+
                    str(this_alignment[0][2])+"\t"+
                    str(len(query.seq))+"\t"+
                    str(this_alignment[0][0])+"\t"+
                    str(this_alignment[0][1])+"\n"
                    )

else:

    if args.ref is not None:
        reference = read_seq_file(args.ref)
    elif args.ref_string is not None:
        reference = [read_seq_string(args.ref_string),]

    if args.query is not None:
        query = read_seq_file(args.query)
    elif args.query_string is not None:
        query = [read_seq_string(args.query_string),]

    for this_ref in reference:
        for this_query in query:
            this_alignment = pairwise2.align.localxx(
                    this_query.seq, this_ref.seq, 
                    #1,0,0,0,
                    penalize_end_gaps=False)
            if args.verbose >= 2:
                print( pairwise2.format_alignment(*this_alignment[0]))
            #print(this_alignment[0])
            sys.stdout.write(
                    this_query.id+"\t"+
                    this_ref.id+"\t"+
                    str(this_alignment[0][2])+"\t"+
                    str(len(this_query.seq))+"\t"+
                    str(this_alignment[0][0])+"\t"+
                    str(this_alignment[0][1])+"\n"
                    )
