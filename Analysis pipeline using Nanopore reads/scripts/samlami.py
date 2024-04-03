#!/usr/bin/env python3

import argparse
import pysam
import sys

parser = argparse.ArgumentParser(
        description='samlami-slicer - cut up a SAM to left or right of alignment'
        )

xorgroup = parser.add_mutually_exclusive_group()
xorgroup.add_argument('-l','--cut-on-left',action="store_true",default=False,
    help='DEFAULT cut the left side of aligment off, and the alignment too')
xorgroup.add_argument('-r','--cut-on-right',action="store_true",default=False,
    help='cut the right side of aligment off, and the alignment too')
xorgroup.add_argument('-m','--modulus',action="store_true",default=False,
    help='orient the sequence so that it starts with the aligned part')

parser.add_argument('-i','--input-sam',default="-",
    help='where read')
parser.add_argument('-o','--output-sam',default="-",
    help='where out')

args = parser.parse_args()

infile  = pysam.AlignmentFile(args.input_sam,"r")
outfile = pysam.AlignmentFile(args.output_sam,"w",template=infile)

def cut_either(record, which):
    cigar_tuples = record.cigartuples

    if which == "modulus":
        start_clip = cigar_tuples[0]
        if start_clip[0] == 4 or start_clip[0] == 5:
            start_pos = start_clip[1]
        elif start_clip[0] == 0:
            start_pos = 0
        record.query_sequence  = record.query_sequence[start_pos:]+\
                                    record.query_sequence[:start_pos]
        try:
            tmp_qualities = record.query_qualities
            record.query_qualities = tmp_qualities[start_pos:]+\
                                    tmp_qualities[:start_pos]
        except:
            record.query_qualities = [40]*len(record.query_sequence)
        record.cigartuples = tuple()
        record.flag = 0
        if len(record.seq) > 1:
            return record

    if which == "left":
        the_op = cigar_tuples[len(cigar_tuples)-1]
        if the_op[0] == 4 or the_op[0] == 5:
            start_pos = len(record.seq) - the_op[1]
            end_pos = -1
    elif which == "right":
        the_op = cigar_tuples[0]
        if the_op[0] == 4 or the_op[0] == 5:
            start_pos = 0
            end_pos = the_op[1]

    try:
        tmp_qualities = record.query_qualities
        record.query_sequence  = record.query_sequence[start_pos:end_pos]
        record.query_qualities = tmp_qualities[start_pos:end_pos]
        record.cigartuples = tuple()
        if len(record.seq) > 1:
            return record
    except:
        raise Exception("uh... the leading or lagging operation is not a "+
                "soft or hard clip")

for i in infile:
    if not i.is_unmapped and not i.is_secondary and not i.is_supplementary:
        try:
            if args.cut_on_left:
                outfile.write(cut_either(i,"left"))
            elif args.cut_on_right:
                outfile.write(cut_either(i,"right"))
            elif args.modulus:
                outfile.write(cut_either(i,"modulus"))
        except:
            print("Welp record "+i.query_name+" did not play nice",
                file=sys.stderr)
