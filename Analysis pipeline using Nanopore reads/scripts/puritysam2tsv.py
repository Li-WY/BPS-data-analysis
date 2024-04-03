#!/usr/bin/env python3

import os
import sys
import multiprocessing
import pysam
import argparse

parser = argparse.ArgumentParser("""
        This takes a "purity-align" SAM file (so each read aligned back to a 
        reference), processes it, and summarizes how each read matches.
        Written in here because it was blowing out the RAM in R.
        """)

parser.add_argument('-i','--input-sam')

args = parser.parse_args()

infile = pysam.AlignmentFile(args.input_sam,'r')

def cs_parser(string):
    # state switching parser, so start in None state
    state = None
    counter = 0
    assemble_number = ''
    stats = { 'bases': {'=': 0,':': 0,'*': 0,'+': 0,'-': 0,'~': 0 } ,
            'states': {'=': 0,':': 0,'*': 0,'+': 0,'-': 0,'~': 0 } }
    # just initialized stuff to hold things, then start going through string
    for i in string:
        if i in ['=',':','*','+','-','~']:
            # then it's a mode change
            if state is not None:
                # check that we're not at the very start, ie we have a populated
                # counter or assemble_number
                if state == '*':
                    # if it was previously a * then it should be two bases 
                    # showing change, but is only ever one mismatch
                    stats['bases'][state] += 1
                if state in [':','~']:
                    # if it was previously one of these, then we should have
                    # an assemble_number populated up there, we should be
                    # able to interpret that as a counter and reset
                    counter = int(assemble_number)
                    assemble_number = ''
                else:
                    # increment by the counter and reset
                    stats['bases'][state] += counter
                    counter = 0
            # switch states
            state = i
            # and record starting a state
            stats['states'][state] += 1
        else:
            if state == '*':
                # if so, then skip because it's only ever two letters
                next
            elif state in [':','~'] :
                try: 
                    # if here, then try and see if it's numeric (int)
                    assemble_number += str(int(i))
                except:
                    # else its a base, so keep going
                    next
            else:
                # otherwise then just increment it along
                counter += 1
    # at end, add on counter to the started state
    stats['bases'][state] += counter
    return(stats)

def get_cigar_out(sam_entry):
    try:
        return( [ { (a): b for a,b in 
                zip(['M',' I',' D',' N',' S',' H',' P',' =',' X',' B',' NM'],i) 
                } for i in sam_entry.get_cigar_stats() ]
            )
    except:
        return(None)

def get_cs_out(sam_entry):
    # here trying two tag names, incase case changes
    try:
        return( cs_parser(sam_entry.get_tag('cs')) )
    except:
        return( cs_parser(sam_entry.get_tag('CS')) )

def purity_parser(sam_entry):
    return('\t'.join( [
            sam_entry.query_name,
            sam_entry.reference_name,
            str(get_cs_out(sam_entry)['bases']['=']),
            str(sam_entry.query_length),
            str(sam_entry.reference_length)
            ] )
        )
# 230731 fyi pairwise aln uses format 
#     ('id','ref','aln_score','length','payload_seq','ref_seq') 

for i in infile:
    if not i.is_unmapped and not i.is_secondary and not i.is_supplementary:
        try:
            print(purity_parser(i))
        except:
            print("Welp record "+i.query_name+" did not process nice",
                file=sys.stderr)

