#!/bin/env/python3

#import Bio as bp
from Bio import SeqIO as si
from Bio import AlignIO as ai
import re
#import glob

import pandas as pd

import sys
#import multiprocessing
#import itertools

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fasta',default=sys.stdin,type=str)
parser.add_argument('--type',choices=['nanopore'],#,'pacbio'],
    default='nanopore',type=str)
parser.add_argument('--verbose','-v',action='count',default=0)
parser.add_argument('--output',default=sys.stdout,type=str)
#parser.add_argument('--cores',default=None,type=int)
parser.add_argument('--delimiter-between-position-and-votes',type=str,default='-')
parser.add_argument('--rq_cutoff',default=None,type=float)
args = parser.parse_args()

def parse_fasta_seqrecord(x,type="nanopore",rq_cutoff=None):
    if type == "nanopore":
        try:
#            ( id, code ) = re.split('_',x.id)
            return {
                'id':x.id,
#                'code':code,
                'seq':list(str(x.seq))
                }
        except:
            return None
#    elif type == "pacbio":
#        try:
#            ( id1,id2,id3,id4,id5, num_passes, quality, multi, 
#                orf, code, cigar ) = re.split('_',x.id)
#            if rq_cutoff is not None:
#                if float(quality) < rq_cutoff:
#                    return None
#            return {
#                'id':x.id,
#                'np':num_passes,
#                'rq':quality,
#                'multi':multi,
#                'orf':orf,
#                'cigar':cigar,
#                'code':code,
#                'seq':list(str(x.seq))
#                }
#        except:
#            return None
    else:
        return None

def tabulate_bases_atcg(bases):
    return [
        bases.count('A'),
        bases.count('T'),
        bases.count('C'),
        bases.count('G'),
        bases.count('-')
        ]

def vote_along_seqs(seqs,weights=None): #,np):
    if weights is None: # Then equally weigh each sequence
        weights = pd.Series([1]*len(seqs))
    def vote(*bases):
        dot_prodz = \
            [ [i*k for k in j] 
                for i,j in zip(weights, map(tabulate_bases_atcg, bases)) 
                ]
        votes = [ sum(map(lambda x: x[i], dot_prodz)) for i in [0,1,2,3,4] ]
        return votes
    votes = list(map(vote,*seqs))
    this_seq = list(map(lambda x: ('A','T','C','G','-')[x.index(max(x))],votes))
    which_seq = list( map(
                lambda x,y:                 # i is the base position
                    [ seqs.iloc[i][x] == y for i in range(len(seqs)) ], 
                range(len(seqs.iloc[0])),   # x is which seq
                this_seq                    # y is which base was voted
                ) )
#    nps = list( map(  # from doing this for compiling number passes from pacbio
#                lambda x:
#                    (   sum([ np.iloc[i] for i,z in enumerate(x) if z]), 
#                        sum([ np.iloc[i] for i,z in enumerate(x) if not z]) ),
#                which_seq
#        )   )
#    more_than_np = [ x[0]-x[1] for x in nps ]
    rqs = list( map( 
                lambda x:
                    (   sum([ weights.iloc[i] for i,z in enumerate(x) if z]), 
                        sum([ weights.iloc[i] for i,z in enumerate(x) if not z]) ),
                which_seq
        )   )
    more_than_rq = [ x[0]-x[1] for x in rqs ]
    return ( 
        re.sub("-","","".join(this_seq)), 
# TODO does this mean the qualities are offset from the removed gaps in the seq?
        #nps[more_than_np.index(min(more_than_np))],
        rqs[more_than_rq.index(min(more_than_rq))],
        )

# for use later, for voting ? if need be?
#with multiprocessing.Pool(args.cores) as p:
#    p.starmap(vote_and_write_pdg,
#        [ i+(output_prefix,args.type) for i in meta_data.groupby('code') ]
#        )

def vote_and_write_pdg(code,group,output,type='nanopore'):
    if type == 'nanopore':
        try:
            this_seq_result = vote_along_seqs(
                    #group.loc[:,'rq'].astype("float") ,
                    group.loc[:,'seq'], 
                    #group.loc[:,'np'].astype("int")
                    )
            if args.verbose >= 3:
                print("Resulting in:")
                print(this_seq_result)
            with open(output,"w") as f:
                f.write(
                    ">"+code+args.delimiter_between_position_and_votes+
                    str(this_seq_result[1][0])+"to"+str(this_seq_result[1][1])+
                    "\n"+
                    this_seq_result[0]+"\n"
                    )
            if args.verbose >= 2:
                print('Printed out to',output)
        except:
            if args.verbose >= 2:
                try:
                    print('Found error with'+
                        ">"+code+args.delimiter_between_position_and_votes+
                        str(this_seq_result[1][0])+"to"+str(this_seq_result[1][1])
                        )
                except:
                    print('Found double error with'+code)

#    elif type == 'pacbio':
#        for multi, subgroup in group.groupby('multi'):
#            this_seq_result = vote_along_seqs(
#                    subgroup.loc[:,'seq'], 
#                    subgroup.loc[:,'rq'].astype("float") ,
#                    #subgroup.loc[:,'np'].astype("int")
#                    )
#            with open(output,"w") as f:
#                f.write(
#                    ">"+code+"_"+
#                    str(this_seq_result[1][0])+"to"+str(this_seq_result[1][1])+"_"+
#                    str(multi)+"\n"+
#                    this_seq_result[0]+"\n"
#                    )
#            if args.verbose >= 2:
#                print('Printed out to',output)


if args.verbose >= 1:
    print("Parsing",args.fasta)

these_fastas = si.parse(args.fasta,"fasta")

list_of_seqs_per_code = list()
for j in these_fastas:
    parsed = parse_fasta_seqrecord(j,rq_cutoff=args.rq_cutoff)
    if parsed is not None:
        parsed['code'] = re.sub(".*\/","",re.sub("\..*","",args.fasta))
        if args.verbose >= 3:
            print('parsed ',parsed['id'],end="")
            print()
        if args.verbose >= 3:
            print(parsed)
        list_of_seqs_per_code.append( parsed )
    else:
        if args.verbose >= 1:
            print("Could not parse that one")

the_datar = pd.DataFrame(
        [ i for i in 
            [k for k in list_of_seqs_per_code if k is not None] 
            ]
    )

vote_and_write_pdg(
    parsed['code'],
    the_datar,
    args.output)


# exec(open("msafasta2consensus.py").read())
