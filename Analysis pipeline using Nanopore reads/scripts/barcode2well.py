#!/bin/env/python3

import re
import rapidfuzz as rf
import sys
import argparse
import multiprocessing
import tqdm

parser = argparse.ArgumentParser()
parser.add_argument('--barcode-obs',default=sys.stdin,type=str)
parser.add_argument('--barcode-ref',type=str,default=None)
parser.add_argument('--cpus',type=int,default=1)
parser.add_argument('--report',type=str,default=None)
parser.add_argument('--max-distance',type=float,default=3)
parser.add_argument('--verbose','-v',action='count',default=0)
args = parser.parse_args()

try:
    with open(args.barcode_obs,'r') as f:
        barcode_obs = [ i.rstrip() for i in f.readlines()]
except:
    barcode_obs = [ i.rstrip() for i in args.barcode_obs.readlines()]

def parse_ref_codes(x):
    return { j[1]: j[0] for j in [i.strip().split("\t") for i in x] }

try:
    with open(args.barcode_ref,'r') as f:
        barcode_ref = parse_ref_codes(f.readlines())
except:
    barcode_ref = parse_ref_codes(args.barcode_ref.readlines())

matches = []
approx_matches = {}


def dist_check_return(z):
    distance = rf.distance.Levenshtein.distance(z[0],z[1],
        score_cutoff=z[2])
    if distance < z[2] + 1 :
        return [z[1],distance]
    else:
        return None

def try_matching(x,refs,max_dist):
    with multiprocessing.Pool(args.cpus) as mp:
        hits = []
        for dist in range(int(max_dist)):
            hits.extend( [ i for i in 
                    mp.map(dist_check_return,[ [x,j,dist+1] for j in refs]) 
                    if i is not None ] )
            if len(hits) > 0:
                break
    return hits
    
for i in tqdm.tqdm(barcode_obs):
    if i in barcode_ref.keys():
        matches.append([i,i,0])
    else:
        approx_matches[i] = try_matching(i,barcode_ref.keys(),args.max_distance)

for i in tqdm.tqdm(approx_matches):
    if len(approx_matches[i]) == 1:
        matches.append([i,approx_matches[i][0][0],approx_matches[i][0][1]])
    elif len(approx_matches[i]) == 0:
        matches.append([i,i,None])
    elif len(approx_matches[i]) > 1:
        matches.append(
            [i,
                "multi_"+"_".join([ x[0] for x in approx_matches[i]]),
                "_".join([ str(x[1]) for x in approx_matches[i]])
                ]
            )

with sys.stdout as f:
    for i in matches:
        if i[1] in barcode_ref.keys():
            f.write("\t".join(i[0:2])+"\t"+barcode_ref[i[1]]+
                "\t"+str(i[2])+"\n")
        else:
            f.write("\t".join(i[0:2])+"\t"+"unknown"+
                "\t"+str(i[2])+"\n")

if args.report is not None:
    with open(args.report,'w') as f:
        for i in approx_matches:
            f.write(i,approx_matches[i][0],approx_matches[i][1])

