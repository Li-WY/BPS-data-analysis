#!/bin/bash

this_basename=$(basename "$1" .fastq )
# ${this_basename}

cat "$1" \
    | awk '{if (NR%4==1) {sub(/^/,">"); print} if (NR%4==2) {print} }' \
    | kalign --format fasta -o "./$2/${this_basename}.fasta" \
        > "./$2/${this_basename}.stdout" \
    || cat "$1" > "./$2/${this_basename}.fasta"

