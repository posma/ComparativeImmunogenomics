#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
from math import ceil
from concurrent.futures import ProcessPoolExecutor


def generateUniformCoverageReadsFromRecords(chroms, label: str, readLen: int, coverage=30):
    os.makedirs("simulated_uniform", exist_ok=True)
    offset = ceil(readLen / coverage)
    out_path = f"simulated_uniform/{label}_{coverage}_{readLen}.fasta"
    with open(out_path, 'w') as ReadsFile:
        print(f"{out_path} created")
        readId = 0
        for chrom in chroms.keys():
            i = 0 
            while i + readLen < len(chroms[chrom]):
                print(f">read_{readId}\n{str(chroms[chrom][i: i + readLen].seq)}", file=ReadsFile)
                readId += 1
                i += offset

def processSubdir(args):
    subdirectory, readLen, coverage = args
    label = os.path.basename(subdirectory)
    if os.path.isdir(subdirectory):
        print(f"Processing: {label}")
        chroms = {}
        for file in os.listdir(subdirectory):
            if file.endswith(".fasta"):
                file_path = os.path.join(subdirectory, file)
                if os.path.isfile(file_path):
                    for record in SeqIO.parse(file_path, "fasta"):
                        chroms[record.id] = record
        generateUniformCoverageReadsFromRecords(chroms, label, readLen, coverage)

def main():
    parser = argparse.ArgumentParser(description="Simulate reads from genome files in subdirectories.")
    parser.add_argument("directory", help="Path to the parent directory with subdirectories")
    parser.add_argument("--coverage", type=int, default=10, help="Desired uniform coverage (default: 10)")
    parser.add_argument("--readlen", type=int, default=150, help="Read length (default: 150)")
    args = parser.parse_args()

    subdirs = [os.path.join(args.directory, d) for d in os.listdir(args.directory)
               if os.path.isdir(os.path.join(args.directory, d))]
    
    with ProcessPoolExecutor() as executor:
        executor.map(processSubdir, [(d, args.readlen, args.coverage) for d in subdirs])

if __name__ == "__main__":
    main()
