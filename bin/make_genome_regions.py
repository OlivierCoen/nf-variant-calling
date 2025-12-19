#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import argparse
import logging
from pathlib import Path

from Bio import SeqIO

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

OUTFILE = "genome_regions.tsv"


#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", type=Path, required=True, help="Genome file")
    parser.add_argument(
        "--chunksize", type=int, required=True, help="Size of the regions"
    )
    return parser.parse_args()


def make_genome_regions(genome_file: Path, chunksize: int):
    regions = []
    with open(genome_file, "r") as f:
        for chrom in SeqIO.parse(f, "fasta"):
            for i in range(0, len(chrom), chunksize):
                # make list chunk
                start = i
                end = min(i + chunksize, len(chrom))
                regions.append((chrom.id, start, end))
    return regions


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()
    regions = make_genome_regions(args.genome, args.chunksize)

    with open(OUTFILE, "w") as f:
        for region in regions:
            f.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")

    logger.info("Done")
