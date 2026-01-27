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
    parser.add_argument(
        "--ratio-overlap-to-chunksize",
        dest="ratio_overlap_to_chunksize",
        type=float,
        required=True,
        help="Ratio of overlap to chunksize",
    )
    return parser.parse_args()


def get_overlap_size(chunksize: int, ratio_overlap_to_chunksize: float) -> int:
    return max(1, int(chunksize * ratio_overlap_to_chunksize))


def make_genome_regions(
    genome_file: Path, chunksize: int, overlap: int
) -> list[tuple[str, int, int]]:
    regions = []
    with open(genome_file, "r") as f:
        for contig in SeqIO.parse(f, "fasta"):
            for i in range(0, len(contig), chunksize):
                # if lower than zero (for instance for the first chunk: sets to 0
                start = max(0, i - overlap)
                # if exceeds the right limit of the cotig, sets to the end of the contig
                end = min(len(contig), i + chunksize + overlap)
                regions.append((contig.id, start, end))
    return regions


#####################################################
#####################################################
# MAIN
#####################################################
#####################################################

if __name__ == "__main__":
    args = parse_args()

    overlap = get_overlap_size(args.chunksize, args.ratio_overlap_to_chunksize)
    logger.info(f"Overlap: {overlap}")

    regions = make_genome_regions(args.genome, args.chunksize, overlap)

    logger.info(f"Total nb of regions: {len(regions)}")

    with open(OUTFILE, "w") as f:
        for region in regions:
            f.write(f"{region[0]}\t{region[1]}\t{region[2]}\n")

    logger.info("Done")
