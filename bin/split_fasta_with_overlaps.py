#!/usr/bin/env python

import argparse
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

RATIO_OVERLAP_TO_CHUNKSIZE = 0.1


def parse_args():
    parser = argparse.ArgumentParser(
        description="Split Fasta file into chunks of equal size with overlaps"
    )
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--chunksize", type=int, required=True)
    parser.add_argument("--outdir", dest="outdir", type=str, required=True)
    return parser.parse_args()


def get_overlap_size(chunksize: int) -> int:
    return max(1, int(chunksize * RATIO_OVERLAP_TO_CHUNKSIZE))


def get_contig_chunks(
    contig: SeqRecord, chunksize: int, overlap: int
) -> list[SeqRecord]:
    contig_chunks = []
    for i in range(0, len(contig), chunksize):
        # if lower than zero (for instance for the first chunk: sets to 0
        chunk_start = max(0, i - overlap)
        # if exceeds the right limit of the cotig, sets to the end of the contig
        chunk_end = min(len(contig), i + chunksize + overlap)
        # create a SeqRecord for the chunk
        chunk_sequence = str(contig.seq)[chunk_start:chunk_end]
        chunk_id = f"{contig.id}_chunk_{i // chunksize}"
        chunk_record = SeqRecord(Seq(chunk_sequence), id=chunk_id)
        contig_chunks.append(chunk_record)
    return contig_chunks


def main():
    args = parse_args()

    logger.info(f"Computing chunks from {args.fasta} with chunksize {args.chunksize}")

    overlap = get_overlap_size(args.chunksize)
    logger.info(f"Overlap: {overlap}")

    contig_chunks = []
    for contig in SeqIO.parse(args.fasta, "fasta"):
        chunks = get_contig_chunks(contig, args.chunksize, overlap)
        contig_chunks += chunks

    logger.info(f"Total chunks: {len(contig_chunks)}")

    # Create output directory if it doesn't exist
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    for chunk in contig_chunks:
        chunk_file = Path(args.outdir) / f"{chunk.id}.fasta"
        with open(chunk_file, "w") as fout:
            SeqIO.write([chunk], fout, "fasta")

    logger.info("Done")


if __name__ == "__main__":
    main()
