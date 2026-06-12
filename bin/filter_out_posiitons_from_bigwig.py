from pathlib import Path
import argparse
import polars as pl
import pybigtools

#####################################################
#####################################################
# FUNCTIONS
#####################################################
#####################################################


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate and filter VCF files")
    parser.add_argument(
        "--bw1",
        type=Path,
        dest="bw_file1",
        required=True,
        help="Path to bigwig file",
    )
    parser.add_argument(
        "--bw2",
        type=Path,
        dest="bw_file2",
        required=True,
        help="Path to bigwig file containing positions to exclude",
    )
    parser.add_argument(
        "--out",
        type=Path,
        dest="output_file",
        required=True,
        help="Path to output file",
    )
    return parser.parse_args()

def parse_vcf_header(vcf_file: Path) -> list:
    header_lines = []
    with open(vcf_file, "rb") as fin:
        for line in fin:
            line = line.decode("utf-8")
            if line.startswith("##"):
                header_lines.append(line)
            else:
                break
    return header_lines

def parse_vcf_data(vcf_file: Path) -> pl.DataFrame:
    return pl.read_csv(
        vcf_file, separator="\t", has_header=True, comment_prefix="##", low_memory=True
    )
    
#####################################################
#####################################################
# MAIN
#####################################################
#####################################################


def main():

    args = parse_args()
    
    bw1 = pybigtools.open(args.bw_file1, "r")
    bw2 = pybigtools.open(args.bw_file2, "r")
    chrom2sizes = bw1.chroms()

    tups = []
    for chrom in chrom2sizes:
        bw2_records = list(bw2.records(chrom))
        for rec1 in bw1.records(chrom):
            found_match = False
            for rec2 in bw2_records:
                if rec1[0] == rec2[0] and rec1[1] == rec2[1]:
                    print("Same positions:", rec1, rec2)
                    found_match = True
                    break
            if not found_match:
                tups.append((chrom, rec1[0], rec1[1], rec1[2]))

    out_bw = pybigtools.open(args.output_file, "w")
    out_bw.write(chrom2sizes, iter(tups))
    out_bw.close()

if __name__ == "__main__":
    main()
