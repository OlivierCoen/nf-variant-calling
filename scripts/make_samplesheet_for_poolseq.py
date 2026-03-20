#!/usr/bin/env python

from pathlib import Path

import pandas as pd

FOLDER_PATTERN = "X204SC24107022-Z01-F001"
FILE_PREFIXED_TO_INCLUDE = ["ST24_M", "ST24_F", "LU26_F", "LU26_M", "PO21_F", "PO21_M"]
OUTFILE = "samplesheet.poolseq.csv"

folders = Path().cwd().glob(f"{FOLDER_PATTERN}_*")

fastq_files = []
for folder in folders:
    data_subfolder = folder / "01.RawData"
    for sample_folder in data_subfolder.iterdir():
        for file in Path(sample_folder).iterdir():
            if file.suffix == ".gz":
                for prefix in FILE_PREFIXED_TO_INCLUDE:
                    if file.name.startswith(prefix):
                        fastq_files.append(file)
                        break

metadatas = []
for file in fastq_files:
    file_metadatas = file.stem.replace(".fq", "").split("_")

    population = file_metadatas[0]
    sex = file_metadatas[1]
    chip = file_metadatas[3]
    lane = file_metadatas[4]
    read_type = file_metadatas[5]

    sample = f"{population}_{sex}"

    # special case for some LU26_F files
    if lane.startswith("S"):
        file_id = file_metadatas[6]
        if file_id[0] == "I":
            lane = 2
        else:
            lane = 3
        read_type = file_id[1]

    complete_lane = f"{chip}_{lane}"

    file_metadata = dict(
        sample=sample,
        lane=complete_lane,
        read_type=read_type,
        file=str(file.absolute()),
    )
    metadatas.append(file_metadata)

df = pd.DataFrame(metadatas)

grouped_df = df.groupby(["sample", "lane"]).agg(list)


def separate_fastq_files(read_types: list[int], files: list[str]):
    if read_types[0] == "1":
        return files[0], files[1]
    else:
        return files[1], files[0]


final_df = (
    grouped_df.apply(
        lambda x: separate_fastq_files(x["read_type"], x["file"]),
        axis=1,
        result_type="expand",
    )
    .reset_index()
    .rename(columns={0: "fastq_1", 1: "fastq_2"})
)

print(final_df)

print("\nWriting samplesheet file...")
final_df.to_csv(OUTFILE, index=False, header=True)
