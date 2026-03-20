#!/usr/bin/env python

from pathlib import Path

import pandas as pd

FOLDER_PATTERN = "X204SC24107022-Z01-F001"
FILE_PREFIXED_TO_INCLUDE = ["CO25_f10", "CO25_f15", "N39", "N40", "N47", "PO8_f02", "PO15_f01", "PO22_f03", "PO22_f04", "PO22_f08", "PO22_f11", "PO22_f12"]
OUTFILE = "samplesheet.individuals.csv"

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

    if len(file_metadatas) > 5:
        population = file_metadatas[0]
        sub_population = file_metadatas[1]
        chip = file_metadatas[3]
        lane = file_metadatas[4]
        read_type = file_metadatas[5]

        sample = f"{population}_{sub_population}"

    else:  # special case for N47, N40 and N39 samples
        sample = file_metadatas[0]
        chip = file_metadatas[2]
        lane = file_metadatas[3]
        read_type = file_metadatas[4]

        population = sample[0]
        sub_population = sample[1:]

    complete_lane = f"{chip}_{lane}"

    file_metadata = dict(
        sample=sample,
        population=population,
        sub_population=sub_population,
        lane=complete_lane,
        read_type=read_type,
        file=str(file.absolute()),
    )
    metadatas.append(file_metadata)

df = pd.DataFrame(metadatas)

grouped_df = df.groupby(["sample", "lane", "population", "sub_population"]).agg(list)


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
