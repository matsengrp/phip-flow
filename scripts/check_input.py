import os
import glob
import pandas
import sys

# single argument in the sample metadata
df = pandas.read_csv(sys.argv[1])

for idx, row in df.iterrows():
    print("########################################")
    exp = row["seq_dir"]
    fqp = row["fastq_filename"]
    matches = glob.glob(os.path.join(exp, fqp))
    print(f"ID: {idx}")
    print(f"exp: {exp}\n")
    print(f"matches: {matches}")
    assert len(matches) == 1
    print("########################################\n\n")

