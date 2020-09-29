import os
import glob
import pandas
import json
import sys

config = json.load(open(sys.argv[1],"r"))
df = pandas.read_csv(config["samples"], header=0, index_col=0)

for idx, row in df.iterrows():
    print("########################################")
    exp = row["seq_dir"]
    fqp = row["fastq_pattern"]
    stype = row["sample_type"]
    base = config["seq_dir"][exp]
    matches = glob.glob(os.path.join(base, fqp))
    print(f"ID: {idx}")
    print(f"Sample type: {stype}")
    print(f"exp: {exp}\n")
    print(f"base: {base}\n")
    print(f"pattern: {fqp}\n")
    print(f"matches: {matches}")
    print(f"num matches: {len(matches)}")
    assert len(matches) > 0 and len(matches) < 3
    print("########################################\n\n")

