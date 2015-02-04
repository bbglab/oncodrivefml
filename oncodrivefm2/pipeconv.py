import argparse
import gzip
import json
import math
import pandas

#"ENSG00000150594": [{"REF": "G", "SAMPLE": "853TD", "CHROMOSOME": "10", "TYPE": "subs", "ALT": "C", "POSITION": 112838418}]


def define_type(row):

    ref = row['REF'].replace("-", "")
    alt = row['ALT'].replace("-", "")
    is_inframe = math.fabs(len(ref) - len(alt)) % 3 == 0

    if is_inframe:
        return "subs"
    else:
        return "indel"


def convert(i, o):
    muts_df = pandas.DataFrame.from_csv(i, sep="\t", index_col=None)

    muts_df = muts_df[["START", "REF", "ALT", "SAMPLE", "CHR", "GENE"]]
    muts_df.rename( columns={"START": "POSITION", "CHR": "CHROMOSOME"}, inplace=True)
    muts_df.drop_duplicates(inplace=True)
    muts_df["TYPE"] = muts_df.apply(lambda row: define_type(row), axis=1)
    muts_df.to_json(o, orient='records')
    grouped = muts_df.groupby(["GENE"])
    data = {k: v.to_dict(orient='records') for k,v in grouped}
    with gzip.open(o, 'wt') as fd:
            json.dump(data, fd)

if __name__ == "__main__":
  # Parse the arguments
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", required=True)
    parser.add_argument("-o", required=True)


    args = parser.parse_args()

    convert(args.i, args.o)
