from argparse import ArgumentParser
import pandas as pd


def main():
    args = get_args()
    df = pd.read_csv(filepath_or_buffer=args.STRs,
                     sep="\t",
                     header=0)
    df["identifier"] = df["chrom"] + ":" + df["start"].astype(str) + "-" + df["repeatunit"]
    pivot = df.pivot(index='identifier',
                     columns='sample',
                     values=['locuscoverage', 'outlier', 'p_adj', 'bpInsertion', 'repeatUnits'])
    for feature in ['locuscoverage', 'outlier', 'bpInsertion', 'repeatUnits']:
        pivot.loc[:, (feature, "max")] = pivot[feature].max(axis="columns")
    pivot.loc[:, ('p_adj', "min")] = pivot[feature].min(axis="columns")

    pivot.columns = ['_'.join(col).strip() for col in pivot.columns.values]
    print(
        df[["identifier", "chrom", "start", "end", "repeatunit", "reflen"]]
        .drop_duplicates()
        .set_index('identifier')
        .join(pivot)
        .sort_values(["chrom", "start", "end"])
        .to_csv(sep="\t")
    )


def get_args():
    parser = ArgumentParser(description="Reformat STRetch STRs.tsv to stdout")
    parser.add_argument("STRs", help="File STRs.tsv generated by STRetch")
    return parser.parse_args()


if __name__ == '__main__':
    main()
