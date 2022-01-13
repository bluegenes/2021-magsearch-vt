"""
Minor cleanup and parsing of MAGsearch results
"""

import os
import sys
import pandas as pd
import argparse



def main(args):

    ms = pd.read_csv(args.magsearch_output, sep=",", quotechar="'", header=0, names=["search_genome", "metagenome", "containment"])
    # Fix names so it's easier to query
    ms['search_genome'] = ms['search_genome'].str.replace(r"'(?P<id>.*)'", lambda m: m.group("id"), regex=True)
    ms['metagenome'] = ms['metagenome'].str.replace(r".*/(?P<id>.*).sig.*", lambda m: m.group("id"), regex=True)
    ms['accession'] = ms['search_genome'].str.split(" ", 1, expand=True)[0]
    ms['sci_name'] = ms['search_genome'].str.split(" ", 1, expand=True)[1]

    num_search_genomes_with_matches = len(ms["search_genome"].unique())

    if not 0 <= args.f_containment_threshold <= 1:
        print(f"containment threshold is provided ({args.f_containment_threshold}), but is not between 0 (keep all results) and 1 (keep only results with 100% containment of at least one genome)\n")
        print("Exiting.\n")
        sys.exit(-1)
        # subset by containment
    ms = ms[ms['containment'] >= args.f_containment_threshold]
    perc_containment = args.f_containment_threshold * 100

    print(f"# Search genomes with metagenome matches: {num_search_genomes_with_matches}\n")
    num_search_genomes_with_matches_above_threshold = len(ms["search_genome"].unique())
    print(f"# Search genomes with metagenome matches above containment threshold ({perc_containment}%): {num_search_genomes_with_matches_above_threshold}\n")

    num_unique_metagenome_matches = len(ms["metagenome"].unique())
    print(f"# Metagenomes with > {perc_containment}% containment of at least one search genome: {num_unique_metagenome_matches}\n")

    # sort by containment
    ms = ms.sort_values(by="containment", ascending=False, ignore_index=True)
    # output results
    ms.to_csv(args.output_csv, index=False)

    # optionally output csv of "search genomes \t comma-separated matched metagenomes"
    if args.matched_metagenomes_tsv:
        mm = ms.groupby("search_genome")["metagenome"].apply(list).reset_index(name='matched_metagenomes')
        # convert matched genome lists to comma separated strings
        mm['matched_metagenomes'] = [','.join(map(str, l)) for l in mm['matched_metagenomes']]
        # make this tab separated bc matched metagenomes are comma sep
        mm.to_csv(args.matched_metagenomes_tsv, index=False, sep='\t')

    # optionally output list of metagenomes that were matched
    if args.metagenome_list:
        mgs = ms["metagenome"].drop_duplicates()
        mgs.to_csv(args.metagenome_list, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("magsearch_output")
    p.add_argument("--output-csv", required=True)
    p.add_argument("--f-containment-threshold", default=0, type=float, \
                   help="fraction containment should be between 0 (keep all results) and 1 (keep only results with 100%% containment of at least one genome)") # default save all results
    p.add_argument("--matched-metagenomes-tsv")
    p.add_argument("--metagenome-list")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
