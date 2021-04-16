"""
First, find genome assemblies of interest using ncbi.datasets api
See example at https://github.com/ncbi/datasets/blob/master/examples/jupyter/ncbi-datasets-pylib/ncbi-datasets-assembly.ipynb
"""

import os
import sys
import pandas as pd
import argparse
from collections import namedtuple

import ncbi.datasets

AssemblyInfo = namedtuple('AssemblyInfo',
                          'accession, display_name, tax_id, parent_tax_id, sci_name, strain, title, submission_date, seq_length, assembly_level')

# ncbi datasets api instance
api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())

def find_genomes(taxon, acc_set):
    taxon_info = []
    genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=taxon,limit='all')
    print(f"Number of assemblies in the group '{taxon}': {genome_summary.total_count}")

    for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
        if not assembly.annotation_metadata:
            continue
        unique_acc = assembly.assembly_accession[3:]
        # want to ignore GCA, GCF. Refseq and genbank should have same assembly, right?!
        if unique_acc in acc_set:
            print(f"\t{unique_acc} has already been added! skipping...")
            continue
        print(f"\tfinding info for {assembly.assembly_accession}")
        #for fileinfo in assembly.annotation_metadata.file:
        #    if fileinfo.type == "PROT_FASTA":
         #       print("proteins found!")
        acc_set.add(unique_acc)
        taxon_info.append(AssemblyInfo(accession=assembly.assembly_accession, \
                                     display_name=assembly.display_name, \
                                     tax_id=assembly.org.tax_id, \
                                     parent_tax_id=assembly.org.parent_tax_id, \
                                     sci_name = assembly.org.sci_name, \
                                     strain = assembly.org.strain, \
                                     title = assembly.org.title, \
                                     submission_date = assembly.submission_date, \
                                     seq_length = assembly.seq_length, \
                                     assembly_level = assembly.assembly_level))

    print(f"Number of unique assemblies for {taxon}: {len(taxon_info)}")
    return taxon_info, acc_set


def main(args):

    taxon_info = args.taxa
    taxon_info = [x.replace("_", " ") for x in taxon_info] # replace underscores with spaces
    all_genome_info=[]
    acc_set=set()
    for taxon in taxon_info:
        tax_info, acc_set = find_genomes(taxon, acc_set)
        all_genome_info += tax_info

    # convert namedtuple to pandas dataframe
    infoDF = pd.DataFrame.from_records(all_genome_info, columns = AssemblyInfo._fields)
    infoDF.to_csv(args.output_csv, index=False)

    print(f"Total of unique assemblies found {len(infoDF.index)}")

    if args.output_complete_genome:
        # select only those with "Complete Genome":
        completeG = infoDF[infoDF["assembly_level"] == "Complete Genome"]
        print(f"Number of unique Complete Genome assemblies {len(completeG.index)}")
        completeG.to_csv(args.output_complete_genome, index=False)

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--taxa", nargs="+", help = "one or more ncbi taxid or scientific names, separated by a space (for scientific names, use underscores instead of spaces within each name)")
    p.add_argument("--output-csv", default = "genome_info.csv")
    p.add_argument("--output-complete-genome") # optionally also output just the accessions with a complete genome
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
