"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  download-and-sketch.snakefile -n
"""

import os
import re
import csv
import pandas as pd

configfile: "prep_magsearch.v2.yml"
out_dir = config['outdir']
logs_dir = os.path.join(out_dir, "logs")

basename = config["basename"]
scaled = config["scaled"]
ksize = config["ksize"]
if not isinstance(scaled, list):
    config["scaled"] = [scaled]
if not isinstance(ksize, list):
    config["ksize"] = [ksize]

acc_info = pd.read_csv(config['accession_info'])
acc_info['species'] = acc_info['gtdb_taxonomy'].str.split('s__', 1, expand=True)[1]
acc_info['signame'] = acc_info["accession"] + ' ' + acc_info["species"]
acc_info['signame'].fillna(acc_info["accession"], inplace=True)
ACCS = acc_info['accession']

acc_info.set_index('accession', inplace=True)



rule all:
    #input: expand(os.path.join(out_dir, "data/genomic/{accession}.fna.gz"), accession = info_csv["accession"]),
    input: os.path.join(out_dir, f"{basename}.siglist.txt")

# download genbank genome details; make an info.csv file for entry.
rule make_genome_info_csv:
    output:
        csvfile = f'{out_dir}/info/{{acc}}.info.csv'
    shell: """
        python -Werror genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
    """

# download actual genomes!
rule download_matching_genome_wc:
     input:
         csvfile = ancient(f'{out_dir}/info/{{acc}}.info.csv')
     output:
         genome = os.path.join(out_dir, 'data', "{acc}_genomic.fna.gz")
     run:
         with open(input.csvfile, 'rt') as infp:
             r = csv.DictReader(infp)
             rows = list(r)
             assert len(rows) == 1
             row = rows[0]
             acc = row['acc']
             assert wildcards.acc.startswith(acc)
             url = row['genome_url']
             name = row['ncbi_tax_name']

             print(f"downloading genome for acc {acc}/{name} from NCBI...",
                   file=sys.stderr)
             with open(output.genome, 'wb') as outfp:
                 with urllib.request.urlopen(url) as response:
                     content = response.read()
                     outfp.write(content)
                     print(f"...wrote {len(content)} bytes to {output.genome}",
                           file=sys.stderr)


def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"

rule sourmash_sketch:
    input:
        #os.path.join(out_dir, "data", "{acc}.fna.gz"),
        ancient(os.path.join(out_dir, "data", "{acc}_genomic.fna.gz")),
    output:
        os.path.join(out_dir, "signatures", "{acc}.sig"),
    params:
        sketch_params=make_param_str(config["ksize"], config["scaled"]),
        signame = lambda w: acc_info.at[w.acc, "signame"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch", "{acc}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch", "{acc}.sketch.benchmark")
    conda: "envs/sourmash4.yml"
    group: "sketch"
    shell:
        """
        sourmash sketch dna -p {params.sketch_params} -o {output} --name {params.signame:q} {input} 2> {log}
        """


localrules: signames_to_file
rule signames_to_file:
    input: 
        sigs=expand(os.path.join(out_dir, "signatures", "{acc}.sig"), acc = ACCS),
    output: os.path.join(out_dir, f"{basename}.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                full_filename = os.path.abspath(str(inF))
                outF.write(full_filename + "\n")

