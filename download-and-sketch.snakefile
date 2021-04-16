"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  download-and-sketch.snakefile -n
"""

import os
import re
import pandas as pd

configfile: "prep_magsearch.yml"
out_dir = config['outdir']
logs_dir = os.path.join(out_dir, "logs")
taxon_info = config["taxa"]
basename = config["basename"]
scaled = config["scaled"]
ksize = config["ksize"]
if not isinstance(scaled, list):
    config["scaled"] = [scaled]
if not isinstance(ksize, list):
    config["ksize"] = [ksize]

# ctb checkpoint code to specify all the outputs
class Checkpoint_MakePattern:
    def __init__(self, pattern):
        self.pattern = pattern

    def get_names_from_csv(self):
        global genome_info
        genome_info = pd.read_csv(f'{out_dir}/data/{basename}.genome-info.csv')
        genome_info["signame"] = genome_info["accession"] + " " + genome_info["sci_name"] + " " + genome_info["strain"]
        accessions = genome_info['accession']
        genome_info.set_index("accession", inplace=True)
        return accessions

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'check_csv'; this will trigger an
        # exception until that rule has been run.
        checkpoints.check_csv.get(**w)

        accs = self.get_names_from_csv()

        pattern = expand(self.pattern, accession=accs, **w)
        return pattern


rule all:
    #input: expand(os.path.join(out_dir, "data/genomic/{accession}.fna.gz"), accession = info_csv["accession"]),
    input: os.path.join(out_dir, f"{basename}.siglist.txt")

rule download_ncbi_datasets_tool:
    output: "scripts/ncbi-datasets"
    shell:
        """
        # linux version
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
        chmod +x {output}
        """

localrules: find_accessions
rule find_accessions:
    output: os.path.join(out_dir, f"data/{basename}.genome-info.csv")
    params:
        taxa = " ".join(config["taxa"])
    log: os.path.join(logs_dir, "find_accessions", f"{basename}.find_accessions.log")
    conda: "envs/ncbi-datasets.yml"
    threads: 1
    shell:
        """
        python find-accessions-of-interest.py --taxa {params.taxa} --output-csv {output} > {log} 2>&1
        """

localrules: check_csv
checkpoint check_csv:
    input:
        accession_info=os.path.join(out_dir, f"data/{basename}.genome-info.csv")
    output: touch(f"{out_dir}/.check_csv")


rule ncbi_datasets_download:
    input: 
        check_csv=f"{out_dir}/.check_csv",
        tool="scripts/ncbi-datasets"
    output:
        genomic=protected(os.path.join(out_dir, "data/{accession}.fna.gz")),
        #protein=protected(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"))
    params:
        tmp= lambda w: os.path.join(out_dir, w.accession) + ".zip"
    threads: 1
    resources:
        mem_mb= 2000,
        runtime= 60
    log: os.path.join(logs_dir, "ncbi_datasets", "{accession}.download")
    shell:
        """
        echo "scripts/ncbi-datasets download genome accession {wildcards.accession} --exclude-rna --exclude-gff3 -f {params.tmp}" > {log} 
        scripts/ncbi-datasets download genome accession {wildcards.accession} --exclude-rna --exclude-gff3 -f {params.tmp}
        unzip -p {params.tmp} ncbi_dataset/data/{wildcards.accession}/*.fna | gzip -9 > {output.genomic}
        unzip {params.tmp} -d nd_{wildcards.accession}
        echo "fna filenames: " >> {log}
        ls -1 nd_{wildcards.accession}/ncbi_dataset/data/{wildcards.accession}/*.fna >> {log}
        rm -rf nd_{wildcards.accession} 
        rm -rf {params.tmp} 
        """
        #unzip -p {params.tmp} ncbi_dataset/data/{wildcards.accession}/*.faa | gzip -9 > {output.protein}

def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"

rule sourmash_sketch:
    input:
        ancient(os.path.join(out_dir, "data", "{accession}.fna.gz")),
    output:
        os.path.join(out_dir, "signatures", "{accession}.sig"),
    params:
        sketch_params=make_param_str(config["ksize"], config["scaled"]),
        signame = lambda w: genome_info.at[w.accession, "signame"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch", "{accession}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch", "{accession}.sketch.benchmark")
    conda: "envs/sourmash4.yml"
    group: "sketch"
    shell:
        """
        sourmash sketch dna -p {params.sketch_params} -o {output} --name {params.signame:q} {input} 2> {log}
        """


localrules: signames_to_file
rule signames_to_file:
    input: 
        csv=f"{out_dir}/.check_csv",
        sigs=Checkpoint_MakePattern(os.path.join(out_dir, "signatures", "{accession}.sig")),
    output: os.path.join(out_dir, f"{basename}.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                full_filename = os.path.abspath(str(inF))
                #outF.write(str(inF) + "\n")
                outF.write(full_filename + "\n")

