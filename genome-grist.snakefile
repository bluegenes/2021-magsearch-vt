"""
Author: N. Tessa Pierce 
Run: snakemake -s genome-grist.snakefile --profile default
About: convenience snakefile for running genome-grist on magsearch results
"""

configfile: "config.grist-vt.yml"
out_dir = config['outdir']
logs_dir = os.path.join(out_dir, "logs")
ksize = config.get("ksize", ["21", "31", "51"])
if not isinstance(ksize, list):
    config["ksize"] = [ksize]

basename = config["basename"]
metagenome_list = config["metagenome_list"]
metagenomes = [x.rstrip() for x in open(metagenome_list, "r")]

# test one to start out
metagenomes = metagenomes[:1]
print(metagenomes)

rule all:
    input: 
        results=os.path.join(out_dir, f"{basename}.results.zip"),

localrules: sigs_to_zipfile
rule sigs_to_zipfile:
    input: config["query_genomes_siglist"]
    output: os.path.join(out_dir, "databases", "{basename}.zip")
    params:
        compression=config.get("gzip_compression", 9)
    log: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.sigs-to-zipfile.log")
    benchmark: os.path.join(logs_dir, "sigs-to-zipfile", "{basename}.sigs-to-zipfile.benchmark")
    conda: "envs/sourmash4.yml"
    shell:
        """
        python sigs-to-zipfile.py --compression {params.compression} --sig-pathlist {input} {output} 2> {log}
        """

rule write_grist_config:
    input: 
        database=os.path.join(out_dir, "databases", "{basename}.zip"),
        metagenomes=config["metagenome_list"],
    output: os.path.join(out_dir, "config.{basename}.grist")
    params:
        metagenome_trim_memory=config.get("metagenome_trim_memory", "1e9"),
        ksize=ksize
    run:
        with open(str(output), 'w') as out:
            out.write(f"outdir: {out_dir}\n")
            out.write(f"metagenome_trim_memory: {params.metagenome_trim_memory}\n")
            out.write(f"sourmash_database_glob_pattern: {input.database}\n") # can this have filepath, or does it need to be the basename only?
            out.write(f"sample:\n")
            for mg in metagenomes:
                out.write(f"  - {mg}\n")
            out.write(f"sourmash_database_ksize:\n")
            for k in params.ksize:
                out.write(f"  - {k}\n")
            out.write(f"sourmash_compute_ksizes:\n")
            for k in params.ksize:
                out.write(f"  - {k}\n")


rule run_genome_grist:
    input:  os.path.join(out_dir, f"config.{basename}.grist")
    output: expand(os.path.join(out_dir, "reports/report-{sample}.html"), sample=metagenomes)
    conda: "envs/genome-grist.yml"
    #threads: 16
    threads: 1
    resources: 
        #mem_mb=lambda wildcards, attempt: attempt*145000
        mem_mb=20000
    shell:
        """
        genome-grist run {input} --resources mem_mb=145000 -j {threads} summarize
        """

localrules: zip_genome_grist
rule zip_genome_grist:
    input:
        config=os.path.join(out_dir, "config.{basename}.grist"),
        results=rules.run_genome_grist.output
    output: os.path.join(out_dir, "{basename}.results.zip")
    params: 
        temp_zip = os.path.join(out_dir, "transfer.zip")
    conda: "envs/genome-grist.yml"
    shell:
        """
        genome-grist run {input.config} zip 
        mv {params.temp_zip} {output}
        """
