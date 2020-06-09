import csv
from Bio import SeqIO
import os
import collections

rule check_cog_db:
    input:
        query = config["query"],
        cog_seqs = config["cog_seqs"],
        metadata = config["cog_metadata"]
    output:
        cog = os.path.join(config["outdir"],"query_in_cog.csv"),
        cog_seqs = os.path.join(config["outdir"],"query_in_cog.fasta"),
        not_cog = os.path.join(config["outdir"],"not_in_cog.csv")
    run:
        query_names = []
        in_cog_metadata = []
        in_cog_names = set()
        with open(input.query,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                query_names.append(row["name"])

        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames
            for row in reader:
                for seq in query_names:
                    seq_name = row["sequence_name"].split("/")
                    if seq in seq_name:
                        print(seq)
                        row["query"]=row["sequence_name"]
                        row["closest"]=row["sequence_name"]
                        in_cog_metadata.append(row)
                        in_cog_names.add(row["sequence_name"])

            print(f"Number of seqs found in metadata: {len(in_cog_metadata)}")
            with open(output.cog, "w") as fw:
                header_names.append("query")
                header_names.append("closest")
                writer = csv.DictWriter(fw, fieldnames=header_names)
                writer.writeheader()
                writer.writerows(in_cog_metadata)

        fw = open(output.cog_seqs, "w")
        for record in SeqIO.parse(input.cog_seqs, "fasta"):
            for name in query_names:
                seq_name = record.id.split("/")
                if name in seq_name:
                    fw.write(f">{record.id}\n{record.seq}\n")
        
        with open(output.not_cog, "w") as fw:
            print("The following sequences were not found in the cog database:\n")

            for query in query_names:
                in_cog = False
                for name in in_cog_names:
                    if query in name:
                        in_cog = True
                if not in_cog:
                    fw.write(query + '\n')
                    print(f"{query}")
            print("If you wish to access sequences in the cog database\nwith your query, ensure you have the correct sequence id.")

rule get_closest_cog:
    input:
        snakefile = os.path.join(workflow.current_basedir,"find_closest_cog.smk"),
        reference_fasta = config["reference_fasta"],
        cog_seqs = config["cog_seqs"],
        cog_metadata = config["cog_metadata"],
        query = config["post_qc_query"],
        not_cog_csv = os.path.join(config["outdir"],"not_in_cog.csv")
    params:
        outdir= config["outdir"],
        # tempdir= config["tempdir"],
        path = workflow.current_basedir,
        cores = workflow.cores,
        force = config["force"],
        fasta = config["fasta"],
        trim_start = config["trim_start"],
        trim_end = config["trim_end"]
    output:
        closest_cog = os.path.join(config["outdir"],"closest_cog.csv")
    run:
        if params.fasta != "":
            print(f"Passing {input.query} into processing pipeline.")
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force} "
                        # "--directory {params.tempdir:q} "
                        "--config "
                        "outdir={params.outdir:q} "
                        # "tempdir={params.tempdir:q} "
                        "not_cog_csv={input.not_cog_csv:q} "
                        "post_qc_query={input.query:q} "
                        "cog_seqs={input.cog_seqs:q} "
                        "trim_start={params.trim_start} "
                        "trim_end={params.trim_end} "
                        "reference_fasta={input.reference_fasta:q} "
                        "cog_metadata={input.cog_metadata:q} "
                        "--cores {params.cores}")
        else:
            shell("touch {output.closest_cog:q}")

rule combine_metadata:
    input:
        closest_cog = rules.get_closest_cog.output.closest_cog,
        in_cog = rules.check_cog_db.output.cog
    output:
        combined_csv = os.path.join(config["outdir"],"combined_metadata.csv")
    run:
        with open(output.combined_csv,"w") as fw:
            with open(input.in_cog, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    fw.write(l + '\n')
            with open(input.closest_cog, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    if "sequence_name" in l:
                        pass
                    else:
                        fw.write(l + '\n')

rule get_remote_lineage_trees:
    input:
        combined_csv = os.path.join(config["outdir"],"combined_metadata.csv")
    params:
        outdir = os.path.join(config["outdir"],"lineage_trees"),
        uun = config["username"]
    output:
        lineage_trees_remote = os.path.join(config["outdir"],"lineage_trees","lineage_tree_summary.remote.txt")
    run:
        lineages = set()
        with open(input.combined_csv, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                lineages.add(row["uk_lineage"])
        fw = open(output.lineage_trees_remote, "w")
        for lineage in lineages:
            try:
                shell(f"""scp {params.uun}@bham.covid19.climb.ac.uk:\
/cephfs/covid/bham/artifacts/published/latest/phylogenetics/trees/uk_lineages/uk_lineage_{lineage}.tree {params.outdir}"""
                )
                fw.write(lineage+ ',pass\n')
            except:
                fw.write(lineage+ ',fail\n')
        fw.close()


rule get_lineage_trees:
    input:
        combined_csv = os.path.join(config["outdir"],"combined_metadata.csv")
    params:
        outdir = os.path.join(config["outdir"],"lineage_trees")
    output:
        lineage_trees = os.path.join(config["outdir"],"lineage_trees","lineage_tree_summary.txt")
    run:
        lineages = set()
        with open(input.combined_csv, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                lineages.add(row["uk_lineage"])
        fw = open(output.lineage_trees, "w")
        for lineage in lineages:
            try:
                shell(f"""scp /cephfs/covid/bham/artifacts/published/latest/phylogenetics/trees/uk_lineages/uk_lineage_{lineage}.tree \
                    {params.outdir}""")
                fw.write(lineage + 'pass\n')
            except:
                fw.write(lineage+ ',fail\n')
        fw.close()

rule make_report:
    input:
        lineage_trees = os.path.join(config["outdir"],"lineage_trees","lineage_tree_summary.txt"),
        query = config["query"],
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        full_cog_metadata = config["cog_metadata"],
        report_template = config["report_template"],
        font = config["font_file"] 
    params:
        tree_dir = os.path.join(config["outdir"],"lineage_trees"),
        outdir = config["outdir"],
        fields = config["fields"]
    output:
        outfile = os.path.join(config["outdir"], "civet_report.pmd")
    shell:
        """
        make_report.py \
        --input-csv {input.query:q} \
        -f {params.fields:q} \
        -t {params.tree_dir:q} \
        --report-template {input.report_template} \
        --filtered-cog-metadata {input.combined_metadata:q} \
        --cog-metadata {input.full_cog_metadata:q} \
        --outfile {output.outfile:q} \
        --outdir {params.outdir:q} \
        --font-file {input.font}
        """


rule remote_report:
    input:
        lineage_trees = os.path.join(config["outdir"],"lineage_trees","lineage_tree_summary.remote.txt"),
        query = config["query"],
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        full_cog_metadata = config["cog_metadata"],
        report_template = config["report_template"],
        font = config["font_file"] 
    params:
        tree_dir = os.path.join(config["outdir"],"lineage_trees"),
        outdir = config["outdir"],
        fields = config["fields"]
    output:
        outfile = os.path.join(config["outdir"], "civet_report.remote.md")
    shell:
        """
        make_report.py \
        --input-csv {input.query:q} \
        -f {params.fields:q} \
        -t {params.tree_dir:q} \
        --report-template {input.report_template} \
        --filtered-cog-metadata {input.combined_metadata:q} \
        --cog-metadata {input.full_cog_metadata:q} \
        --outfile {output.outfile:q} \
        --outdir {params.outdir:q} \
        --font-file {input.font}
        """

