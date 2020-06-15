import csv
from Bio import SeqIO
import os
import collections

rule check_cog_db:
    input:
        query = config["query"],
        cog_seqs = config["cog_seqs"],
        metadata = config["cog_metadata"]
    params:
        field_to_match = config["search_field"]
    output:
        cog = os.path.join(config["outdir"],"query_in_cog.csv"),
        cog_seqs = os.path.join(config["outdir"],"query_in_cog.fasta"),
        not_cog = os.path.join(config["outdir"],"not_in_cog.csv")
    run:
        query_names = []
        in_cog_metadata = []
        in_cog_names = {}
        column_to_match = params.field_to_match
        with open(input.query,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                query_names.append(row["name"])

        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames
            for row in reader:
                for seq in query_names:

                    cog_id = row[column_to_match]
                    if seq == cog_id:
                        print(seq)
                        row["query_id"]=seq
                        row["cog_id"] = row[column_to_match]
                        row["query"]=row["sequence_name"]
                        row["closest"]=row["sequence_name"]
                        in_cog_metadata.append(row)
                        in_cog_names[column_to_match] = row["sequence_name"]

            print(f"Number of seqs found in metadata: {len(in_cog_metadata)}")
            with open(output.cog, "w") as fw:
                header_names.append("query_id")
                header_names.append("cog_id")
                header_names.append("query")
                header_names.append("closest")
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()
                writer.writerows(in_cog_metadata)

        with open(output.cog_seqs, "w") as fw:
            for record in SeqIO.parse(input.cog_seqs, "fasta"):
                for name in in_cog_names:
                    sequence_name = in_cog_names[name]
                    if sequence_name==record.id:
                        fw.write(f">{record.id} status=in_cog\n{record.seq}\n")
        
        with open(output.not_cog, "w") as fw:
            print("The following sequences were not found in the cog database:\n")

            for query in query_names:
                in_cog = False
                for name in in_cog_names:
                    if query == name:
                        in_cog = True
                if not in_cog:
                    fw.write(query + '\n')
                    print(f"{query}")
            print("If you wish to access sequences in the cog database\nwith your query, ensure you have the correct sequence id.")

rule check_cog_all:
    input:
        not_in_cog = rules.check_cog_db.output.not_cog,
        cog_seqs = config["all_cog_seqs"],
        all_cog_seqs = config["all_cog_metadata"]
    params:
        column_to_match = config["search_field"]
    output:
        cog = os.path.join(config["outdir"],"query_in_all_cog.csv"),
        cog_seqs = os.path.join(config["outdir"],"query_in_all_cog.fasta"),
        not_cog = os.path.join(config["outdir"],"not_in_all_cog.csv")
    run:
        not_cog = []
        in_all_cog_names = []
        in_all_cog_metadata = []
        with open(input.not_in_cog, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                not_cog.append(l)

        print("Checking for sequences in whole COG database")
        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames
            for row in reader:
                for seq in not_cog:

                    cog_id = row[column_to_match]
                    if seq == cog_id:
                        print(seq)
                        row["query_id"]=seq
                        row["cog_id"] = row[column_to_match]
                        row["query"]=row["sequence_name"]
                        row["closest"]=row["sequence_name"]
                        in_all_cog_metadata.append(row)
                        in_all_cog_names[column_to_match] = row["sequence_name"]

            print(f"The following sequences were found in COG-UK put hadn't passed the QC.\nLowering QC and adding them in to analysis now.")
            with open(output.cog, "w") as fw:
                header_names.append("query_id")
                header_names.append("cog_id")
                header_names.append("query")
                header_names.append("closest")
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()
                writer.writerows(in_all_cog_metadata)

        with open(output.cog_seqs, "w") as fw:
            
            for record in SeqIO.parse(input.cog_seqs, "fasta"):
                for name in in_all_cog_names:
                    sequence_name = in_all_cog_names[name]
                    if sequence_name==record.id:
                        fw.write(f">{record.id} status=in_all_cog\n{record.seq}\n")
        
        with open(output.not_cog, "w") as fw:
            c = 0
            print("The following sequences were not found in the cog database:\n")

            for query in not_cog:
                in_cog = False
                for name in in_all_cog_names:
                    if query == name:
                        in_cog = True
                if not in_cog:
                    c+=1
                    fw.write(query + '\n')
                    print(f"{query}")
        
            print(f"{c} sequences remaining not in COG, will find nearest COG sequence.")


rule get_closest_cog:
    input:
        snakefile = os.path.join(workflow.current_basedir,"find_closest_cog.smk"),
        reference_fasta = config["reference_fasta"],
        cog_seqs = config["cog_seqs"],
        cog_metadata = config["cog_metadata"],
        not_cog_csv = os.path.join(config["outdir"],"not_in_all_cog.csv")
    params:
        outdir= config["outdir"],
        # tempdir= config["tempdir"],
        path = workflow.current_basedir,
        cores = workflow.cores,
        force = config["force"],
        fasta = config["fasta"],
        query = config["post_qc_query"],
        quiet_mode = config["quiet_mode"],
        trim_start = config["trim_start"],
        trim_end = config["trim_end"]
    output:
        closest_cog = os.path.join(config["outdir"],"closest_cog.csv")
    run:
        if params.fasta != "":
            print(f"Passing {params.query} into processing pipeline.")
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force} "
                        "{params.quiet_mode} "
                        # "--directory {params.tempdir:q} "
                        "--config "
                        "outdir={params.outdir:q} "
                        # "tempdir={params.tempdir:q} "
                        "not_cog_csv={input.not_cog_csv:q} "
                        "post_qc_query={params.query:q} "
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

rule prune_out_catchments:
    input:
        tree = config["cog_tree"],
        metadata = rules.combine_metadata.output.combined_csv
    params:
        outdir = os.path.join(config["outdir"],"catchment_trees")
    output:
        txt = os.path.join(config["outdir"],"catchment_trees","catchment_tree_summary.txt")
    shell:
        """
        clusterfunk find_catchments -i {input.tree:q} \
        -o {params.outdir:q} \
        --metadata {input.metadata} \
        --index-column closest \
        --threshold 1 \
        --branch-count \
        --out-format newick \
        && touch {output.txt} 

        """

rule process_catchments:
    input:
        snakefile_collapse_after = os.path.join(workflow.current_basedir,"process_catchment_trees.smk"),
        snakefile_collapse_before = os.path.join(workflow.current_basedir,"process_collapsed_trees.smk"),
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        catchment_placeholder = os.path.join(config["outdir"],"catchment_trees","catchment_tree_summary.txt"),
        all_cog_seqs = config["all_cog_seqs"],
        in_all_cog_fasta = os.path.join(config["outdir"],"in_all_cog.fasta"),
        not_cog_csv = os.path.join(config["outdir"],"not_in_all_cog.csv")
    params:
        outdir= config["outdir"],
        # tempdir= config["tempdir"],
        path = workflow.current_basedir,
        query = config["post_qc_query"],
        cores = workflow.cores,
        delay_collapse = config["delay_collapse"],
        force = config["force"],
        fasta = config["fasta"],
        tree_dir = os.path.join(config["outdir"],"catchment_trees"),
        quiet_mode = config["quiet_mode"]
    output:
        tree_summary = os.path.join(config["outdir"],"combined_trees","collapse_report.txt")
    run:
        catchment_trees = []
        for r,d,f in os.walk(params.tree_dir):
            for fn in f:
                if fn.endswith(".tree"):
                    file_stem = ".".join(fn.split(".")[:-1])
                    catchment_trees.append(file_stem)
        catchment_str = ",".join(catchment_trees)

        num_in_all_cog = 0
        for record in SeqIO.parse(input.in_all_cog_fasta,"fasta"):
            num_in_all_cog +=1

        if params.fasta != "" or num_in_all_cog !=0:
            if params.delay_collapse==False:
                print(f"Passing {params.query} into processing pipeline.")
                shell("snakemake --nolock --snakefile {input.snakefile_collapse_before:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            # "--directory {params.tempdir:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={params.outdir:q} "
                            # "tempdir={params.tempdir:q} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "in_all_cog_fasta={input.in_all_cog_fasta:q} "
                            "post_qc_query={params.query:q} "
                            "all_cog_seqs={input.all_cog_seqs:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {params.cores}")
            else:
                print(f"Passing {params.query} into processing pipeline.")
                shell("snakemake --nolock --snakefile {input.snakefile_collapse_after:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            # "--directory {params.tempdir:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={params.outdir:q} "
                            # "tempdir={params.tempdir:q} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "in_all_cog_fasta={input.in_all_cog_fasta:q} "
                            "post_qc_query={params.query:q} "
                            "all_cog_seqs={input.all_cog_seqs:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {params.cores}")

        else:
            shell("touch {output.tree_summary:q}")

rule make_report:
    input:
        lineage_trees = os.path.join(config["outdir"],"combined_trees","collapse_report.txt"),
        query = config["query"],
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        full_cog_metadata = config["cog_metadata"],
        report_template = config["report_template"],
        font = config["font_file"],
        polytomy_figure = config["polytomy_figure"]
    params:
        treedir = os.path.join(config["outdir"],"restored_trees"),
        outdir = config["rel_outdir"],
        fields = config["fields"],
        figdir = os.path.join("./figures")
    output:
        poly_fig = os.path.join(config["outdir"],"figures","polytomies.png"),
        outfile = os.path.join(config["outdir"], "civet_report.md")
    shell:
        """
        cp {input.polytomy_figure:q} {output.poly_fig} &&
        make_report.py \
        --input-csv {input.query:q} \
        -f {params.fields:q} \
        --figdir {params.figdir:q} \
        --treedir {params.treedir:q} \
        --report-template {input.report_template:q} \
        --filtered-cog-metadata {input.combined_metadata:q} \
        --cog-metadata {input.full_cog_metadata:q} \
        --outfile {output.outfile:q} \
        --outdir {params.outdir:q} \
        --font-file {input.font:q}
        """
