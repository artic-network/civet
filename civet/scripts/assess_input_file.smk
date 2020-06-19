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
        cog = os.path.join(config["tempdir"],"query_in_cog.csv"),
        cog_seqs = os.path.join(config["tempdir"],"query_in_cog.fasta"),
        not_cog = os.path.join(config["tempdir"],"not_in_cog.csv")
    run:
        found = []
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
                        found.append(seq)
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
            print("\nThe following sequences were not found in the cog database:")

            for query in query_names:
                if query not in found:
                    fw.write(query + '\n')
                    print(f"{query}")
            print("If you wish to access sequences in the cog database\nwith your query, ensure you have the correct sequence id.")

rule check_cog_all:
    input:
        not_in_cog = rules.check_cog_db.output.not_cog,
        cog_seqs = config["all_cog_seqs"],
        metadata = config["all_cog_metadata"]
    params:
        field_to_match = config["search_field"]
    output:
        cog = os.path.join(config["tempdir"],"query_in_all_cog.csv"),
        cog_seqs = os.path.join(config["tempdir"],"query_in_all_cog.fasta"),
        not_cog = os.path.join(config["tempdir"],"not_in_all_cog.csv")
    run:
        found = []
        not_cog = []
        in_all_cog_names = {}
        in_all_cog_metadata = []
        with open(input.not_in_cog, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                not_cog.append(l)
        column_to_match = params.field_to_match
        print("Checking for sequences in whole COG database")
        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames
            for row in reader:
                for seq in not_cog:

                    cog_id = row[column_to_match]
                    if seq == cog_id:
                        print(seq)
                        found.append(seq)
                        row["query_id"]=seq
                        row["cog_id"] = row[column_to_match]
                        row["query"]=row["sequence_name"]
                        row["closest"]=row["sequence_name"]
                        in_all_cog_metadata.append(row)
                        in_all_cog_names[column_to_match] = row["sequence_name"]

            
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
        if found != []:
            print(f"The following sequences were found in COG-UK put hadn't passed the QC.\nLowering QC and adding them in to analysis now.")
            for i in found:
                print(f"    - {i}")

        with open(output.not_cog, "w") as fw:
            c = 0
            print("The following sequences were not found in the cog database:")
            for query in not_cog:
                if query not in found:
                    c+=1
                    fw.write(query + '\n')
                    print(f"{query}")

rule get_closest_cog:
    input:
        snakefile = os.path.join(workflow.current_basedir,"find_closest_cog.smk"),
        reference_fasta = config["reference_fasta"],
        cog_seqs = config["cog_seqs"],
        cog_metadata = config["cog_metadata"],
        seq_db = config["seq_db"],
        not_cog_csv = os.path.join(config["tempdir"],"not_in_all_cog.csv"),
        in_all_cog_metadata = os.path.join(config["tempdir"],"query_in_all_cog.csv"),
        in_all_cog_seqs = os.path.join(config["tempdir"],"query_in_all_cog.fasta")
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        path = workflow.current_basedir,
        cores = workflow.cores,
        force = config["force"],
        fasta = config["fasta"],
        search_field = config["search_field"],
        query = config["post_qc_query"],
        quiet_mode = config["quiet_mode"],
        stand_in_query = os.path.join(config["tempdir"], "temp.fasta"),
        trim_start = config["trim_start"],
        trim_end = config["trim_end"]
    output:
        closest_cog = os.path.join(config["tempdir"],"closest_cog.csv"),
        not_cog_query = os.path.join(config["tempdir"],"not_in_all_cog.fasta"),
        combined_query = os.path.join(config["tempdir"],"to_find_closest.fasta"),
        not_processed = os.path.join(config["tempdir"], "no_seq_to_process.csv")
    run:
        query_with_no_seq = []
        if params.fasta != "":
            not_cog = []
            with open(input.not_cog_csv, "r") as f:
                for l in f:
                    l = l.rstrip('\n')
                    not_cog.append(l)
            num_seqs = 0
            not_cog_with_seqs = []
            with open(output.not_cog_query, "w") as fw:
                for record in SeqIO.parse(params.query, "fasta"):
                    if record.id in not_cog:
                        num_seqs +=1
                        not_cog_with_seqs.append(record.id)
                        fw.write(f">{record.id} status=not_cog\n{record.seq}\n")
            for q in not_cog:
                if q not in not_cog_with_seqs:
                    query_with_no_seq.append(q)
            if num_seqs !=0:
                print(f"Passing {num_seqs} sequences into nearest COG search pipeline.")
                for i in not_cog_with_seqs:
                    print(f"    - {i}")
                shell("snakemake --nolock --snakefile {input.snakefile:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            "--directory {params.tempdir:q} "
                            "--config "
                            "outdir={params.outdir:q} "
                            "tempdir={params.tempdir:q} "
                            "seq_db={input.seq_db} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "post_qc_query={output.not_cog_query:q} "
                            "in_all_cog_metadata={input.in_all_cog_metadata} "
                            "in_all_cog_seqs={input.in_all_cog_seqs} "
                            "search_field={params.search_field} "
                            "cog_seqs={input.cog_seqs:q} "
                            "trim_start={params.trim_start} "
                            "trim_end={params.trim_end} "
                            "reference_fasta={input.reference_fasta:q} "
                            "cog_metadata={input.cog_metadata:q} "
                            "--cores {params.cores}")
        else:
            num_seqs = 0
            all_cog_seqs = []

            not_cog = []
            with open(input.not_cog_csv, "r") as f:
                for l in f:
                    l = l.rstrip('\n')
                    not_cog.append(l)

            for record in SeqIO.parse(input.in_all_cog_seqs,"fasta"):
                num_seqs +=1
                all_cog_seqs.append(record.id)

            for q in not_cog:
                if q not in all_cog_seqs:
                    query_with_no_seq.append(q)

            if num_seqs != 0:
                shell("touch {params.stand_in_query:q}")
                print(f"Passing {num_seqs} sequences into nearest COG search pipeline.")
                for i in all_cog_seqs:
                    print(f"    - {i}")
                shell("snakemake --nolock --snakefile {input.snakefile:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            "--directory {params.tempdir:q} "
                            "--config "
                            "outdir={params.outdir:q} "
                            "tempdir={params.tempdir:q} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "post_qc_query={params.stand_in_query:q} "
                            "in_all_cog_metadata={input.in_all_cog_metadata} "
                            "in_all_cog_seqs={input.in_all_cog_seqs} "
                            "search_field={params.search_field} "
                            "cog_seqs={input.cog_seqs:q} "
                            "trim_start={params.trim_start} "
                            "trim_end={params.trim_end} "
                            "reference_fasta={input.reference_fasta:q} "
                            "cog_metadata={input.cog_metadata:q} "
                            "--cores {params.cores}")
            else:
                with open(input.not_cog_csv, "r") as f:
                    for l in f:
                        l = l.rstrip('\n')
                        query_with_no_seq.append(l)
                    
                shell("touch {output.closest_cog:q} && touch {output.not_cog_query} && touch {output.combined_query}")
        with open(output.not_processed, "w") as fw:
            for q in list(set(query_with_no_seq)):
                fw.write(f"{q},fail=no sequence provided\n")

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
        outdir = os.path.join(config["tempdir"],"catchment_trees"),
        distance = config["distance"]
    output:
        txt = os.path.join(config["tempdir"],"catchment_trees","catchment_tree_summary.txt")
    shell:
        """
        clusterfunk find_catchments -i {input.tree:q} \
        -o {params.outdir:q} \
        --metadata {input.metadata} \
        --index-column closest \
        --threshold {params.distance} \
        --branch-count \
        --out-format newick \
        && touch {output.txt} 

        """

rule process_catchments:
    input:
        snakefile_collapse_after = os.path.join(workflow.current_basedir,"process_catchment_trees.smk"),
        snakefile_collapse_before = os.path.join(workflow.current_basedir,"process_collapsed_trees.smk"),
        snakefile_just_collapse = os.path.join(workflow.current_basedir,"just_collapse_trees.smk"),
        combined_metadata = rules.combine_metadata.output.combined_csv,
        catchment_placeholder = rules.prune_out_catchments.output.txt,
        all_cog_seqs = config["all_cog_seqs"],
        not_cog_query_seqs = rules.get_closest_cog.output.combined_query,
        not_cog_csv = rules.check_cog_all.output.not_cog
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        path = workflow.current_basedir,
        cores = workflow.cores,
        delay_collapse = config["delay_collapse"],
        force = config["force"],
        fasta = config["fasta"],
        tree_dir = os.path.join(config["tempdir"],"catchment_trees"),
        quiet_mode = config["quiet_mode"]
    output:
        tree_summary = os.path.join(config["outdir"],"local_trees","collapse_report.txt")
    run:
        catchment_trees = []
        for r,d,f in os.walk(params.tree_dir):
            for fn in f:
                if fn.endswith(".tree"):
                    file_stem = ".".join(fn.split(".")[:-1])
                    catchment_trees.append(file_stem)
        catchment_str = ",".join(catchment_trees)

        num_not_cog_query_seqs = 0
        for record in SeqIO.parse(input.not_cog_query_seqs,"fasta"):
            num_not_cog_query_seqs +=1

        if params.fasta != "" or num_not_cog_query_seqs !=0:
            if params.delay_collapse==False:
                print(f"Passing {input.not_cog_query_seqs} into processing pipeline.")
                shell("snakemake --nolock --snakefile {input.snakefile_collapse_before:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            "--directory {params.tempdir:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={params.outdir:q} "
                            "tempdir={params.tempdir:q} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "post_qc_query={input.not_cog_query_seqs:q} "
                            "all_cog_seqs={input.all_cog_seqs:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {params.cores}")
            else:
                print(f"Passing {input.not_cog_query_seqs} into processing pipeline.")
                shell("snakemake --nolock --snakefile {input.snakefile_collapse_after:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            "--directory {params.tempdir:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={params.outdir:q} "
                            "tempdir={params.tempdir:q} "
                            "not_cog_csv={input.not_cog_csv:q} "
                            "post_qc_query={input.not_cog_query_seqs:q} "
                            "all_cog_seqs={input.all_cog_seqs:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {params.cores}")

        else:
            print(f"No new sequences to add in, just collapsing trees.")
            shell("snakemake --nolock --snakefile {input.snakefile_just_collapse:q} "
                            "{params.force} "
                            "{params.quiet_mode} "
                            "--directory {params.tempdir:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={params.outdir:q} "
                            "tempdir={params.tempdir:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {params.cores}")

rule make_report:
    input:
        lineage_trees = rules.process_catchments.output.tree_summary,
        query = config["query"],
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        cog_global_metadata = config["cog_global_metadata"],
        report_template = config["report_template"],
        polytomy_figure = config["polytomy_figure"],
        no_seq = rules.get_closest_cog.output.not_processed
    params:
        treedir = os.path.join(config["outdir"],"local_trees"),
        outdir = config["rel_outdir"],
        fields = config["fields"],
        sc_source = config["sequencing_centre"],
        sc = config["sequencing_centre_file"],
        sc_flag = config["sequencing_centre_flag"],
        rel_figdir = os.path.join(".","figures"),
        figdir = os.path.join(config["outdir"],"figures"),
        failure = config["qc_fail"]
    output:
        poly_fig = os.path.join(config["outdir"],"figures","polytomies.png"),
        outfile = os.path.join(config["outdir"], "civet_report.md")
    run:
        if params.sc != "":
            shell("cp {params.sc_source} {params.sc}")
        shell(
        """
        cp {input.polytomy_figure:q} {output.poly_fig} 
        make_report.py \
        --input-csv {input.query:q} \
        -f {params.fields:q} \
        --figdir {params.rel_figdir:q} \
        {params.sc_flag} \
        {params.failure} \
        --no-seq-provided {input.no_seq} \
        --treedir {params.treedir:q} \
        --report-template {input.report_template:q} \
        --filtered-cog-metadata {input.combined_metadata:q} \
        --cog-metadata {input.cog_global_metadata:q} \
        --outfile {output.outfile:q} \
        --outdir {params.outdir:q} 
        """)
