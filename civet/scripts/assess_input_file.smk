import csv
from Bio import SeqIO
import os
import collections

#print(config)
rule check_cog_db:
    input:
        query = config["query"],
        cog_seqs = config["cog_seqs"],
        cog_metadata = config["cog_metadata"]
    params:
        field_to_match = config["search_field"]
    output:
        cog = os.path.join(config["tempdir"],"query_in_cog.csv"),
        cog_seqs = os.path.join(config["tempdir"],"query_in_cog.fasta"),
        not_cog = os.path.join(config["tempdir"],"not_in_cog.csv")
    shell:
        """
        check_cog_db.py --query {input.query:q} \
                        --cog-seqs {input.cog_seqs:q} \
                        --cog-metadata {input.cog_metadata:q} \
                        --field {params.field_to_match} \
                        --in-metadata {output.cog:q} \
                        --in-seqs {output.cog_seqs:q} \
                        --not-in-cog {output.not_cog:q}
        """
        
rule check_cog_all:
    input:
        not_in_cog = rules.check_cog_db.output.not_cog,
        cog_seqs = config["all_cog_seqs"],
        cog_metadata = config["all_cog_metadata"]
    params:
        field_to_match = config["search_field"]
    output:
        cog = os.path.join(config["tempdir"],"query_in_all_cog.csv"),
        cog_seqs = os.path.join(config["tempdir"],"query_in_all_cog.fasta"),
        not_cog = os.path.join(config["tempdir"],"not_in_all_cog.csv")
    shell:
        """
        check_cog_db.py --query {input.not_in_cog:q} \
                        --cog-seqs {input.cog_seqs:q} \
                        --cog-metadata {input.cog_metadata:q} \
                        --field {params.field_to_match} \
                        --in-metadata {output.cog:q} \
                        --in-seqs {output.cog_seqs:q} \
                        --not-in-cog {output.not_cog:q}
        """

rule get_closest_cog:
    input:
        snakefile = os.path.join(workflow.current_basedir,"find_closest_cog.smk"),
        reference_fasta = config["reference_fasta"],
        cog_seqs = config["cog_seqs"],
        cog_metadata = config["cog_metadata"],
        seq_db = config["seq_db"],
        not_cog_csv = rules.check_cog_all.output.not_cog, #use
        in_all_cog_metadata = rules.check_cog_all.output.cog,
        in_all_cog_seqs = rules.check_cog_all.output.cog_seqs #use 
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        path = workflow.current_basedir,
        cores = workflow.cores,
        force = config["force"],
        fasta = config["fasta"], #use
        search_field = config["search_field"],
        query = config["post_qc_query"], #use
        stand_in_query = os.path.join(config["tempdir"], "temp.fasta"),
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        quiet_mode = config["quiet_mode"]
    output:
        closest_cog = os.path.join(config["tempdir"],"closest_cog.csv"),
        combined_query = os.path.join(config["tempdir"],"to_find_closest.fasta"),
        aligned_query = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta"),
        not_processed = os.path.join(config["tempdir"], "no_seq_to_process.csv")
    run:
        query_with_no_seq = []
        to_find_closest = {}

        for record in SeqIO.parse(input.in_all_cog_seqs,"fasta"):
            to_find_closest[record.id] = ("in_all_cog",record.seq)

        not_cog = []
        with open(input.not_cog_csv, newline = "") as f: # getting list of non-cog queries
            reader = csv.DictReader(f)
            for row in reader:
                not_cog.append(row["name"])

        if params.fasta != "":
             # get set with supplied sequences
                print("Not in COG but have a sequence supplied:")
                for record in SeqIO.parse(params.query, "fasta"):
                    if record.id in not_cog:
                        to_find_closest[record.id] = ("not_cog",record.seq) # overwrites with supplied seq if found in all cog

        with open(output.combined_query, "w") as fw:
            for seq in to_find_closest:
                fw.write(f">{seq} status={to_find_closest[seq][0]}\n{to_find_closest[seq][1]}\n")

        for query in not_cog: # get set with no sequences supplied
            if query not in to_find_closest:
                query_with_no_seq.append(query)

        if to_find_closest != {}:
            print(f"Passing {len(to_find_closest)} sequences into nearest COG search pipeline:")
            for seq in to_find_closest:
                print(f"    - {seq}    {to_find_closest[seq][0]}")
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force} "
                        "{params.quiet_mode} "
                        "--directory {params.tempdir:q} "
                        "--config "
                        "tempdir={params.tempdir:q} "
                        "seq_db={input.seq_db:q} "
                        "to_find_closest={output.combined_query} "
                        "search_field={params.search_field} "
                        "trim_start={params.trim_start} "
                        "trim_end={params.trim_end} "
                        "reference_fasta={input.reference_fasta:q} "
                        "cog_metadata={input.cog_metadata:q} "
                        "--cores {params.cores}")

        else:
            shell("touch {output.closest_cog:q} && touch {output.aligned_query} ")

        with open(output.not_processed, "w") as fw:
            for query in list(set(query_with_no_seq)):
                fw.write(f"{query},fail=no sequence provided\n")

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
        --metadata {input.metadata:q} \
        --index-column closest \
        --threshold {params.distance} \
        --branch-count \
        --out-format newick \
        && touch {output.txt:q} 

        """

rule process_catchments:
    input:
        snakefile_collapse_after = os.path.join(workflow.current_basedir,"process_catchment_trees.smk"), #alternative snakefiles
        snakefile_collapse_before = os.path.join(workflow.current_basedir,"process_collapsed_trees.smk"),
        snakefile_just_collapse = os.path.join(workflow.current_basedir,"just_collapse_trees.smk"),
        combined_metadata = rules.combine_metadata.output.combined_csv, 
        query_seqs = rules.get_closest_cog.output.aligned_query, #datafunk-processed seqs
        catchment_prompt = rules.prune_out_catchments.output.txt,
        all_cog_seqs = config["all_cog_seqs"],
        outgroup_fasta = config["outgroup_fasta"],
        cog_global_seqs = config["cog_global_seqs"]
        # not_cog_csv = rules.check_cog_all.output.not_cog
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        path = workflow.current_basedir,
        threshold = config["threshold"],
        delay_collapse = config["delay_collapse"],
        
        fasta = config["fasta"],
        tree_dir = os.path.join(config["tempdir"],"catchment_trees"),

        cores = workflow.cores,
        force = config["force"],
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
        catchment_str = ",".join(catchment_trees) #to pass to snakemake pipeline

        query_seqs = 0
        for record in SeqIO.parse(input.query_seqs,"fasta"):
            query_seqs +=1

        if query_seqs !=0:
            if params.delay_collapse==False:
                snakefile = input.snakefile_collapse_before
            else:
                snakefile = input.snakefile_collapse_after

            snakestring = f"'{snakefile}' "
            print(f"Passing {input.query_seqs} into processing pipeline.")
            shell(f"snakemake --nolock --snakefile {snakestring}"
                        "{params.force} "
                        "{params.quiet_mode} "
                        "--directory {params.tempdir:q} "
                        "--config "
                        f"catchment_str={catchment_str} "
                        "outdir={params.outdir:q} "
                        "tempdir={params.tempdir:q} "
                        "outgroup_fasta={input.outgroup_fasta:q} "
                        "aligned_query_seqs={input.query_seqs:q} "
                        "all_cog_seqs={input.all_cog_seqs:q} "
                        "cog_global_seqs={input.cog_global_seqs:q} "
                        "combined_metadata={input.combined_metadata:q} "
                        "threshold={params.threshold} "
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

rule regional_mapping:
    input:
        query = config['query'],
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        cog_global_metadata = config["cog_global_metadata"]
    params:
        local_lineages = config["local_lineages"],
        daterestrict = config["date_restriction"],
        datestart = config["date_range_start"],
        dateend = config["date_range_end"],
        datewindow = config["date_window"],
        outdir = config["rel_outdir"],
        tempdir = config['tempdir'],
        figdir = os.path.join(config["outdir"],'figures', "regional_mapping")
    output:
        central = os.path.join(config["tempdir"], "central_map_ukLin.vl.json"),
        neighboring = os.path.join(config["tempdir"], "neighboring_map_ukLin.vl.json"),
        region = os.path.join(config["tempdir"], "region_map_ukLin.vl.json")
    run:
        if params.local_lineages:
            shell("""
            mkdir {params.figdir} -p
            """)
            shell(
            """
        local_scale_analysis.py \
        --uk-map {params.mapfile:q} \
        --hb-translation {params.hbtrans:q} \
        --date-restriction {params.daterestrict:q} \
        --date-pair-start {params.datestart:q} \
        --date-pair-end {params.dateend:q} \
        --date-window {params.datewindow:q} \
        --cog-meta-global {input.cog_global_metadata:q} \
        --user-sample-data {input.query:q} \
        --output-base-dir {params.figdir:q} \
        --output-temp-dir {params.tempdir:q}
            """)
        else:
            shell("touch {output.central:q}")
            shell("touch {output.neighbouring:q}")
            shell("touch {output.region:q}")

rule regional_map_rendering:
    input:
        central = os.path.join(config["tempdir"], "central_map_ukLin.vl.json"),
        neighboring = os.path.join(config["tempdir"], "neighboring_map_ukLin.vl.json"),
        region = os.path.join(config["tempdir"], "region_map_ukLin.vl.json")
    params:
        outdir = config["rel_outdir"],
        local_lineages = config["local_lineages"],
        central = os.path.join(config["tempdir"], "central_map_ukLin.vg.json"),
        neighboring = os.path.join(config["tempdir"], "neighboring_map_ukLin.vg.json"),
        region = os.path.join(config["tempdir"], "region_map_ukLin.vg.json")
    output:
        central = os.path.join(config["outdir"], 'figures', "regional_mapping", "central_map_ukLin.png"),
        neighboring = os.path.join(config["outdir"], 'figures', "regional_mapping", "neighboring_map_ukLin.png"),
        region = os.path.join(config["outdir"], 'figures', "regional_mapping", "region_map_ukLin.png")
    run:
        if params.local_lineages:
            shell(
            """
            npx -p vega-lite vl2vg {input.central} {params.central}
            npx -p vega-cli vg2png {params.central} {output.central}
            """)
            shell(
            """
            npx -p vega-lite vl2vg {input.neighboring} {params.neighboring}
            npx -p vega-cli vg2png {params.neighboring} {output.neighboring}
            """)
            shell(
            """
            npx -p vega-lite vl2vg {input.region} {params.region}
            npx -p vega-cli vg2png {params.region} {output.region}
            """)
        else:
            shell("touch {output.central}")
            shell("touch {output.neighboring}")
            shell("touch {output.region}")


rule make_report:
    input:
        lineage_trees = rules.process_catchments.output.tree_summary,
        query = config["query"],
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        cog_global_metadata = config["cog_global_metadata"],
        report_template = config["report_template"],
        polytomy_figure = config["polytomy_figure"],
        footer = config["footer"],
        no_seq = rules.get_closest_cog.output.not_processed,
        clean_locs = config["clean_locs"],
        uk_map = config["uk_map"],
        channels_map = config["channels_map"],
        ni_map = config["ni_map"],
        local_lineages_png = rules.regional_map_rendering.output.outpng,

    params:
        treedir = os.path.join(config["outdir"],"local_trees"),
        outdir = config["rel_outdir"],
        fields = config["fields"],
        sc_source = config["sequencing_centre"],
        sc = config["sequencing_centre_file"],
        sc_flag = config["sequencing_centre_flag"],
        rel_figdir = os.path.join(".","figures"),
        figdir = os.path.join(config["outdir"],"figures"),
        failure = config["qc_fail"],
        local_lineages = config["local_lineages"],
        local_lineages_tables = os.path.join(config["outdir"], 'figures', "regional_mapping", "{wildcards.location}_lineageTable.md")
        #localLinMaps = {input.local_lineages_png},
        #localLinTables = {input.local_lineages_tables}
    output:
        poly_fig = os.path.join(config["outdir"],"figures","polytomies.png"),
        footer_fig = os.path.join(config["outdir"], "figures", "footer.png"),
        outfile = os.path.join(config["outdir"], "civet_report.md")
    run:
        if params.sc != "":
            shell("cp {params.sc_source:q} {params.sc:q}")
        if params.local_lineages:
            #change the relevant bits here
            print({input.local_lineages_png})
            print({input.local_lineages_tables})
            #localLinMaps = ';'.join(list({input.local_lineages_png}))
            #localLinTables = ';'.join({input.local_lineages_tables})
            shell(
            """
                cp {input.polytomy_figure:q} {output.poly_fig:q}
                cp {input.footer:q} {output.footer_fig:q}
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
                --clean-locs {input.clean_locs} \
                --uk-map {input.uk_map} \
                --channels-map {input.channels_map} \
                --ni-map {input.ni_map} \
                --local_lineages

                --outfile {output.outfile:q} \
                --outdir {params.outdir:q} 
                """)
                            #--local_lin_maps {localLinMaps:q} \
                #--local_lin_tables {localLinTables:q} \
        else:
            shell(
            """
                cp {input.polytomy_figure:q} {output.poly_fig:q}
                cp {input.footer:q} {output.footer_fig:q}
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
                --clean-locs {input.clean_locs} \
                --uk-map {input.uk_map} \
                --channels-map {input.channels_map} \
                --ni-map {input.ni_map} \
                --outfile {output.outfile:q} \
                --outdir {params.outdir:q} 
                """)

rule launch_grip:
    input:
        mdfile = os.path.join(config["outdir"], "civet_report.md")
    output:
        out_file = os.path.join(config["outdir"],"civet_report.html")
    run:
        shell("grip {input.mdfile:q} --export")
        for i in range(8000, 8100):
            try:
                shell("grip {input.mdfile:q} -b {i}")
                break
            except:
                print("Trying next port")
