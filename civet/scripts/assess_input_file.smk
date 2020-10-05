import csv
from Bio import SeqIO
import os
import collections
import snakemake
import sys
import yaml
from reportfunk.funks import custom_logger as custom_logger


from reportfunk.funks import io_functions as qcfunk
import civetfunks as cfunk

rule check_cog_db:
    input:
        query = config["query"],
        background_seqs = config["background_seqs"],
        background_metadata = config["background_metadata"]
    output:
        cog = os.path.join(config["tempdir"],"query_in_cog.csv"),
        background_seqs = os.path.join(config["tempdir"],"query_in_cog.fasta"),
        not_cog = os.path.join(config["tempdir"],"not_in_cog.csv")
    shell:
        """
        check_cog_db.py --query {input.query:q} \
                        --cog-seqs {input.background_seqs:q} \
                        --cog-metadata {input.background_metadata:q} \
                        --field {config[data_column]} \
                        --in-metadata {output.cog:q} \
                        --in-seqs {output.background_seqs:q} \
                        --input-column {config[input_column]} \
                        --not-in-cog {output.not_cog:q} 
        """

rule get_closest_cog:
    input:
        snakefile = os.path.join(workflow.current_basedir,"find_closest_cog.smk"),
        reference_fasta = config["reference_fasta"],
        background_seqs = config["background_seqs"],
        background_metadata = config["background_metadata"]
    output:
        closest_cog = os.path.join(config["tempdir"],"closest_cog.csv"),
        aligned_query = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")
    run:
        if config["fasta"] != "":
            if config["num_seqs"] != 0:
                num_seqs = config["num_seqs"]
                print(qcfunk.green(f"Passing {num_seqs} sequences into search pipeline:"))

                for record in SeqIO.parse(config["post_qc_query"], "fasta"):
                    print(f"    - {record.id}")

                shell("snakemake --nolock --snakefile {input.snakefile:q} "
                            "{config[force]} "
                            "{config[log_string]}"
                            "--directory {config[tempdir]:q} "
                            "--config "
                            "tempdir={config[tempdir]:q} "
                            "background_metadata={input.background_metadata:q} "
                            "background_seqs={input.background_seqs:q} "
                            "to_find_closest='{config[post_qc_query]}' "
                            "data_column={config[data_column]} "
                            "trim_start={config[trim_start]} "
                            "trim_end={config[trim_end]} "
                            "reference_fasta={input.reference_fasta:q} "
                            "--cores {workflow.cores}")
        else:
            shell("touch {output.closest_cog:q} && touch {output.aligned_query:q} && echo 'No closest sequences to find'")

rule combine_metadata:
    input:
        closest_cog = rules.get_closest_cog.output.closest_cog,
        in_cog = rules.check_cog_db.output.cog
    output:
        combined_csv = os.path.join(config["outdir"],"combined_metadata.csv")
    run:
        c = 0
        with open(input.in_cog, newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames

            with open(output.combined_csv, "w") as fw:
                header_names.append("SNPdistance")
                header_names.append("SNPs")
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()
            
                for row in reader:
                    c +=1
                    new_row = row
                    new_row["SNPdistance"]="0"
                    new_row["SNPs"]= ""

                    writer.writerow(new_row)

                with open(input.closest_cog, newline="") as fc:
                    readerc = csv.DictReader(fc)
                    for row in readerc:
                        writer.writerow(row)
                        c+=1
        if c ==0:
            sys.stderr.write(qcfunk.cyan(f'Error: no valid querys to process\n'))
            sys.exit(-1)

rule prune_out_catchments:
    input:
        tree = config["background_tree"],
        metadata = rules.combine_metadata.output.combined_csv
    params:
        outdir = os.path.join(config["outdir"],"catchment_trees")
    output:
        summary = os.path.join(config["outdir"],"catchment_trees", "tree_collapsed_nodes.csv")
    run:
        shell(
        """
        jclusterfunk context \
        -i "{input.tree}" \
        -o "{params.outdir}" \
        --max-parent {config[up_distance]} \
        --max-child {config[down_distance]} \
        -f newick \
        -p tree_ \
        --ignore-missing \
        -m "{input.metadata}" \
        --id-column closest 
        """)

rule process_catchments:
    input:
        snakefile_collapse_before = os.path.join(workflow.current_basedir,"process_collapsed_trees.smk"),
        snakefile_just_collapse = os.path.join(workflow.current_basedir,"just_collapse_trees.smk"),
        combined_metadata = rules.combine_metadata.output.combined_csv, 
        query_seqs = rules.get_closest_cog.output.aligned_query, #datafunk-processed seqs
        collapse_summary = rules.prune_out_catchments.output.summary,
        background_seqs = config["background_seqs"],
        outgroup_fasta = config["outgroup_fasta"]
    params:
        tree_dir = os.path.join(config["tempdir"],"catchment_trees")
    output:
        tree_summary = os.path.join(config["outdir"],"local_trees","collapse_report.txt")
    run:
        catchment_trees = []
        for r,d,f in os.walk(params.tree_dir):
            for fn in f:
                if fn.endswith(".newick"):
                    file_stem = ".".join(fn.split(".")[:-1])
                    catchment_trees.append(file_stem)
        catchment_str = ",".join(catchment_trees) #to pass to snakemake pipeline

        query_seqs = 0
        for record in SeqIO.parse(input.query_seqs,"fasta"):
            query_seqs +=1

        if query_seqs !=0:
            
            snakestring = f"'{input.snakefile_collapse_before}' "
            print(f"Processing catchment trees")
            shell(f"snakemake --nolock --snakefile {snakestring}"
                        "{config[force]} "
                        "{config[log_string]} "
                        "--directory {config[tempdir]:q} "
                        "--config "
                        f"catchment_str={catchment_str} "
                        "outdir={config[outdir]:q} "
                        "tempdir={config[tempdir]:q} "
                        "outgroup_fasta={input.outgroup_fasta:q} "
                        "aligned_query_seqs={input.query_seqs:q} "
                        "background_seqs={input.background_seqs:q} "
                        "combined_metadata={input.combined_metadata:q} "
                        "collapse_summary={input.collapse_summary:q} "
                        "collapse_threshold={config[collapse_threshold]} "
                        "--cores {workflow.cores}")
        else:
            print(f"No new sequences to add in, just collapsing trees")
            shell("snakemake --nolock --snakefile {input.snakefile_just_collapse:q} "
                            "{config[force]} "
                            "{config[log_string]} "
                            "--directory {config[tempdir]:q} "
                            "--config "
                            f"catchment_str={catchment_str} "
                            "outdir={config[outdir]:q} "
                            "tempdir={config[tempdir]:q} "
                            "collapse_threshold={config[collapse_threshold]} "
                            "collapse_summary={input.collapse_summary:q} "
                            "combined_metadata={input.combined_metadata:q} "
                            "--cores {workflow.cores}")

rule find_snps:
    input:
        tree_summary = os.path.join(config["outdir"],"local_trees","collapse_report.txt"),
        snakefile = os.path.join(workflow.current_basedir,"find_snps.smk"),
        query_seqs = rules.get_closest_cog.output.aligned_query, #datafunk-processed seqs
        background_seqs = config["background_seqs"],
        outgroup_fasta = config["outgroup_fasta"]
    params:
        tree_dir = os.path.join(config["outdir"],"local_trees"),
    output:
        genome_graphs = os.path.join(config["tempdir"],"gather_prompt.txt") 
    run:
        local_trees = []
        for r,d,f in os.walk(params.tree_dir):
            for fn in f:
                if fn.endswith(".tree"):
                    file_stem = ".".join(fn.split(".")[:-1])
                    local_trees.append(file_stem)
        local_str = ",".join(local_trees) #to pass to snakemake pipeline
        temp_output = os.path.join(config["outdir"],"gather_prompt.txt") 
        shell("snakemake --nolock --snakefile {input.snakefile:q} "
                            "{config[force]} "
                            "{config[log_string]} "
                            "--directory {config[tempdir]:q} "
                            "--config "
                            f"catchment_str={local_str} "
                            "outdir={config[outdir]:q} "
                            "tempdir={config[tempdir]:q} "
                            "outgroup_fasta={input.outgroup_fasta:q} "
                            "aligned_query_seqs={input.query_seqs:q} "
                            "background_seqs={input.background_seqs:q} "
                            "collapse_threshold={config[collapse_threshold]} "
                            "--cores {workflow.cores} ")
        if not config["no_temp"]:
            shell(f"&& mv '{temp_output}' "
            "{config[tempdir]:q}")

rule regional_mapping:
    input:
        query = config['query'],
        combined_metadata = os.path.join(config["outdir"],"combined_metadata.csv"),
        background_metadata = config["background_metadata"]
    params:
        figdir = os.path.join(config["outdir"],'figures'),
    output:
        central = os.path.join(config["tempdir"], "central_map_ukLin.vl.json"),
        neighboring = os.path.join(config["tempdir"], "neighboring_map_ukLin.vl.json"),
        region = os.path.join(config["tempdir"], "region_map_ukLin.vl.json")
    run:
        if config["local_lineages"] == True:
            shell("""
        local_scale_analysis.py \
        --uk-map {config[uk_map_d3]:q} \
        --hb-translation {config[HB_translations]:q} \
        --date-restriction {config[date_restriction]:q} \
        --date-pair-start {config[date_range_start]:q} \
        --date-pair-end {config[date_range_end]:q} \
        --date-window {config[date_window]:q} \
        --cog-meta-global {input.background_metadata:q} \
        --user-sample-data {input.query:q} \
        --combined-metadata {input.combined_metadata:q} \
        --input-name {config[input_column]:q} \
        --output-base-dir {params.figdir:q} \
        --output-temp-dir {config[tempdir]:q}
            """)
        else:
            shell("touch {output.central:q}")
            shell("touch {output.neighboring:q}")
            shell("touch {output.region:q}")

rule regional_map_rendering:
    input:
        central = os.path.join(config["tempdir"], "central_map_ukLin.vl.json"),
        neighboring = os.path.join(config["tempdir"], "neighboring_map_ukLin.vl.json"),
        region = os.path.join(config["tempdir"], "region_map_ukLin.vl.json")
    params:
        central = os.path.join(config["tempdir"], "central_map_ukLin.vg.json"),
        neighboring = os.path.join(config["tempdir"], "neighboring_map_ukLin.vg.json"),
        region = os.path.join(config["tempdir"], "region_map_ukLin.vg.json")
    output:
        central = os.path.join(config["outdir"], 'figures', "central_map_ukLin.png"),
        neighboring = os.path.join(config["outdir"], 'figures', "neighboring_map_ukLin.png"),
        region = os.path.join(config["outdir"], 'figures', "region_map_ukLin.png")
    run:
        if config["local_lineages"] == True:
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
        background_metadata = config["background_metadata"],
        snp_figure_prompt = rules.find_snps.output.genome_graphs,
        genome_graphs = rules.find_snps.output.genome_graphs, 
        central = os.path.join(config["outdir"], 'figures', "central_map_ukLin.png"),
        neighbouring = os.path.join(config["outdir"], 'figures', "neighboring_map_ukLin.png"),
        region = os.path.join(config["outdir"], 'figures', "region_map_ukLin.png")
    output:
        poly_fig = os.path.join(config["outdir"],"figures","polytomies.png"),
        footer_fig = os.path.join(config["outdir"], "figures", "footer.png"),
        yaml = os.path.join(config["outdir"],"config.final.yaml"),
        outfile = os.path.join(config["outdir"], "civet_report.md")
    run:
        
        shell("cp {config[sequencing_centre_source]:q} {config[sequencing_centre_dest]:q}")

        cfunk.local_lineages_to_config(input.central, input.neighbouring, input.region, config)

        config["figdir"] = os.path.join(".","figures") #changed from rel_figdir
        config["treedir"] = os.path.join(config["outdir"],"local_trees")
        config["outfile"] = os.path.join(config["outdir"], "civet_report.md")
        config["summary_dir"] = os.path.join(config["outdir"], "summary_files")
        config["filtered_background_metadata"] = input.combined_metadata
        config["name_stem"] = "civet_report" #may need to change this if we still have different templates - if we don't then we can hard
        qcfunk.get_tree_name_stem(config["treedir"],config)

        with open(output.yaml, 'w') as fw:
            yaml.dump(config, fw) #so at the moment, every config option gets passed to the report

        shell("""
        cp {config[polytomy_figure]:q} {output.poly_fig:q} &&
        cp {config[footer]:q} {output.footer_fig:q} """)
        
        shell("make_report.py --config {output.yaml:q} ")

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