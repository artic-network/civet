import os
import csv

catchments = [f"catchment_{i}" for i in range(1,config["catchment_count"]+1)]

rule all:
    input:
        expand(os.path.join(config["outdir"],"catchments","{catchment}.tree"), catchment=catchments)

rule iqtree:
    input:
        aln = os.path.join(config["data_outdir"],"catchments","{catchment}.fasta")
    output:
        tree = os.path.join(config["tempdir"],"catchments","{catchment}.fasta.treefile")
    log:
        os.path.join(config["tempdir"],"catchments","{catchment}.log")
    shell:
        """
        iqtree  -s {input.aln:q} \
                -m HKY \
                -blmin  0.0000000001 \
                -nt 1 \
                -redo \
                --fast \
                -o outgroup \
                -mem {config[max_memory]}G \
                -quiet &> {log:q}
        """

rule prune_outgroup:
    input:
        tree = rules.iqtree.output.tree,
        prune = config["outgroup_fasta"]
    output:
        tree = os.path.join(config["tempdir"],"catchments","{catchment}.og_prune.tree")
    shell:
        """
        jclusterfunk prune  -i {input.tree:q} \
                            -o {output.tree:q} \
                            -t outgroup \
                            -f newick 
        """

rule expand_hash:
    input:
        tree = rules.prune_outgroup.output.tree,
        csv = config["csv"]
    output:
        tree = os.path.join(config["tempdir"],"catchments","{catchment}.expanded.tree")
    shell:
        """
        jclusterfunk insert --destination-column hash \
                            -c {config[report_column]} \
                            -i {input.tree:q} \
                            -m {input.csv:q} \
                            --ignore-missing \
                            -f newick \
                            -o {output.tree}
        """

rule prune_hashed_seqs:
    input:
        tree = rules.expand_hash.output.tree,
    output:
        tree = os.path.join(config["tempdir"],"catchments","{catchment}.hashed_prune.tree")
    run:
        hash_strings = []
        with open(config["csv"],"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                hash_strings.append(row["hash"])
        hash_strings = ' '.join(list(set(hash_strings)))

        shell("jclusterfunk prune  -i {input.tree:q} "
                           " -o {output.tree:q} "
                           f" -t '{hash_strings}' "
                           " --ignore-missing "
                           " -c hash "
                           " -f newick ")

rule clump:
    input:
        tree = rules.prune_hashed_seqs.output.tree,
        csv = config["csv"]
    params:
        prefix = "{catchment}",
        outdir = os.path.join(config["data_outdir"],"catchments")
    output:
        tree = os.path.join(config["data_outdir"],"catchments","{catchment}_tree.nexus")
    shell:
        """
        jclusterfunk sample -c {config[background_column]} \
                            --collapse-by {config[location_column]} \
                            --min-clumped 4 \
                             -i {input.tree:q} \
                             -m {input.csv:q} \
                             -p {params.prefix:q} \
                             -f nexus \
                             -o {params.outdir:q} \
                            --ignore-missing
        """

rule annotate:
    input:
        tree = rules.clump.output.tree,
        csv = config["csv"]
    output:
        tree = os.path.join(config["outdir"],"catchments","{catchment}.tree")
    shell:
        """
        jclusterfunk annotate -c {config[background_column]} \
                             -i {input.tree:q} \
                             -m {input.csv:q} \
                             --tip-attributes query_boolean {config[tree_annotations]} \
                             -f nexus \
                             -o {output.tree:q} \
                            --ignore-missing
        """