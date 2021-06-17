import os

catchments = [f"catchment_{i}" for i in range(1,config["catchment_count"]+1)]

rule all:
    input:
        expand(os.path.join(config["data_outdir"],"catchments","{catchment}.tree"), catchment=catchments)

rule iqtree:
    input:
        aln = os.path.join(config["data_outdir"],"catchments","{catchment}.fasta")
    output:
        tree = os.path.join(config["tempdir"],"catchments","{catchment}.fasta.treefile")
    log:
        os.path.join(config["tempdir"],"catchments","{catchment}.log")
    shell:
        """
        iqtree -s {input.aln:q} \
        -m HKY \
        -blmin  0.0000000001 \
        -nt 1 \
        -redo \
        -o outgroup \
        -quiet &> {log:q}
        """

rule prune_outgroup:
    input:
        tree = rules.iqtree.output.tree,
        prune = config["outgroup_fasta"]
    output:
        tree = os.path.join(config["data_outdir"],"catchments","{catchment}.pruned.tree")
    shell:
        """
        jclusterfunk prune -i {input.tree:q} \
        -o {output.tree:q} \
        -t outgroup \
        -f newick 
        """

rule expand_hash:
    input:
        tree = rules.prune_outgroup.output.tree,
        csv = config["csv"]
    output:
        tree = os.path.join(config["data_outdir"],"catchments","{catchment}.tree")
    shell:
        """
        jclusterfunk insert -c hash \
                            -n {config[input_column]} \
                            -i {input.tree:q} \
                            -m {input.csv:q} \
                            --ignore-missing \
                            -f newick \
                            -o {output.tree}
        """
