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
        iqtree  -s {input.aln:q} \
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
        tree = os.path.join(config["data_outdir"],"catchments","{catchment}.tree")
        # tree = os.path.join(config["tempdir"],"catchments","{catchment}.pruned.tree")
    shell:
        """
        jclusterfunk prune  -i {input.tree:q} \
                            -o {output.tree:q} \
                            -t outgroup \
                            -f newick 
        """

# rule expand_hash:
#     input:
#         tree = rules.prune_outgroup.output.tree,
#         csv = config["csv"]
#     output:
#         tree = os.path.join(config["tempdir"],"catchments","{catchment}.expanded.tree")
#     shell:
#         """
#         jclusterfunk insert --destination-column hash \
#                             -c {config[report_column]} \
#                             -i {input.tree:q} \
#                             -m {input.csv:q} \
#                             --ignore-missing \
#                             -f newick \
#                             -o {output.tree}
#         """

# rule prune_hashed_seqs:
#     input:
#         tree = rules.expand_hash.output.tree,
#         prune = config["csv"]
#     output:
#         tree = os.path.join(config["data_outdir"],"catchments","{catchment}.tree")
#     run:
#         """
#         jclusterfunk prune  -i {input.tree:q} \
#                             -o {output.tree:q} \
#                             -m {input.prune:q} \
#                             -c hash \
#                             -f newick 
#         """
