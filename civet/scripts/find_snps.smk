from Bio import Phylo
from Bio import SeqIO
import csv
import collections

config["tree_stems"] = config["catchment_str"].split(",")

rule all:
    input:
        expand(os.path.join(config["outdir"],"figures","genome_graph_{tree}.png"), tree=config["tree_stems"]),
        os.path.join(config["outdir"],"gather_prompt.txt")

rule extract_taxa:
    input:
        collapsed_tree = os.path.join(config["outdir"],"local_trees","{tree}.tree")
    output:
        tree_taxa = os.path.join(config["tempdir"], "local_trees_seqs","{tree}_taxon_names.txt")
    shell:
        "clusterfunk get_taxa -i {input.collapsed_tree:q} --in-format nexus -o {output.tree_taxa:q} --out-format newick"

rule get_sequence_names:
    input:
        tree_taxa = rules.extract_taxa.output.tree_taxa,
        combined_metadata = config["combined_metadata"],
        query = config["query"]
    output:
        seq_names = os.path.join(config["tempdir"], "local_trees_seqs", "{tree}_filtered_to_query.txt")
    run:
        display_name = config["display_name"]
        input_column = config["input_column"]
        data_column = config["data_column"]


        # in input file find input_column and display name
        queries = {}
        with open(input.query,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                queries[row[input_column]] = row[display_name]

        # in combined metadata, match the input column to the data column
        # add sequence name -> display name mapping
        name_maps = {}
        with open(input.combined_metadata, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row[data_column] in queries:
                    name_maps[row["sequence_name"]] = queries[row[data_column]]
        
        # if an input column not found in combined metadata,
        # add the input_column -> display name mapping
        # these would be the ones not in the tree
        for query in queries:
            if query not in name_maps.values():
                name_maps[query] = queries[query]

        # get all the taxa from a given local tree
        taxa = []
        with open(input.tree_taxa, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                taxa.append(l)

        with open(output.seq_names, "w") as fw:
            fw.write("name,label\n")
            fw.write("Reference,Reference\n")
            for name in name_maps:
                if name in taxa:

                    fw.write(f"{name},{name_maps[name]}\n")

rule gather_fasta_seqs:
    input:
        aligned_query_seqs = config["aligned_query_seqs"],
        background_seqs = config["background_seqs"],
        outgroup_fasta = config["outgroup_fasta"],
        seq_names = rules.get_sequence_names.output.seq_names
    output:
        aln = os.path.join(config["tempdir"], "seqs_for_snps","{tree}.fasta")
    run:
        taxa = []
        with open(input.seq_names, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                taxa.append(row["name"])

        with open(output.aln, "w") as fw:
            for record in SeqIO.parse(input.outgroup_fasta, "fasta"):
                fw.write(f">Reference\n{record.seq}\n")

            for record in SeqIO.parse(input.aligned_query_seqs, "fasta"):
                if record.id in taxa:
                    fw.write(f">{record.description}\n{record.seq}\n")

            for record in SeqIO.parse(input.background_seqs,"fasta"):
                if record.id in taxa:
                    fw.write(f">{record.description}\n{record.seq}\n")

rule make_snp_figure:
    input:
        aln = rules.gather_fasta_seqs.output.aln,
        names = rules.get_sequence_names.output.seq_names
    params:
        out_stem = os.path.join(config["outdir"],"figures","genome_graph_{tree}")
    output:
        os.path.join(config["outdir"],"figures","genome_graph_{tree}.png")
    shell:
        """
        snipit {input.aln:q} -r "Reference" -o {params.out_stem} -l {input.names}
        """

rule gather_graphs:
    input:
        expand(os.path.join(config["outdir"],"figures","genome_graph_{tree}.png"), tree=config["tree_stems"])
    output:
        os.path.join(config["outdir"],"gather_prompt.txt")
    shell:
        "touch {output}"
    