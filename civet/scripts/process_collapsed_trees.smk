"""
Passed into config:

catchment_str=tree_1,tree_2,...,tree_X
outdir=path/to/outdir
not_cog_csv
in_all_cog_fasta
post_qc_query
cog_seqs
combined_metadata
"""
from Bio import Phylo
from Bio import SeqIO
import csv
import collections

config["tree_stems"] = config["catchment_str"].split(",")

rule all:
    input:
        expand(os.path.join(config["outdir"], "collapsed_trees","{tree}.tree"), tree = config["tree_stems"]),
        os.path.join(config["outdir"],"combined_trees","collapse_report.txt"),
        expand(os.path.join(config["outdir"],"restored_trees","{tree}.tree"), tree = config["tree_stems"])

rule summarise_polytomies:
    input:
        tree = os.path.join(config["outdir"], "catchment_trees","{tree}.tree"),
        metadata = config["combined_metadata"]
    params:
        tree_dir = os.path.join(config["outdir"],"catchment_trees")
    output:
        collapsed_tree = os.path.join(config["outdir"],"collapsed_trees","{tree}.tree"),
        collapsed_information = os.path.join(config["outdir"],"collapsed_trees","{tree}.collapsed_info.txt")
    shell:
        """
        clusterfunk focus -i {input.tree:q} \
        -o {output.collapsed_tree:q} \
        --metadata {input.metadata} \
        --index-column closest \
        --in-format newick \
        --out-format newick \
        --output-tsv {output.collapsed_information}
        """

rule get_collapsed_representative:
    input:
        cog_seqs = config["all_cog_seqs"],
        collapsed_information = rules.summarise_polytomies.output.collapsed_information
    params:
        tree_dir = os.path.join(config["outdir"],"collapsed_trees")
    output:
        representative_seq = os.path.join(config["outdir"],"collapsed_trees","{tree}_representatives.fasta"),
    run:
        collapsed = {}
        collapsed_seqs = collections.defaultdict(list)
        
        with open(input.collapsed_information, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                collapsed_name,taxa = l.split('\t')
                collapsed[collapsed_name] = taxa.split(",")
        for record in SeqIO.parse(input.cog_seqs,"fasta"):
            for node in collapsed:
                if record.id in collapsed[node]:
                    collapsed_seqs[node].append(record)

        with open(output.representative_seq, "w") as fw:
            for node in collapsed_seqs:
                records = collapsed_seqs[node]
                sorted_with_amb = []
                for record in records:
                    amb_count = 0
                    for base in record.seq:
                        if base.upper() not in ["A","T","C","G","-"]:
                            amb_count +=1
                    amb_pcent = (100*amb_count) / len(record.seq)
                    sorted_with_amb.append((record.id, amb_pcent, record.seq))
                sorted_with_amb = sorted(sorted_with_amb, key = lambda x : x[1])
                rep = sorted_with_amb[0]
                fw.write(f">{node} representative={rep[0]} ambiguity={rep[1]}\n{rep[2]}\n")
        
rule extract_taxa:
    input:
        collapsed_tree = os.path.join(config["outdir"],"collapsed_trees","{tree}.tree")
    output:
        tree_taxa = os.path.join(config["outdir"], "collapsed_trees","{tree}_taxon_names.txt")
    shell:
        "clusterfunk get_taxa -i {input.collapsed_tree} --in-format newick -o {output.tree_taxa} --out-format newick"

rule gather_fasta_seqs:
    input:
        collapsed_nodes = os.path.join(config["outdir"],"collapsed_trees","{tree}_representatives.fasta"),
        post_qc_query = config["post_qc_query"],
        in_all_cog_fasta = config["in_all_cog_fasta"],
        cog_seqs = config["all_cog_seqs"],
        combined_metadata = config["combined_metadata"],
        tree_taxa = rules.extract_taxa.output.tree_taxa
    output:
        aln = os.path.join(config["outdir"], "catchment_aln","{tree}.query.aln.fasta")
    run:
        taxa = []
        with open(input.tree_taxa, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                taxa.append(l)

        queries = {}
        with open(input.combined_metadata,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row["closest"] in taxa:
                    if row["closest"] != row["query"]:
                        queries[row["query"]] = row["query_id"]

        added_seqs = []
        with open(output.aln, "w") as fw:

            for record in SeqIO.parse(input.post_qc_query, "fasta"):
                if record.id in queries.values() or record.id in queries.keys():
                    iqtree_friendly = record.id.replace("/","_")
                    fw.write(f">{record.id}\n{record.seq}\n")
                    added_seqs.append(record.id)

            for record in SeqIO.parse(input.in_all_cog_fasta, "fasta"):
                if record.id in queries.values() or record.id in queries.keys():
                    iqtree_friendly = record.id.replace("/","_")
                    fw.write(f">{record.id}\n{record.seq}\n")
                    added_seqs.append(record.id)

            for record in SeqIO.parse(input.collapsed_nodes, "fasta"):
                iqtree_friendly = record.id.replace("/","_")
                fw.write(f">{record.id}\n{record.seq}\n")
                added_seqs.append(record.id)

            for record in SeqIO.parse(input.cog_seqs,"fasta"):
                if record.id in taxa:
                    iqtree_friendly = record.id.replace("/","_")
                    fw.write(f">{record.id}\n{record.seq}\n")
                    added_seqs.append(record.id)
        not_added = []
        for seq in taxa:
            if seq not in added_seqs:
                not_added.append(seq)
        print(f"Seqs not added are:")
        for seq in not_added:
            print(f"- {seq}")

rule hash_for_iqtree:
    input:
        aln = rules.gather_fasta_seqs.output.aln
    output:
        hash = os.path.join(config["outdir"], "renamed_trees","{tree}.hash_for_iqtree.csv"),
        hashed_aln = os.path.join(config["outdir"], "renamed_trees","{tree}.query.aln.fasta")
    run:
        fw = open(output.hash, "w")
        fw.write("taxon,iqtree_hash,cluster_hash\n")
        hash_count = 0
        with open(output.hashed_aln, "w") as fseq:
            for record in SeqIO.parse(input.aln, "fasta"):
                hash_count +=1
                without_str = record.id.rstrip("'").lstrip("'")
                fw.write(f"{without_str},_taxon_{hash_count}_,taxon_{hash_count}\n")
                fseq.write(f">'taxon_{hash_count}'\n{record.seq}\n")
        fw.close()

rule hash_tax_labels:
    input:
        tree=os.path.join(config["outdir"],"collapsed_trees","{tree}.tree"),
        hash = rules.hash_for_iqtree.output.hash
    output:
        tree = os.path.join(config["outdir"],"renamed_trees","{tree}.tree")
    shell:
        """
        clusterfunk relabel_tips -i {input.tree} \
        -o {output[0]} \
        --in-metadata {input.hash} \
        --index-column taxon \
        --trait-columns cluster_hash \
        --replace \
        --in-format newick \
        --out-format newick
        """

rule iqtree_catchment:
    input:
        aln = rules.hash_for_iqtree.output.hashed_aln,
        guide_tree = rules.hash_tax_labels.output.tree,
        taxa = rules.extract_taxa.output.tree_taxa
    output:
        tree = os.path.join(config["outdir"], "renamed_trees","{tree}.query.aln.fasta.treefile")
    run:
        taxa = 0
        aln_taxa = 0
        with open(input.taxa, "r") as f:
            for l in f:
                l  = l.rstrip ("\n")
                taxa +=1
        for record in SeqIO.parse(input.aln, "fasta"):
            aln_taxa +=1 
        if taxa != aln_taxa:
            shell("iqtree -s {input.aln:q} -bb 1000 -au -alrt 1000 -m HKY -nt 1 -redo")
        else:
            with open(output.tree,"w") as fw:
                with open(input.guide_tree, "r") as f:
                    for l in f:
                        l = l.rstrip("\n")
                        l = l.replace("'","_")
                        fw.write(l)

rule restore_tip_names:
    input:
        tree = rules.iqtree_catchment.output.tree,
        hash = rules.hash_for_iqtree.output.hash
    output:
        os.path.join(config["outdir"],"restored_trees","{tree}.tree")
    shell:
        """
        clusterfunk relabel_tips -i {input.tree} \
        -o {output[0]} \
        --in-metadata {input.hash} \
        --index-column iqtree_hash \
        --trait-columns taxon \
        --replace \
        --in-format newick \
        --out-format newick
        """

rule summarise_processing:
    input:
        collapse_reports = expand(os.path.join(config["outdir"],"collapsed_trees","{tree}.collapsed_info.txt"), tree=config["tree_stems"])
    output:
        report = os.path.join(config["outdir"],"combined_trees","collapse_report.txt")
    run:
        with open(output.report, "w") as fw:
            for report in input.collapse_reports:
                fn = os.path.basename(report)
                with open(report, "r") as f:
                    for l in f:
                        l = l.rstrip("\n")
                        new_l = f"{fn}\t{l}\n"
                        fw.write(new_l)
