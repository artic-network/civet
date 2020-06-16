import os
from Bio import SeqIO
rule all:
    input:
        os.path.join(config["tempdir"],"closest_cog.csv"),
        os.path.join(config["tempdir"],"to_find_closest.fasta")

rule combine_in_all_cog_and_query:
    input:
        fasta = config["post_qc_query"],
        in_all_cog_seqs = config["in_all_cog_seqs"]
    output:
        fasta = os.path.join(config["tempdir"],"to_find_closest.fasta")
    run:
        num_query = 0
        num_all_cog = 0
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input.fasta, "fasta"):
                num_query +=1
                fw.write(f">{record.id} status=query_sequence\n{record.seq}\n")
            for record in SeqIO.parse(input.in_all_cog_seqs, "fasta"):
                num_all_cog +=1
                fw.write(f">{record.id} status=in_all_cog\n{record.seq}\n")
        print(f"{num_query} sequences from the query fasta\n{num_all_cog} found un-analysed in COG database\n")


rule non_cog_minimap2_to_reference:
    input:
        fasta = rules.combine_in_all_cog_and_query.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = os.path.join(config["tempdir"],"post_qc_query.reference_mapped.sam")
    shell:
        """
        minimap2 -a -x asm5 {input.reference:q} {input.fasta:q} > {output.sam:q}
        """

rule non_cog_remove_insertions_and_trim_and_pad:
    input:
        sam = rules.non_cog_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = os.path.join(config["tempdir"],"post_qc_query.insertions.txt")
    output:
        fasta = os.path.join(config["tempdir"],"post_qc_query.alignment.trimmed.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam:q} \
          -r {input.reference:q} \
          -o {output.fasta:q} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts &
        mv insertions.txt {params.insertions}
        """

rule minimap2_against_cog:
    input:
        query_seqs = rules.combine_in_all_cog_and_query.output.fasta,
        cog_seqs = config["cog_seqs"]
    output:
        paf = os.path.join(config["tempdir"],"post_qc_query.cog_mapped.paf")
    shell:
        """
        minimap2 -x asm5 --secondary=no --paf-no-hit {input.cog_seqs:q} {input.query_seqs:q} > {output.paf:q}
        """

rule parse_paf:
    input:
        paf = rules.minimap2_against_cog.output.paf,
        metadata = config["cog_metadata"],
        fasta = config["cog_seqs"]
    params:
        search_field = config["search_field"]
    output:
        fasta = os.path.join(config["tempdir"],"closest_cog.fasta"),
        csv = os.path.join(config["tempdir"],"closest_cog.csv")
    shell:
        """
        parse_paf.py \
        --paf {input.paf:q} \
        --metadata {input.metadata:q} \
        --seqs {input.fasta:q} \
        --csv-out {output.csv:q} \
        --seqs-out {output.fasta:q} \
        --search-field {params.search_field}
        """
