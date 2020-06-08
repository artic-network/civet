import os

rule all:
    input:
        os.path.join(config["outdir"],"closest_cog.csv")

rule non_cog_match_fasta:
    input:
        fasta = config["post_qc_query"],
        not_cog = config["not_cog_csv"]
    output:
        fasta = os.path.join(config["outdir"],"not_in_cog.fasta")
    run:
        not_cog = []
        with open(input.not_cog, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                not_cog.append(l)
        in_fasta = []
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input.fasta, "fasta"):
                if record.id in not_cog:
                    fw.write(f">{record.id}\n{record.seq}\n")
                    in_fasta.append(record.id)
                else:
                    in_fasta.append(record.id)
        print("The following sequences queried not in COG \nand not in fasta file so cannot be analysed\n")
        for record in not_cog:
            if record not in in_fasta: 
                print(record)
                
rule non_cog_minimap2_to_reference:
    input:
        fasta = rules.non_cog_match_fasta.output.fasta,
        reference = config["reference_fasta"]
    output:
        sam = os.path.join(config["outdir"],"post_qc_query.reference_mapped.sam")
    shell:
        """
        minimap2 -a -x asm5 {input.reference} {input.fasta} > {output.sam}
        """

rule non_cog_remove_insertions_and_trim_and_pad:
    input:
        sam = rules.non_cog_minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = os.path.join(config["outdir"],"post_qc_query.insertions.txt")
    output:
        fasta = os.path.join(config["outdir"],"post_qc_query.alignment.trimmed.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam} \
          -r {input.reference} \
          -o {output.fasta} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts &
        mv insertions.txt {params.insertions}
        """

rule minimap2_against_cog:
    input:
        query_seqs = config["post_qc_query"],
        cog_seqs = config["cog_seqs"]
    output:
        paf = os.path.join(config["outdir"],"post_qc_query.cog_mapped.paf")
    shell:
        """
        minimap2 -x asm5 --secondary=no --paf-no-hit {input.cog_seqs} {input.query_seqs} > {output.paf}
        """

rule parse_paf:
    input:
        paf = rules.minimap2_against_cog.output.paf,
        metadata = config["cog_metadata"],
        fasta = config["cog_seqs"]
    output:
        fasta = os.path.join(config["outdir"],"closest_cog.fasta"),
        csv = os.path.join(config["outdir"],"closest_cog.csv")
    shell:
        """
        parse_paf.py \
        --paf {input.paf:q} \
        --metadata {input.metadata:q} \
        --seqs {input.fasta} \
        --csv-out {output.csv:q} \
        --seqs-out {output.fasta:q}
        """
