
import os
from Bio import SeqIO
import csv

rule all:
    input:
        os.path.join(config["tempdir"],"closest_cog.csv"),
        os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")


rule non_cog_minimap2_to_reference:
    input:
        fasta = config["to_find_closest"],
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
        insertions = os.path.join(config["tempdir"],"post_qc_query.insertions.txt")
    output:
        fasta = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam:q} \
          -r {input.reference:q} \
          -o {output.fasta:q} \
          -t [{config[trim_start]}:{config[trim_end]}] \
          --pad \
          --log-inserts 
        """

rule minimap2_against_cog:
    input:
        query_seqs = rules.non_cog_remove_insertions_and_trim_and_pad.output.fasta,
        cog_seqs = config["seq_db"]
    output:
        paf = os.path.join(config["tempdir"],"post_qc_query.cog_mapped.paf")
    shell:
        """
        minimap2 -x asm5  -t {workflow.cores} --secondary=no --paf-no-hit {input.cog_seqs:q} {input.query_seqs:q} > {output.paf:q}
        """

rule parse_paf:
    input:
        paf = rules.minimap2_against_cog.output.paf,
        metadata = config["cog_metadata"],
        fasta = config["seq_db"]
    output:
        fasta = os.path.join(config["tempdir"],"closest_cog.fasta"),
        csv = os.path.join(config["tempdir"],"closest_cog.no_snps.csv")
    shell:
        """
        parse_paf.py \
        --paf {input.paf:q} \
        --metadata {input.metadata:q} \
        --seqs {input.fasta:q} \
        --csv-out {output.csv:q} \
        --seqs-out {output.fasta:q} \
        --search-field {config[search_field]}
        """

rule snp_diff:
    input:
        closest_fasta = os.path.join(config["tempdir"],"closest_cog.fasta"),
        input_fasta = rules.non_cog_remove_insertions_and_trim_and_pad.output.fasta,
        csv = os.path.join(config["tempdir"],"closest_cog.no_snps.csv")
    output:
        csv = os.path.join(config["tempdir"],"closest_cog.csv")
    run:
        input_seqs = {}
        source_dict = {}
        for record in SeqIO.parse(input.input_fasta, "fasta"):
            input_seqs[record.id] = record
            desc = record.description.split(" ")
            source = ""
            for info in desc:
                if info.startswith("status="):
                    source = info.split("=")[1]
            source_dict[record.id]=source
        
        closest_map = {}
        for record in SeqIO.parse(input.closest_fasta, "fasta"):
            desc = record.description.split(" ")
            for info in desc:
                if info.startswith("query="):
                    queries = info.split("=")[1]
                    for q in queries.split(","):
                        closest_map[q] = record
        non_amb = ["a","t","g","c"]
        with open(input.csv, newline="") as f:
            reader = csv.DictReader(f)
            header_names = reader.fieldnames

            with open(output.csv, "w") as fw:
                header_names.append("closest_distance")
                header_names.append("snps")
                header_names.append("source")
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()
            
                for row in reader:
                    query = row["query_id"]
                    snps =[]
                    q_record = input_seqs[query]
                    c_record = closest_map[query]

                    for i in range(len(q_record.seq)):
                        bases = [q_record.seq[i],c_record.seq[i]]
                        if bases[0] != bases[1]:
                            if bases[0].lower() in non_amb and bases[1].lower() in non_amb:
                                
                                snp = f"{i+1}{bases[0]}{bases[1]}"
                                snps.append(snp)

                    print(query, c_record.id)
                    print(snps)
                    
                    new_row = row
                    new_row["closest_distance"]= f"{len(snps)}"
                    if len(snps) ==0:
                        new_row["snps"]= ""
                    else:
                        snp_str = ";".join(snps)
                        new_row["snps"]= snp_str
                    new_row["source"] = source_dict[query]
                    writer.writerow(new_row)
