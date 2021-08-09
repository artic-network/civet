from Bio import SeqIO

fasta = "/home/mwatson/COVID-19/outbreak/P.1_Potential_July_28_2021/focal_aln.fasta"

fasta_sequences = SeqIO.parse(open(fasta), 'fasta')

first = "yes"
ref_name = "MN908947"

record_names = []
with open("focal_sequence_names.txt", "w") as handle:
    handle.write(f"name,label\n")
    for record in fasta_sequences:
        record_names.append(record.name)
        if record.name in ref_name and first != "yes":
            handle.write(f"{record.name},Reference\n")
        elif first == "yes":
            handle.write(f"{record_names[0]},Reference\n")
            first = "no"
        else:
            handle.write(f"{record.name},{record.name}\n")
