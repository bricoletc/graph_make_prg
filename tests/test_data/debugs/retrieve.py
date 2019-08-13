from Bio import SeqIO

with open("AMA_Plasmodium_cycle.fasta","w") as f:

    for record in SeqIO.parse("AMA_Plasmodium.fasta","fasta"):
        f.write(">" + record.id + "\n")
        f.write(str(record.seq[277:285]) + "\n")


