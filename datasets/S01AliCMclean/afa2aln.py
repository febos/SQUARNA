from Bio import SeqIO
import glob

for file in glob.glob("afa/*"):
    alnfile = file.replace("afa/","aln/").replace(".afa",".aln")
    records = SeqIO.parse(file, "fasta")
    count = SeqIO.write(records, alnfile, "clustal")
    print("Converted {}".format(file))
