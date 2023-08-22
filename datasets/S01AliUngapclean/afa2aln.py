from Bio import SeqIO
import glob

for file in glob.glob("afa_clean/*"):
    alnfile = file.replace("afa_clean/","aln/").replace(".afa",".aln")
    records = SeqIO.parse(file, "fasta")
    count = SeqIO.write(records, alnfile, "clustal")
    print("Converted {}".format(file))
