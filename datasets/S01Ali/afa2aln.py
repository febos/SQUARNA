from Bio import SeqIO
import glob

for file in glob.glob("afa_ungap/*"):
    alnfile = file.replace("afa_ungap/","aln_ungap/").replace(".afa",".aln")
    records = SeqIO.parse(file, "fasta")
    count = SeqIO.write(records, alnfile, "clustal")
    print("Converted {}".format(file))
