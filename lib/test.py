import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
file1=sys.argv[1]
def get_chr_length(file_name):
    fa_dict = SeqIO.to_dict(SeqIO.parse(file_name, "fasta", IUPAC.unambiguous_dna))
    for i in fa_dict.keys():
        chr_length=i+'\t'+str(len(fa_dict[i]))
    return chr_length
get_chr_length(file1)