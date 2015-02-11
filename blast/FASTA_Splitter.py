'''
    This file is to break down a fasta file into a number of subfiles so that it can be processed 
    in parallel
'''
from Bio import SeqIO
import os

def fasta_chunker(n, fname):
    records = list(SeqIO.parse(fname, "fasta"))
    length = len(records)
    csize = length/n
    rem = length%n
    start = 0
    end = csize + rem
    i = 1
    while (end <= length):
        chunk = records[start:end]
        output_handle = open(os.path.join("BLAST/Mouse_FASTA","part"+str(i)+".fasta"), "w")
        SeqIO.write(chunk, output_handle, "fasta")
        output_handle.close()
        start = end
        end = start + csize
        i += 1

if __name__ == "__main__":
    fastafile = "BLAST/mouse.fasta"
    chunks = 10
    fasta_chunker(chunks, fastafile)
    print "Finish!"
