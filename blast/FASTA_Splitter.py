"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module is to break down a fasta file into a number of subfiles so that it can be processed in parallel
"""
from Bio import SeqIO
import os

def fasta_chunker(n, fname):
    """ Split a fasta file into n files that can be processed in parallel using PBS
    and store the chunked files in BLAST/Mouse_FASTA folder
    
    Args:
        n: The number of files that we desire to split the original fasta file into
        fname: The name and location of the input big fasta file
    
    Returns:
        None
    """
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
    """ The main function that contains the number of chunks and the function call for the fasta
    file splitter
    
    Example of the call:
    
    fastafile = "BLAST/mouse.fasta"
    chunks = 10
    fasta_chunker(chunks, fastafile)
    """
