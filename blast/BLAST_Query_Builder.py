"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module creates fasta files from a list of proteins
"""
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import cPickle as cp


def create_FASTA_file(proteins):
    """ Create a fasta file from a dictionary of mouse proteins and their corresponding sequences and store
    it on disk under the name mouse.fasta
    
    Args:
        proteins: A dictionary that maps the mouse protein ids to their sequences
    
    Returns:
        None
    """
    frecords = []
    count = 0
    for prot in proteins:
        if proteins[prot][1] != "":
            prot_id = prot
            sequence = str(proteins[prot][1])
            frec = SeqRecord(Seq(sequence, IUPAC.protein), id=prot_id, description="")
            frecords.append(frec)
        else:
            count += 1
    output_handle = open("BLAST/mouse.fasta", "w")
    SeqIO.write(frecords, output_handle, "fasta")
    output_handle.close()
    print "count = ", count

if __name__ == "__main__":
    """ This is the main function that loads the mouse protein dictionary into memory and calls the fasta
    file creating function

    mouse_prots = cp.load(open("mouse_proteins.pik"))
    create_FASTA_file(mouse_prots)
    """