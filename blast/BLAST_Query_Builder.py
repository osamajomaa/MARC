from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline
import os
import cPickle as cp


def create_FASTA_file(proteins):
    
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
    
    mouse_prots = cp.load(open("mouse_proteins.pik"))
    create_FASTA_file(mouse_prots)
