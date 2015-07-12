"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module creates a BLAST database, creates a BLAST query and executes 
"""

from Bio.Blast.Applications import NcbiblastpCommandline
from subprocess import call
import sys
import os

# A constant variable tha represents the e-value for BLAST queries
E_VALUE = 10**-8

def build_BLAST_db(fasta_file, blast_db):
    """ Creates a blast database from a fasta file
    
    Args:
        fasta_file: The name and location of the fasta file in a string format
        blast_db: The desired name and location for the newly created BLAST database in string format
    
    Returns:
        None
    """
    cmd= "makeblastdb -in " + fasta_file + " -dbtype 'prot' -out " + blast_db
    call(cmd, shell=True)
    
def make_BLAST_query(query_file, db_file, out_file):
    """ Creates and executes a BLAST query which queries the query_file against the BLAST database in db_file
    and stores the blast hits in out_file
    
    Args:
        query_file: The name and location of the fasta file that contains the proteins and their sequences in query
        db_file: The name and location of the BLAST datbase
        out_file: The name and location of the output BLAST files
    
    Returns
        None
    """
    cline = NcbiblastpCommandline(query=query_file, db=db_file, outfmt=5, evalue=E_VALUE, out=out_file)
    print (cline)
    cline()

if __name__ == "__main__":
    """ This is the main function that contains the function calls to create BLAST database for 
    the mouse and human proteins and quering these database using a protein fasta files
    
    Example of creating the human BLAST database from a human fasta file:
    
    fasta_file = "BLAST/human.fasta"
    blast_db = "BLAST/DB/humanDB"
    build_BLAST_db(fasta_file, blast_db)
    
    Example of using a mouse fasta file to query a human BLAST DB and storing the BLAST results in
    results.xml file:
    
    query_file = "BLAST/mouse.fasta"
    db_file = "BLAST/DB/humanDB"
    out_file = "results.xml"
    make_BLAST_query(query_file, db_file, out_file)

    Example of using command line arguments to query a BLAST DB instead of hard-coded arguments
    in the above example:
    
    query_file = sys.argv[1]
    db_file = sys.argv[2]
    out_file = sys.argv[3]
    make_BLAST_query(query_file, db_file, out_file)
    """
    