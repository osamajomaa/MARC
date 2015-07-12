"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This is a parser for the blast XML hit files. It extracts the homologs for the query mouse proteins,
and store them in a dictionary that maps each mouse protein id to it's list of homologs. 
"""
from Bio.Blast import NCBIXML
import cPickle
import os

RESULTS_PATH = "BLAST/Results"

def get_blast_hits():
    """ Parses a group of blast hit XML files in the BLAST/Results folder to extract the list of 
    homologs for each queried mouse protein and stores them in a dictionary and writes it to a file
    under the name mouse_homologs.pik
    
    Args:
        None
    
    Returns:
        None
    """
    files = os.listdir(RESULTS_PATH)
    homologs = {}
    count = 1
    for fname in files:
        result_handle = open(os.path.join(RESULTS_PATH, fname))
        records = NCBIXML.parse(result_handle)
        for rec in records:
            query = rec.query.split()[0]
            hits = []
            if rec.alignments:
                rank = 1
                for alignment in rec.alignments:
                    protein = alignment.title.split('|')[2].split()[1]
                    hits.append((protein, rank))
                    rank += 1
            print "Done with ", count
            homologs[query] = hits
    fHandler = open("mouse_homologs.pik", 'w')
    cPickle.dump(homologs, fHandler)
    fHandler.close()

if __name__ == "__main__":
    """ This is the main function that calls the hit parser for the list of blast hit XML files
    
    get_blast_hits()
    """