"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module dumps all the protein, paper and mesh data I have in flat files into MARC DB
"""

import cPickle
import DB
import sys
sys.path.append("../objects")
import MeshTerm
DATABASE_NAME = "MarcDB"

def dump_proteins(proteinList, db):
    """ Dump all the protein data I have in proteinList into protein collection in MARC DB using
    DB module
    
    Args:
        proteinList: List of proteins to dump
    """
    for protein in proteinList:
        if proteinList[protein].sequence != "" and type(proteinList[protein].sequence is not str):
                    proteinList[protein].sequence = proteinList[protein].sequence._data
        DB.InsertProtein(proteinList[protein], db)

def dump_papers(paperList, db):
    """ Dump all the paper data I have in paperList into paper collection in MARC DB using
    DB module
    
    Args:
        paperList: List of papers to dump
    """
    for paper in paperList:
        DB.InsertPaper(paperList[paper], db)
        
def dump_MeshTerms(termList, db):
    """ Dump all the mesh data I have in termList into mesh collection in MARC DB using
    DB module
    
    Args:
        termList: List of mesh terms to dump
    """
    for term in termList:
        DB.InsertMeshTerm(termList[term], db)


if __name__ == "__main__":
    """ Main function that shows examples of how to connect to the database and dump protein, papers
    and mesh data:
    
    - First, load the pickle files of the human and mouse proteins, papers, and mesh terms
    
    huProts = cPickle.load(open("data/human_proteins.pik"))
    muProts = cPickle.load(open("data/mouse_proteins.pik"))
    papers  = cPickle.load(open("data/papers.pik"))
    terms   = cPickle.load(open("../data/mesh_data_final.pik"))
    
    - Connect to the database:
    
    dbClient = DB.Connect()
    db = DB.GetDatabase(dbClient, DATABASE_NAME)
    
    - Dump the files:
    
    dump_MeshTerms(terms, db)
    dump_proteins(huProts, db)
    dump_proteins(muProts, db)
    dump_papers(papers, db)
    
    - Load paper citations:
    
    paper_cits = cPickle.load(open("../data/paper_cits.pik"))
    
    - Add citations to database:
    
    DB.update(paper_cits, db)
    """
    