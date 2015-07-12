"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module includes a function that maps each paper to the papers that cite it
"""

import cPickle as cp

def getCitedBy(paperList):
    """ Map each paper in the paper list to the papers that cite it and write the resulting 
    dictionary to paper_cits.pik file
    
    Args:
        paperList: List of papers
    
    Returns:
        None
    """
    paper_cits = {}
    for paper in paperList:
        paper_cits[paper] = []
    for paper in paperList:
        for ref in paperList[paper].citation:
            if ref in paper_cits:
                paper_cits[ref].append(paper)
    output_handle = open("../data/paper_cits.pik", "w")
    cp.dump(paper_cits, output_handle)
    output_handle.close()


if __name__ == "__main__":
    """ This is the main function that contains example function calls
    
    papers = cp.load(open("../data/papers.pik"))
    getCitedBy(papers)
    """