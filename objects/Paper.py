"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module contains a class that represents a mouse or human paper. It contains 5 data members:
    1) The paper's pubmed id (pmid)
    2) The organism that the paper studies (organism)
    3) The list of mesh terms that the paper is annotated to in pubmed (meshHeadings)
    4) The publication type of the paper, for example: journal, review, conference proceedings, etc (pubTypes)
"""
class Paper:
    def __init__(self, pmid, organism, citations, meshHeadings, pubTypes):
        self.pmid = pmid
        self.organism = organism
        self.citation = citations
        self.meshHeadings = meshHeadings
        self.pubTypes = pubTypes
    
    