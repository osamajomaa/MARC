"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module contains the class that represents a protein. It contains 6 data members:
    1) The protein's id (pid)
    2) The protein's sequence (seq)
    3) The protein's organism (organism)
    4) The protein's list of GO terms that it's annotated to (go_terms)
    5) The list of papers that study this protein (papers)
    6) The list of the protein's homologs (homologs)
"""


class Protein:
    def __init__(self, pid, organism, seq, go_terms, papers, homologs):
        self.pid = pid
        self.organism = organism
        self.sequence = seq
        self.go_terms = go_terms
        self.papers = papers
        self.homologs = homologs
    
    