"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module contains a class that represents the mesh term. It contains 5 data members:
    1) The term's list of identifiers, where an id is just a concatenation of the ids of the terms acestors (tid)
    2) The term's name which uniquily identifies that term (termName)
    3) The term's list of categories, where a category is one of the 16 top level categories that the term falls under (categroy)
    4) The list of the term's parent terms (parent) 
    5) The list of the term's ancestors (ancestors)
"""

class MeshTerm:
    def __init__(self, tid, termName, category, parent, ancestors):
        self.tid = tid
        self.name = termName
        self.category = [category]
        self.parent = [parent]
        self.ancestors = [ancestors]
    
    