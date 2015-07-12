""" 
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module contains functions necessary to connect to MARC DB.
"""


import pymongo
from sets import Set

def Connect():
    """ Connects to a mongo database by reading the following parameters from utils/db_creds.txt file:
        1) Server's address (1st line)
        2) Server's port (2nd line)
        3) Database name (3rd line)
        4) Database Username (4th line)
        5) Database password (5th line)
    
    Args:
        None
    
    Returns:
        Database client instance, which can be used to get an instance of the database
    """
    f = open("../utils/db_creds.txt")
    server = f.readline().strip()
    port = f.readline().strip()
    dbName = f.readline().strip()
    username = f.readline().strip()
    password = f.readline().strip()
    uri = "mongodb://"+username+":"+password+"@"+server+":"+port+"/"+dbName
    return pymongo.MongoClient(uri)

def GetDatabase(Client, dbName):
    """ Get an instance of MARC DB
    
    Args:
        Client: The database client returned from the Connect() function
        dbName: The name of the database (MarcDB)
    
    Returns:
        Database instance
    """
    return Client[dbName]

def InsertProtein(protein, db):
    """ Inserts a protein into the Protein collection in MARC DB
    
    Args:
        protein: an object of the Protein class
        db: an instance of MARC DB
    
    
    Returns:
        None
    """
    protColl = db.Protein
    homologDoc = []
    termDoc = []
    for homo in protein.homologs:
        homologDoc.append({"Homolog": homo[0], "Rank": homo[1]})
    for term in protein.go_terms:
        termDoc.append({"GOID": term[0], "Namespace": term[1], "Name": term[2], "Description":term[3]})
    protDoc = {"PID":       protein.pid, 
               "Organism":  protein.organism,
                "Sequence": protein.sequence,
               "GO_Terms":  termDoc,
               "Papers":    protein.papers, 
               "Homologs":  homologDoc}
    protColl.insert(protDoc)
    

def InsertPaper(paper, db):
    """ Inserts a paper into the Paper collection in MARC DB
    
    Args:
        paper: an object of the paper class
        db: an instance of MARC DB
    
    
    Returns:
        None
    """
    paperColl = db.Paper
    meshDoc = []
    for heading in paper.meshHeadings:
        meshDoc.append({"Descriptor": heading[0], "Qualifiers":heading[1]})
    paperDoc = {"PMID":         paper.pmid,
                "Organism":     paper.organism, 
                "Citations":    paper.citation,
                "MeshHeadings": meshDoc,
                "PubTypes":     paper.pubTypes}
    paperColl.insert(paperDoc)


def InsertMeshTerm(term, db):
    """ Inserts a mesh term into the MeshTerm collection in MARC DB
    
    Args:
        term: an object of the MeshTerm class
        db: an instance of MARC DB
    
    
    Returns:
        None
    """
    termColl = db.MeshTerm
    catDoc= []
    i = 1
    for cat in term.category:
        catDoc.append({"Branch": i, "Category":cat})
        i += 1
    i = 1
    parentDoc= []
    for parent in term.parent:
        parentDoc.append({"Branch": i, "Parent":parent})
        i += 1
    i = 1
    ancestorDoc= []
    for ancs in term.ancestors:
        for anc in ancs:
            ancestorDoc.append({"Branch": i, "Ancestors":anc})
        i += 1
    i = 1
    tids = []
    for tid in term.tid:
        tids.append({"Branch":i, "TID":tid})
        i += 1
        
    termDoc = {"TID":                   tids, 
               "Name":                  term.name,
               "Categories":            catDoc,
                "Parents":              parentDoc,
               "AncestorsList":        ancestorDoc
               }
    termColl.insert(termDoc)

def getAll(projects, collName, db):
    """ Get all the documents in a database collection
    
    Args:
        collName: The name of the collection to get the documents from
        db: An instance of the database
    
    Returns:
        A list of all the documents in the collectiuon collName
    """
    coll = db[collName]
    data = coll.find({})
    return data


def find(projects, collName, query, db):
    """ Find a document in a collection according to some query
    
    Args:
        projects: A dictionary that contains the projections, i.e. the fields to be returns from the documents
        collName: The name of the collection to search in
        query: The NoSQL query which is a python dictionary form the search criteria
        db: An instance of the database
    
    Returns:
        A list of the documents that matches the search query in the collection collName
    """
    data = []
    coll = db[collName]
    rowData = coll.find(query, projects)
    for rowDatum in rowData:
        datum = {}
        for proj in projects:
            datum[proj] = rowDatum[proj]
        data.append(datum)
    

def update(paper_cits, db):
    """ Update the list of citations for each paper in the Paper collection
    
    Args:
        paper_cits: A dictionary that maps each paper with its list of citations
        db: An instance of the database
    
    Returns:
        None
    """
    coll = db["Paper"]
    for paper in paper_cits.keys():
        coll.update({"PMID":paper}, {"$set": {"CitedBy":paper_cits[paper]}})


def findByAncestor(anc, db):
    """ Get the list of descendants for a given ancestor anc
    
    Args:
        anc: The mesh term for the ancestor
        db: An instance of the database
    
    Returns:
        A list of the descendant mesh terms of the ancestor anc
    """
    coll = db["MeshTerm"]
    data = coll.find({"Parents.Parent":{"$in":anc}})
    terms = Set()
    for datum in data:
        terms.add(datum["Name"])
    return terms


def findTermsAtDepth(db, depth, cat):
    """ Get the list of mesh terms at a depth greater than or equal to a given depth for a given mesh category
    
    Args:
        depth: The minimum depth that the mesh terms have to at least have
        cat: The name of the mesh category
        db: An instance of the database
    
    Returns:
        The list of mesh terms at a depth greater than or equal to the depth (depth) and for category (cat)
    """
    coll = db["MeshTerm"]
    terms = coll.find({})
    data = Set()
    for term in terms:
        for tid in term["TID"]:
            if tid["TID"][:3] == cat and tid["TID"].count(".") >= depth-1:
                data.add(term["Name"])
    return data



    