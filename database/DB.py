import pymongo
from sets import Set

def Connect():
    f = open("../utils/db_creds.txt")
    server = f.readline().strip()
    port = f.readline().strip()
    dbName = f.readline().strip()
    username = f.readline().strip()
    password = f.readline().strip()
    uri = "mongodb://"+username+":"+password+"@"+server+":"+port+"/"+dbName
    return pymongo.MongoClient(uri)

def GetDatabase(Client, dbName):
    return Client[dbName]

def InsertProtein(protein, db):
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
    coll = db[collName]
    data = coll.find({})
    return data


def find(projects, collName, query, db):
    data = []
    coll = db[collName]
    rowData = coll.find(query, projects)
    for rowDatum in rowData:
        datum = {}
        for proj in projects:
            datum[proj] = rowDatum[proj]
        data.append(datum)
    

def update(paper_cits, db):
    coll = db["Paper"]
    i = 1
    for paper in paper_cits.keys():
        coll.update({"PMID":paper}, {"$set": {"CitedBy":paper_cits[paper]}})
        print i
        i += 1

def findByAncestor(anc, db):
    coll = db["MeshTerm"]
    data = coll.find({"Parents.Parent":{"$in":anc}})
    terms = Set()
    for datum in data:
        terms.add(datum["Name"])
    return terms

def findTermsAtDepth(db, depth, cat):
    coll = db["MeshTerm"]
    terms = coll.find({})
    data = Set()
    for term in terms:
        for tid in term["TID"]:
            if tid["TID"][:3] == cat and tid["TID"].count(".") >= depth-1:
                data.add(term["Name"])
    return data



    