import cPickle
import DB
DATABASE_NAME = "MarcDB"

def dump_proteins(proteinList, db):
    for protein in proteinList:
        '''
        if i == 25:
            print "i == 25"
        seq = proteinList[protein][1]
        if seq != "" and type(proteinList[protein][1] is not str):
            seq = proteinList[protein][1]._data
        else:
            print seq
        p = Protein.Protein(protein,
                            proteinList[protein][0],
                            seq,
                            proteinList[protein][2],
                            proteinList[protein][3],
                            proteinList[protein][4])
        '''
        if proteinList[protein].sequence != "" and type(proteinList[protein].sequence is not str):
                    proteinList[protein].sequence = proteinList[protein].sequence._data
        DB.InsertProtein(proteinList[protein], db)

def dump_papers(paperList, db):
    for paper in paperList:
        DB.InsertPaper(paperList[paper], db)
        
def dump_MeshTerms(termList, db):
    for term in termList:
        DB.InsertMeshTerm(termList[term], db)


if __name__ == "__main__":
    #huProts = cPickle.load(open("data/human_proteins.pik"))
    #mesh = cPickle.load(open("../data/mesh_data2.pik"))
    #huProts = cPickle.load(open("data/human_proteins.pik"))
    #papers  = cPickle.load(open("data/papers.pik"))
    #terms   = cPickle.load(open("data/mesh_data.pik"))
    
    dbClient = DB.Connect()
    db = DB.GetDatabase(dbClient, DATABASE_NAME)
    paper_cits = cPickle.load(open("../data/paper_cits.pik"))
    DB.update(paper_cits, db)
    #dump_MeshTerms(mesh, db)
    #dump_proteins(huProts, db)
    #dump_papers(papers, db)
    #dump_MeshTerms(terms, db)
    