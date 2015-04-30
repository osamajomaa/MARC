import sys
sys.path.append("../objects")
sys.path.append("../database")
import Paper
import operator
import DB as dbase
from collections import OrderedDict

from sets import Set
import os
from xml.etree.ElementTree import Element, SubElement, tostring
import cPickle


DATABASE_NAME = "MarcDB"

def getMeshAncestors(term, db):
    projects = {"AncestorsList":1, "Categories":1}
    query = {"Name":term}
    data = dbase.find(projects, "MeshTerm", query, db)
    categories = []
    ancestors = []
    for row in data:
        for cat in row["Categories"]:
            categories.append(cat["Category"])
        for anc in row["AncestorsList"]:
            for a in anc["Ancestors"]:
                projects = {"Name":1}
                query = {"TID":a}
                name = ""
                names = dbase.find(projects, "MeshTerm", query, db)
                for n in names:
                    name = n["Name"]
                ancestors.append(name)
    return categories, ancestors




def createXMLElement(parent, title, text, attributes):
    element = SubElement(parent, title, attributes)
    if text != '':
        element.text = text
    return element

def createGEXFHeader():
    root = Element('gexf')
    root.set('version', '1.0')
    root.set('encoding', 'UTF-8')
    meta = SubElement(root, 'meta')
    creator = SubElement(meta, 'creator')
    creator.text = "Osama Jomaa"
    description = SubElement(meta, 'description')
    description.text = "MARC Papers Graph"    
    return root

def getMeshCats():
    meshCats = {}
    meshCats["Anatomy"] = "0"
    meshCats["Organisms"] = "1"
    meshCats["Diseases"] = "2"
    meshCats["Chemicals and Drugs"] = "3"
    meshCats["Analytical,Diagnostic and Therapeutic Techniques and Equipment"] = "4"
    meshCats["Psychiatry and Psychology"] = "5"
    meshCats["Phenomena and Processes"] = "6"
    meshCats["Disciplines and Occupations"] = "7"
    meshCats["Anthropology,Education,Sociology and Social Phenomena"] = "8"
    meshCats["Technology,Industry,Agriculture"] = "9"
    meshCats["Humanities"] = "10"
    meshCats["Information Science"] = "11"
    meshCats["Named Groups"] = "12"
    meshCats["Health Care"] = "13"
    meshCats["Publication Characteristics"] = "14"
    meshCats["Geographicals"] = "15"
    return meshCats

def getpathCondsCats(ancs):
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    res = dbase.findByAncestor(ancs, db)
    return res

def getTermsAtDepth(depth, cat):
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    res = dbase.findTermsAtDepth(db, depth, cat)
    return res
    

def writeXML(root, output):
    
    #xml = minidom.parseString(tostring(root))
    f = open(output, 'w')
    #f.write(xml.toprettyxml())
    f.write(tostring(root))
    f.close()  


def reHash(meshData):
    data = {}
    for term in meshData:
        data[term] = meshData[term]
    return data


def cited(paper_id, pmid):
    for idd in paper_id:
        if pmid in paper_id[idd][1]:
            return True
    return False



def createCatFile(depth, cat, fileName):
    meshData = cPickle.load(open("../data/mesh_data2.pik"))
    paper_cats = OrderedDict()
    pathConds = getTermsAtDepth(depth, cat)
    print len(pathConds)
    for cond in pathConds:
        paper_cats[cond] = Set()
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    papers = dbase.getAll({},"Paper", db)
    count = 0
    for paper in papers:
        #descSet = Set()
        if "MeshHeadings" in paper:
            for term in paper["MeshHeadings"]:
                if term["Descriptor"] in meshData and term["Descriptor"] in pathConds:
                    paper_cats[term["Descriptor"]].add(paper["PMID"])
        print count
        count += 1;
    
    paper_count = OrderedDict()
    for cat in paper_cats:
        paper_count[cat] = len(paper_cats[cat])
    
    paper_count = OrderedDict(sorted(paper_count.iteritems(), key=operator.itemgetter(1), reverse=True))
    
    fHandler = open(fileName, 'w')
    cPickle.dump(paper_count, fHandler)
    fHandler.close()
    print "Finish!"
    

def createGEXFWithDepth(catFile, top, fileName):
    meshData = cPickle.load(open("../data/mesh_data2.pik"))
    root = createGEXFHeader()
    cfile = cPickle.load(open(catFile))
    pathConds = cPickle.load(open(catFile)).keys()[:top]
    graph = createXMLElement(root, 'graph', '', {'defaultedgetype':'directed'})
    attributes = createXMLElement(graph, 'attributes', '', {'class':'node'})
    createXMLElement(attributes, 'attribute', '', {'id':'2', 'title':'Descriptors', 'type':'string'})
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    papers = dbase.getAll({},"Paper", db)
    nodes = createXMLElement(graph, 'nodes', '', {})
    count = 0
    paper_id = {}
    for paper in papers:
        descStr= ""
        if "MeshHeadings" in paper:
            for term in paper["MeshHeadings"]:
                if term["Descriptor"] in meshData and term["Descriptor"] in pathConds:
                    descStr += term["Descriptor"] + '&'
        if len(descStr) > 0:
            descStr = descStr[:-1]
            if '&' not in descStr:
                paper_id[paper["PMID"]] = (str(count), paper["Citations"], "sin", descStr)
            else:
                paper_id[paper["PMID"]] = (str(count), paper["Citations"], "mul", descStr)
            node = createXMLElement(nodes, 'node', '', {"id":str(count), "label":paper["PMID"]})
            attvalues = createXMLElement(node, 'attvalues', '', {})
            createXMLElement(attvalues, 'attvalue', '', {"for":"2", "value": descStr})
            print count
            count = count + 1
        
    edges = createXMLElement(graph, 'edges', '', {})
    count = 0
    for pmid in paper_id:
        source = paper_id[pmid][0]
        for cit in paper_id[pmid][1]:
            if cit in paper_id:
                target = paper_id[cit][0]
                if paper_id[pmid][3] == paper_id[cit][3] and paper_id[pmid][2]=="sin" and paper_id[cit][2]=="sin":                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                    createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":target, "weight":"1.0"})
                elif paper_id[pmid][3] != paper_id[cit][3] and paper_id[pmid][2]=="sin" and paper_id[cit][2]=="sin":                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                    createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":target, "weight":"2.0"})
                else:
                    createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":target, "weight":"3.0"})
                count += 1
            
    writeXML(root,fileName)
    

def createPaperDiseaseGEXF(catFile, top):
    meshData = cPickle.load(open("../data/mesh_data2.pik"))
    root = createGEXFHeader()
    pathConds = getpathCondsCats(["C04"])
    graph = createXMLElement(root, 'graph', '', {'defaultedgetype':'directed'})
    attributes = createXMLElement(graph, 'attributes', '', {'class':'node'})
    createXMLElement(attributes, 'attribute', '', {'id':'1', 'title':'PubTypes', 'type':'string'})
    createXMLElement(attributes, 'attribute', '', {'id':'2', 'title':'Descriptors', 'type':'string'})
    createXMLElement(attributes, 'attribute', '', {'id':'3', 'title':'Qualifiers', 'type':'string'})
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    papers = dbase.getAll({},"Paper", db)
    nodes = createXMLElement(graph, 'nodes', '', {})
    count = 0
    paper_id = {}
    for paper in papers:
        #print paper
        pubTypes = Set()
        if "PubTypes" in paper:
            for ptype in paper["PubTypes"]:
                pubTypes.add(ptype)
        descSet = Set()
        qualSet = Set()
        if "MeshHeadings" in paper:
            for term in paper["MeshHeadings"]:
                if term["Descriptor"] in meshData:
                    descSet.add(str(term["Descriptor"]))
                    ancD = Set(meshData[term["Descriptor"]].ancestors)
                    [descSet.add(each) for each in ancD]
                for qual in term["Qualifiers"]:
                    qualSet.add(str(qual))
        
        interDesc = descSet & pathConds
        if len(interDesc) > 0:
            descStr = ""
            pubTypeStr = ""
            qualStr= ""
            for desc in interDesc:
                descStr += desc + '&'
            for qual in qualSet:
                qualStr += qual + '&'
            for pub in pubTypes:
                pubTypeStr += pub + '&'
            
            if len(descStr) > 0:
                descStr = descStr[:-1]
            if len(pubTypeStr) > 0:
                pubTypeStr = pubTypeStr[:-1];
            if len(qualStr) > 0:
                qualStr = qualStr[:-1]

            paper_id[paper["PMID"]] = (str(count), paper["Citations"])
            node = createXMLElement(nodes, 'node', '', {"id":str(count), "label":paper["PMID"]})
            attvalues = createXMLElement(node, 'attvalues', '', {})
            createXMLElement(attvalues, 'attvalue', '', {"for":"1", "value":pubTypeStr})
            createXMLElement(attvalues, 'attvalue', '', {"for":"2", "value":descStr})
            createXMLElement(attvalues, 'attvalue', '', {"for":"3", "value":qualStr})
            print count
            count = count + 1
        
    edges = createXMLElement(graph, 'edges', '', {})
    count = 0
    for pmid in paper_id:
        source = paper_id[pmid][0]
        for cit in paper_id[pmid][1]:
            if cit in paper_id:                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":paper_id[cit][0]})
                count += 1
            
    writeXML(root,'../data/Graphs/neoplasms.gexf')




'''
def createPaperDiseaseGEXF():
    meshData = cPickle.load(open("../data/mesh_data2.pik"))
    all_papers = cPickle.load(open("../data/papers.pik")).keys()
    #meshData = reHash(meshData)
    root = createGEXFHeader()
    meshCats = getMeshCats()
    graph = createXMLElement(root, 'graph', '', {'defaultedgetype':'directed'})
    attributes = createXMLElement(graph, 'attributes', '', {'class':'node'})
    for categ in meshCats:
        attribute = createXMLElement(attributes, 'attribute', '', {'id':meshCats[categ], 'title':categ, 'type':'boolean'})
        createXMLElement(attribute, 'default', 'false', {})
    createXMLElement(attributes, 'attribute', '', {'id':'16', 'title':'PubTypes', 'type':'string'})
    createXMLElement(attributes, 'attribute', '', {'id':'17', 'title':'Descriptors', 'type':'string'})
    createXMLElement(attributes, 'attribute', '', {'id':'18', 'title':'Qualifiers', 'type':'string'})
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    papers = dbase.getAll({},"Paper", db)
    nodes = createXMLElement(graph, 'nodes', '', {})
    count = 0
    
    paper_id = {}
    for paper in papers:
        cits = False
        for cit in paper["Citations"]:
            if cit in all_papers:
                cits = True
                break
        
        if len(paper["CitedBy"]) == 0 and cits == False:
            continue
        print "Citations = " + str(cits)
        print "CitedBy = " + str(len(paper["CitedBy"]))
        pubTypes = ";"
        descriptors = ";"
        qualifiers = ";"
        if "PubTypes" in paper:
            for ptype in paper["PubTypes"]:
                pubTypes += ptype+";"
        catSet = Set()
        descSet = Set()
        qualSet = Set()
        if "MeshHeadings" in paper:
            for term in paper["MeshHeadings"]:
                #for x in meshData:
                 #   print meshData[x]
                if term["Descriptor"] in meshData:
                    ancD = Set(meshData[term["Descriptor"]].ancestors)
                    catsD = Set(meshData[term["Descriptor"]].category)
                    #catsD, ancD = getMeshAncestors(term["Descriptor"], db)
                    [catSet.add(each) for each in catsD]
                    [descSet.add(each) for each in ancD]
                for qual in term["Qualifiers"]:
                    qualSet.add(qual)
                    
        for desc in descSet:
            if desc != '':
                descriptors += desc+";"
        for qual in qualSet:
            if qual != '':
                qualifiers += qual+";"
        
        if "Neoplasms" in descriptors:
            paper_id[paper["PMID"]] = (str(count), paper["Citations"])
            node = createXMLElement(nodes, 'node', '', {"id":str(count), "label":paper["PMID"]})
            attvalues = createXMLElement(node, 'attvalues', '', {})
            createXMLElement(attvalues, 'attvalue', '', {"for":"16", "value":pubTypes})
            createXMLElement(attvalues, 'attvalue', '', {"for":"17", "value":descriptors})
            createXMLElement(attvalues, 'attvalue', '', {"for":"18", "value":qualifiers})
            for cat in catSet:
                createXMLElement(attvalues, 'attvalue', '', {"for":meshCats[cat], "value":"true"})
            print count
            count = count + 1
        
    edges = createXMLElement(graph, 'edges', '', {})
    count = 0
    for pmid in paper_id:
        source = paper_id[pmid][0]
        for cit in paper_id[pmid][1]:
            if cit in paper_id:                                                                                                                                                                                                                                                                                                                                                                                                                                                             
                createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":paper_id[cit][0]})
                count += 1
            
    writeXML(root,'No_Singletons/neoplasms.gexf')


def createPapersGEXF():
    meshData = cPickle.load(open("../data/mesh_data2.pik"))
    #meshData = reHash(meshData)
    root = createGEXFHeader()
    meshCats = getMeshCats()
    graph = createXMLElement(root, 'graph', '', {'defaultedgetype':'directed'})
    attributes = createXMLElement(graph, 'attributes', '', {'class':'node'})
    for categ in meshCats:
        attribute = createXMLElement(attributes, 'attribute', '', {'id':meshCats[categ], 'title':categ, 'type':'boolean'})
        createXMLElement(attribute, 'default', 'false', {})
    createXMLElement(attributes, 'attribute', '', {'id':'16', 'title':'PubTypes', 'type':'string'})
    createXMLElement(attributes, 'attribute', '', {'id':'17', 'title':'Descriptors', 'type':'string'})
    createXMLElement(attributes, 'attribute', '', {'id':'18', 'title':'Qualifiers', 'type':'string'})
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    papers = dbase.getAll({},"Paper", db)
    nodes = createXMLElement(graph, 'nodes', '', {})
    count = 0
    
    paper_id = {}
    for paper in papers:
        paper_id[paper["PMID"]] = (str(count), paper["Citations"])
        pubTypes = ";"
        descriptors = ";"
        qualifiers = ";"
        for ptype in paper["PubTypes"]:
            pubTypes += ptype+";"
        catSet = Set()
        descSet = Set()
        qualSet = Set()
        for term in paper["MeshHeadings"]:
            #for x in meshData:
             #   print meshData[x]
            if term["Descriptor"] in meshData:
                ancD = Set(meshData[term["Descriptor"]].ancestors)
                catsD = Set(meshData[term["Descriptor"]].category)
                #catsD, ancD = getMeshAncestors(term["Descriptor"], db)
                [catSet.add(each) for each in catsD]
                [descSet.add(each) for each in ancD]
            for qual in term["Qualifiers"]:
                qualSet.add(qual)

        for desc in descSet:
            if desc != '':
                descriptors += desc+";"
        for qual in qualSet:
            if qual != '':
                qualifiers += qual+";"
        
        node = createXMLElement(nodes, 'node', '', {"id":str(count), "label":paper["PMID"]})
        attvalues = createXMLElement(node, 'attvalues', '', {})
        createXMLElement(attvalues, 'attvalue', '', {"for":"16", "value":pubTypes})
        createXMLElement(attvalues, 'attvalue', '', {"for":"17", "value":descriptors})
        createXMLElement(attvalues, 'attvalue', '', {"for":"18", "value":qualifiers})
        for cat in catSet:
            createXMLElement(attvalues, 'attvalue', '', {"for":meshCats[cat], "value":"true"})
        print count
        count = count + 1
        
    edges = createXMLElement(graph, 'edges', '', {})
    count = 0
    for pmid in paper_id:
        source = paper_id[pmid][0]
        for cit in paper_id[pmid][1]:
            if cit in paper_id:                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
                createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":paper_id[cit][0]})
                count += 1
            ../data/
    writeXML(root,'output.xml')
'''

if __name__ == "__main__":
    #print os.getcwd()
    #createCatFile(3, "C16", "../data/FamousCats/cogenitalOrdered.pik")
    #createGEXFWithDepth("../data/FamousCats/cogenitalOrdered.pik", 10, "../data/FamousCats/GraphFiles/cogenitalDisease.gexf")
    createGEXFWithDepth("../data/FamousCats/digestiveOrdered.pik", 10, "../data/FamousCats/GraphFiles/digestiveDisease.gexf")
    #createGEXFWithDepth("../data/FamousCats/neuroSystemOrdered.pik", 10, "../data/FamousCats/GraphFiles/neuroSystemDisease.gexf")
    #createGEXFWithDepth("../data/FamousCats/nutritionalOrdered.pik", 10, "../data/FamousCats/GraphFiles/nutritionalDisease.gexf")
    
    