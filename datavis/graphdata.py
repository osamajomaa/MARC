"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module includes functions that create graph files in GEXF format. Cytoscape does not accept files
in this format, therefore these files can be opened using Gefi software, converted to graphml format, 
and then rendered using Cytoscape.
"""

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

#Constant variable which is the name of the database
DATABASE_NAME = "MarcDB"


def getMeshAncestors(term, db):
    """ Returns the list of ancestor terms and their categories for an imput mesh term
    
    Args:
        term (String): The name of the mesh term
        db (Database): The database class instance
        
    Returns:
        List of Strings: List of the names of the mesh categories that the term ancestors fall under
        List of Strings: List of the names of the mesh ancestors
    """
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
    """ Creates an XML tag element in the GEXF file
    
    Args:
        parent (Element): The parent tag that contains the new tag element
        title (String): The title of the new tag element
        text (String): The text of the new tag element, if any.
        attributes (List of Attribute): List of the attribute pairs (name, value) that the new tag element has
        
    Returns:
        Element: The newly created XML tag element
    """
    element = SubElement(parent, title, attributes)
    if text != '':
        element.text = text
    return element

def createGEXFHeader():
    """ Creates the Header tag for XML file which contains version, encoding and metadata information
    including the name of the graph creator and a brief description.
    
    Args:
        None
    
    Returns:
        Element: The root element of the XML file
    """
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
    """ Generates a dictionary of the 16 top-level mesh categories with their corresponding number values
    
    Args:
        None
    
    Returns:
        Dictionary: mapping of the mesh 16 top-level terms with their number values
    """
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
    """ Returns a list of the mesh terms that are at depth <= depth
    
    Args:
        depth (int): a numeric value of the max depth
        cat (String): the name of the mesh cat
    
    Returns:
        List of Strings: List of the names of the mesh terms at depth <= depth
    """
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    res = dbase.findTermsAtDepth(db, depth, cat)
    return res
    

def writeXML(root, output):
    """Writes an XML instance object to file
    
    Args:
        root (Element): The instance object that contains the whole XML structure in memory
        output (String): The name of the output file to store the XML structure
    """

    f = open(output, 'w')
    f.write(tostring(root))
    f.close()  


def reHash(meshData):
    data = {}
    for term in meshData:
        data[term] = meshData[term]
    return data


def cited(paper_id, pmid):
    """ Check if a a paper which id is pmid is cited by the paper which has the paper_id
    
    Args:
        paper_id (String): The pmid of the citing paper
        pmid (String): The pmid of the potential cited paper
    
    Returns:
        bool: True if pmid is cited by paper_id, false otherwise.
    """
    for idd in paper_id:
        if pmid in paper_id[idd][1]:
            return True
    return False



def createCatFile(depth, cat, fileName):
    """Generates a dictionary that maps the mesh categories that are at depth <= depth with the list 
    of papers that are annotated to them and writes the dictionary to a pickled file named fileName.
    
    Args:
        depth (int): minimum depth that the subcategory must have from the its parent category cat
        cat (String): The name of the parent category
        fileName: The name of the file to write the dictionary to
    
    Returns:
        None
    """
    meshData = cPickle.load(open("../data/mesh_data2.pik"))
    paper_cats = OrderedDict()
    pathConds = getTermsAtDepth(depth, cat)
    print len(pathConds)
    for cond in pathConds:
        paper_cats[cond] = Set()
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    papers = dbase.getAll({},"Paper", db)
    for paper in papers:
        if "MeshHeadings" in paper:
            for term in paper["MeshHeadings"]:
                if term["Descriptor"] in meshData and term["Descriptor"] in pathConds:
                    paper_cats[term["Descriptor"]].add(paper["PMID"])
    
    paper_count = OrderedDict()
    for cat in paper_cats:
        paper_count[cat] = len(paper_cats[cat])
    
    paper_count = OrderedDict(sorted(paper_count.iteritems(), key=operator.itemgetter(1), reverse=True))
    
    fHandler = open(fileName, 'w')
    cPickle.dump(paper_count, fHandler)
    fHandler.close()
    

def createGEXFWithDepth(catFile, top, fileName, edge_file):
    """ Creates a GEXF graph file that contains the citation networks. This graph contains:
    1) Metadata at the root element
    2) Attribute List
    3) Node List
    4) Edge List
    
    This function also creates a dictionary that maps each category in the catFile with a list of tuples
    (5 element structure) which contains:
    (0, 1, 2, 3, 4) = (total number of citations in category cat,
                                        number of intra citations in category cat,
                                        number of inter-citing citations between category cat and other category,
                                        number of inter-cited citations between category cat and other category,
                                        number of isolated nodes annotated to category cat,
                                        total number of papers annotated to category cat)
    It also writes this dictionary on file named edge_file
    
    Args:
        catFile: The name of the mesh categories-papers dictionary file
        top: The number of the categories that have the largest number of papers annotated to them
        fileName: The name of the GEXF graph file to be stored on disk
        edge_file: The name of the output file that contains the categories stats
    
    Returns:
        None
    """
    meshData = cPickle.load(open("../data/mesh_data2.pik"))
    root = createGEXFHeader()
    cfile = cPickle.load(open(catFile))
    pathConds = cPickle.load(open(catFile)).keys()[:top]
    graph = createXMLElement(root, 'graph', '', {'defaultedgetype':'directed'})
    attributes = createXMLElement(graph, 'attributes', '', {'class':'node'})
    createXMLElement(attributes, 'attribute', '', {'id':'2', 'title':'Descriptors', 'type':'string'})
    createXMLElement(attributes, 'attribute', '', {'id':'3', 'title':'MeshCat', 'type':'string'})
    dbClient = dbase.Connect()
    db = dbase.GetDatabase(dbClient, DATABASE_NAME)
    papers = dbase.getAll({},"Paper", db)
    nodes = createXMLElement(graph, 'nodes', '', {})
    count = 0
    paper_id = {}
    cited_nodes = Set()
    for paper in papers:
        descStr= ""
        if "MeshHeadings" in paper:
            for term in paper["MeshHeadings"]:
                if term["Descriptor"] in meshData and term["Descriptor"] in pathConds:
                    descStr += term["Descriptor"] + '&'
        if len(descStr) > 0:
            descStr = descStr[:-1]
            cited_nodes.update(paper["Citations"])
            if '&' not in descStr:
                paper_id[paper["PMID"]] = (str(count), paper["Citations"], "sin", descStr)
            else:
                paper_id[paper["PMID"]] = (str(count), paper["Citations"], "mul", descStr)
            node = createXMLElement(nodes, 'node', '', {"id":str(count), "label":paper["PMID"]})
            attvalues = createXMLElement(node, 'attvalues', '', {})
            createXMLElement(attvalues, 'attvalue', '', {"for":"2", "value": descStr})
            if '&' in descStr:
                createXMLElement(attvalues, 'attvalue', '', {"for":"3", "value": "Many"})
            else:
                createXMLElement(attvalues, 'attvalue', '', {"for":"3", "value": descStr})
            print count
            count = count + 1
        
    edges = createXMLElement(graph, 'edges', '', {})
    count = 0
    cat_record = {}
    for cat in pathConds:
        cat_record[cat] = [0,0,0,0,0,0]
    for pmid in paper_id:
        source = paper_id[pmid][0]
        isolated = True
        for cit in paper_id[pmid][1]:
            if cit in paper_id:
                isolated = False
                target = paper_id[cit][0]
                for cat in pathConds:
                    if cat in paper_id[pmid][3] or cat in paper_id[cit][3]:
                        if cat in paper_id[pmid][3] and cat in paper_id[cit][3]:
                            cat_record[cat][1] += 1
                        elif cat in paper_id[pmid][3]:
                            cat_record[cat][2] += 1
                        elif cat in paper_id[cit][3]:
                            cat_record[cat][3] += 1
                        cat_record[cat][0] += 1                        
                if paper_id[pmid][3] == paper_id[cit][3] and paper_id[pmid][2]=="sin" and paper_id[cit][2]=="sin":                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                    createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":target, "weight":"1.0"})
                elif paper_id[pmid][3] != paper_id[cit][3] and paper_id[pmid][2]=="sin" and paper_id[cit][2]=="sin":                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                    createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":target, "weight":"2.0"})
                else:
                    createXMLElement(edges, 'edge', '', {"id":str(count), "source":source, "target":target, "weight":"3.0"})
                count += 1
        for cat in pathConds:
            if cat in paper_id[pmid][3]:
                cat_record[cat][5] += 1
                if isolated and pmid not in cited_nodes:
                    cat_record[cat][4] += 1 
    writeXML(root,fileName)
    fHandler = open(edge_file, 'w')
    cPickle.dump(cat_record, fHandler)
    fHandler.close()
    

if __name__ == "__main__":
    """ This main function contain example function calls to create a category, graph and
    edge stats files for neoplasm category
    
    createCatFile(3, "C04", "../data/FamousCats/neoplasmOrdered.pik")
    createGEXFWithDepth("../data/FamousCats/neoplasmOrdered.pik", 10, "../data/FamousCats/GraphFiles/Neoplasm.gexf", "../data/Edge_Desc/neoplasm_edge_desc.pik")
    
    Other examples for creating GEXF files for other disease sub-categories:
    
    #createGEXFWithDepth("../data/FamousCats/animalDiseaseOrdered.pik", 10, "../data/FamousCats/GraphFiles/animalDisease.gexf", "../data/Edge_Desc/animal_edge_desc.pik")
    #createGEXFWithDepth("../data/FamousCats/cogenitalOrdered.pik", 10, "../data/FamousCats/GraphFiles/cogenitalDisease.gexf", "../data/Edge_Desc/congenital_edge_desc.pik")
    #createGEXFWithDepth("../data/FamousCats/digestiveOrdered.pik", 10, "../data/FamousCats/GraphFiles/digestiveDisease.gexf", "../data/Edge_Desc/digestive_edge_desc.pik")
    #createGEXFWithDepth("../data/FamousCats/neuroSystemOrdered.pik", 10, "../data/FamousCats/GraphFiles/neuroSystemDisease.gexf", "../data/Edge_Desc/neuro_edge_desc.pik")
    #createGEXFWithDepth("../data/FamousCats/nutritionalOrdered.pik", 10, "../data/FamousCats/GraphFiles/nutritionalDisease.gexf", "../data/Edge_Desc/nutritional_edge_desc.pik")
"""
    