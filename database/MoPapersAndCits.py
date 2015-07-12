"""
Author: Osama Jomaa

Date: 2014-2015

Version: 1.0

This module parses the gene association files and fasta files to extract data about mouse and human 
papers and proteins and store them in flat files
"""


from Bio.UniProt import GOA
import go_obo_parser as gop
from Bio import Entrez
from Bio import SeqIO
import cPickle as cp
import sys
sys.path.append("../objects")
import Mesh_Categories as MC
import MeshTerm as MT
from sets import Set

Entrez.email = "jomaao@miamioh.edu"

def GetMouseProteinData(gaf_file, go_obo_file, uniprot_prots, mouse_homologs):
    """ Extract the papers, go terms, organism and sequence for each mouse protein in a gene association
     file and writes them into mouse_proteins flat file on disk
    
    Args:
        gaf_file: Uniprot_GOA association file in gaf format
        go_obo_file: The file that contains the whole gene ontology in obo format
        uniprot_prots: mouse protein sequences in uniprot database
        mouse_homologs: Mouse homologous proteins from BLAST
    
    Returns:
        None
    """
    unigoa_file = open(gaf_file)
    GO = GetGO(go_obo_file)
    prot_pmid_go = {}
    count = 0
    protCount = 0
    for inrec in GOA.gafiterator(unigoa_file):
        prot = inrec['DB_Object_ID']
        term = inrec['GO_ID']
        namespace = inrec['Aspect']
        seq = ""
        organism = "Mouse"
        go_terms = []
        pmids = []
        homologs = []
        if prot in mouse_homologs:
            homologs = mouse_homologs[prot]
        if prot not in prot_pmid_go:
            prot_pmid_go[prot] = ()
            if prot in uniprot_prots:
                seq = uniprot_prots[prot]
                if seq != "" and type(uniprot_prots[prot] is not str):
                    seq = uniprot_prots[prot]._data
            else:
                protCount += 1
                print prot
        else:
            organism = prot_pmid_go[prot].organism
            seq = prot_pmid_go[prot].sequence
            go_terms = prot_pmid_go[prot].go_terms
            pmids = prot_pmid_go[prot].papers
            homologs = prot_pmid_go[prot].homologs
        if term not in go_terms:
                if term in GO:
                    go_terms.append((term, namespace, GO[term][0], GO[term][1]))
                else:
                    go_terms.append((term, namespace, "", ""))
                    count += 1
        for dbref in inrec['DB:Reference']:
            if dbref[:4] == 'PMID':
                pmid = dbref[5:]
                if pmid not in pmids:
                    pmids.append(pmid)
        prot_pmid_go[prot] = PT.Protein(prot, organism, seq, go_terms, pmids, homologs)
    print "Number of unfound proteins = ", protCount
    output_handle = open("data/mouse_proteins.pik", "w")
    cp.dump(prot_pmid_go, output_handle)
    output_handle.close()
    return prot_pmid_go


def GetGO(go_obo_file):
    go = {}
    for go_term in gop.parseGOOBO(go_obo_file):
        if not go_term.has_key('is_obsolete'):
            go[go_term["id"][0]] = (go_term["name"][0], go_term["def"][0]) 
    return go


def GetHumanProteinData(hu_prots, hu_go, go_obo_file):
    """ Build the human protein flat file from human fasta file and writes it to disk with the name
    human_proteins.pik
    
    Args:
        hu_prots: Human protein fasta file
        go_obo_file: The file that contains the whole gene ontology in obo format
        hu_go: Dictionary that maps human proteins to GO terms
    
    Returns:
        None
    """
    GO = GetGO(go_obo_file)
    prot_pmid_go = {}
    for seq_record in SeqIO.parse(hu_prots, "fasta"):
        protein = seq_record.id.split()[0]
        pmids = seq_record.description.split()[1].split(';')
        seq = seq_record.seq
        go_terms = []
        organism = "Human"
        if protein in hu_go:
            for term in hu_go[protein][1]:
                if term in GO:
                    go_terms.append((term, hu_go[protein][1][term].split(":")[0], GO[term][0], GO[term][1]))
                else:
                    go_terms.append((term, hu_go[protein][1][term].split(":")[0], "", ""))
        prot_pmid_go[protein] = PT.Protein(protein, organism, seq, go_terms, pmids, [])
    output_handle = open("data/human_proteins.pik", "w")
    cp.dump(prot_pmid_go, output_handle)
    output_handle.close()
    return prot_pmid_go


def GetOrganism(desc_names):
    """ Get the list of organisms (mouse, human or both) that a paper studies from the
     PubMed XML description file
    
    Args:
        desc_names: the list of mesh terms that describe a paper
    
    Returns:
        The names of the organism concatenated by ampersands
    """"
    
    names = []
    for name in desc_names:
        names.append(name[0])
        names.extend(name[1])
    if 'Humans' in names and 'Mice' in names:
        return "Human&Mouse"
    elif 'Humans' in names:
        return "Human"
    elif 'Mice' in names:
        return "Mouse"
    else:
        return "Other"

def GetSinglePaperData(pmid):
    """ Retrieves a paper data by its pmid from entrez which is a pubmed API to access publication data.
    The data to retrieve is: mesh headings, publication types and organisms
    
    Args:
        pmid: The pubmed id of the paper
    
    Returns:
        Tuple: (List of mesh headings, list of publication types, organisms)
    """
    
    pmidHeadings = []
    handle = None
    pubtypes = []
    try:
        handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
        records = Entrez.parse(handle)
        for pubmed_rec in records:
            headings = pubmed_rec['MedlineCitation']['MeshHeadingList']
            pubtypes = [str(ptype) for ptype in pubmed_rec['MedlineCitation']['Article']['PublicationTypeList']]
            for head in headings:
                pmidHeadings.append((str(head['DescriptorName']), [str(qualifier) for qualifier in head['QualifierName']]))
    except:
        return ([], [], "other")
    return (pmidHeadings, pubtypes, GetOrganism(pmidHeadings))
        

def GetPaperData(pmid_refs):
    """ Create a dictionary the maps each paper pubmed id to the data the describe it including:
    mesh headings, publication types and organisms.
    
    Args:
        pmid_refs: dictionary that maps each paper to its list of citations
    
    Returns:
        Dictionary: Dictionary that maps each paper to a tuple that has three elements: mesh headings,
        publication types and organisms
    """
    
    pmid_data = {}
    for pmid in pmid_refs:
        paperData = GetSinglePaperData(pmid)
        pmid_data[pmid] = PP.Paper(pmid, paperData[2], pmid_refs[pmid], paperData[0], paperData[1])
        for ref_pmid in pmid_refs[pmid]:
            if ref_pmid not in pmid_data and ref_pmid not in pmid_refs:
                paperData = GetSinglePaperData(ref_pmid)
                pmid_data[ref_pmid] = PP.Paper(ref_pmid, paperData[2], [], paperData[0], paperData[1])
    output_handle = open("data/papers.pik", "w")
    cp.dump(pmid_data, output_handle)
    output_handle.close()
    return pmid_data


def ParseMouseFasta(uniprot_mouse):
    """ Parse mouse protein fasta file and store data in dictionary on file under mouse_uniprot_proteins
    
    Args:
        uniprot_mouse: The name of the mouse protein fasta file
    
    Returns:
        None
    """
    
    mouse_prots = {}
    count = 0
    for seq_record in SeqIO.parse(uniprot_mouse, "fasta"):
        protein = seq_record.id.split('|')[1]
        seq = seq_record.seq
        mouse_prots[protein] = seq
        count += 1
    print "Count = ", count
    fHandler = open("mouse_uniprot_proteins", 'wb')
    cp.dump(mouse_prots, fHandler)
    fHandler.close()


def GetMeshData(mesh_tree):
    """Parses mesh tree bin file and stores the mesh terms along with their ids, parents, ancestors 
    and categories on file mesh_data_final.pik
    
    Args:
        mesh_tree: The name of the mesh bin flat file
    
    Returns:
        Dictionary: Dictionar that maps each mesh term to a tuple of its ids, parents, categories
        and ancestors
    """
    meshData = {}
    MD = {}
    with open(mesh_tree) as f:
        for line in f:
            line = line.strip()
            sc = line.index(';')
            i = len(line)-1
            ancestors = []
            while i > sc:
                if line[i] == '.':
                    ancestors.append(line[sc+1:i])
                i -= 1
            meshTerm = line[:sc]
            termID = line[sc+1:]
            category = MC.Mesh_Cats[line[sc+1]]
            parent = ""
            if len(ancestors) > 0:
                parent = ancestors[0]
            
            MD[termID] = meshTerm
            if meshTerm not in meshData:
                meshData[meshTerm] = MT.MeshTerm([termID], meshTerm, category, parent, ancestors)
            else:
                meshData[meshTerm].category.append(category)
                meshData[meshTerm].parent.append(parent)
                meshData[meshTerm].tid.append(termID)
                meshData[meshTerm].ancestors.append(ancestors)
                
    
    for term in meshData:
        ancList = []
        for ancs in meshData[term].ancestors:
            ancSet = Set()
            for anc in ancs:
                ancSet.add(MD[anc])
            ancList.append(list(ancSet))
        meshData[term].ancestors = ancList
    output_handle = open("../data/mesh_data_final.pik", "w")
    cp.dump(meshData, output_handle)
    output_handle.close()
    return meshData
            


if __name__ == "__main__":
    """ Main function that shows examples of function calls in this module:
    
    - Parse mouse protein fasta file and store data on disk:
    
    uniprot_mouse = "uniprot-mouse.fasta"
    ParseMouseFasta(uniprot_mouse)
    
    
    - Get the mouse protein data and store it in flat file:
    
    mouse_goa = "gene_association.goa_mouse"
    obo_file = "go.obo"
    uniprot_prots = cp.load(open("mouse_uniprot_proteins"))
    homologs = cp.load(open("data/mouse_homologs.pik"))
    data = GetMouseProteinData(mouse_goa, obo_file, uniprot_prots, homologs)
    
    
    - Get the human protein data and store it in flat file:
    
    human_prots = "human_proteins.fasta"
    human_go = cp.load(open("human_go.pik"))
    obo_file = "go.obo"
    homolog_prots = cp.load(open("mohuHits.pik"))
    data = GetHumanProteinData(human_prots, human_go, obo_file, homolog_prots)
    
    
    - Get the papers data:
    
    mouse_cits = cp.load(open("pmid_pmid_mouse"))
    data = GetPaperData(mouse_cits)
    
    - Parse the mesh tree bin file and store the structure on disk
    
    mesh_tree = "../data/mesh_tree.bin"
    GetMeshData(mesh_tree)
    
    """
