from Bio.UniProt import GOA
import go_obo_parser as gop
from Bio import Entrez
from Bio import SeqIO
import cPickle as cp
import Mesh_Categories as MC
import Protein
import Paper
import MeshTerm

Entrez.email = "jomaao@miamioh.edu"

def GetMouseProteinData(gaf_file, go_obo_file, uniprot_prots, mouse_homologs):
    """
    This method retrieves the papers, go terms, organism and sequence for each protein a gene association file.
    @param gaf_file: Uniprot_GOA association file in gaf format.
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
        prot_pmid_go[prot] = Protein.Protein(prot, organism, seq, go_terms, pmids, homologs)
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


def GetHumanProteinData(hu_prots, hu_go, go_obo_file, mohuHits):
    GO = GetGO(go_obo_file)
    prot_pmid_go = {}
    for seq_record in SeqIO.parse(hu_prots, "fasta"):
        protein = seq_record.id.split()[0]
        pmids = seq_record.description.split()[1].split(';')
        seq = seq_record.seq
        go_terms = []
        #homologs = []
        organism = "Human"
        '''
        if protein in mohuHits:
            homologs = mohuHits[protein]
        '''
        
        if protein in hu_go:
            for term in hu_go[protein][1]:
                if term in GO:
                    go_terms.append((term, hu_go[protein][1][term].split(":")[0], GO[term][0], GO[term][1]))
                else:
                    go_terms.append((term, hu_go[protein][1][term].split(":")[0], "", ""))
        prot_pmid_go[protein] = Protein.Protein(protein, organism, seq, go_terms, pmids, [])
    output_handle = open("data/human_proteins.pik", "w")
    cp.dump(prot_pmid_go, output_handle)
    output_handle.close()
    return prot_pmid_go


def GetOrganism(desc_names):
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
    pmid_data = {}
    for pmid in pmid_refs:
        paperData = GetSinglePaperData(pmid)
        pmid_data[pmid] = Paper.Paper(pmid, paperData[2], pmid_refs[pmid], paperData[0], paperData[1])
        for ref_pmid in pmid_refs[pmid]:
            if ref_pmid not in pmid_data and ref_pmid not in pmid_refs:
                paperData = GetSinglePaperData(ref_pmid)
                pmid_data[ref_pmid] = Paper.Paper(ref_pmid, paperData[2], [], paperData[0], paperData[1])
    output_handle = open("data/papers.pik", "w")
    cp.dump(pmid_data, output_handle)
    output_handle.close()
    return pmid_data


def ParseMouseFasta(uniprot_mouse):
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
    meshData = {}
    with open(mesh_tree) as f:
        for line in f:
            line = line.strip()
            meshTerm = line.split(';')[0]
            level = len(line.split(';')[1].split('.'))
            termID = line.split(';')[1].split('.')[level-1]
            category = MC.Mesh_Cats[line.split(';')[1][0]]
            ancestors = line.split(';')[1].split('.')[0:level-1]
            ancesCount = len(ancestors)
            parent = ""
            if ancesCount > 0:
                parent = ancestors[ancesCount-1]
            if meshTerm not in meshData:
                meshData[meshTerm] = MeshTerm.MeshTerm(termID, meshTerm, category, parent, ancestors)
            else:
                meshData[meshTerm].category.append(category)
                meshData[meshTerm].parent.append(parent)
                meshData[meshTerm].ancestors.append(ancestors)
    output_handle = open("data/mesh_data.pik", "w")
    cp.dump(meshData, output_handle)
    output_handle.close()
    return meshData
            


if __name__ == "__main__":
    
    '''
    mouse_goa = "gene_association.goa_mouse"
    obo_file = "go.obo"
    uniprot_prots = cp.load(open("mouse_uniprot_proteins"))
    homologs = cp.load(open("data/mouse_homologs.pik"))
    data = GetMouseProteinData(mouse_goa, obo_file, uniprot_prots, homologs)
    '''
    
    '''
    human_prots = "human_proteins.fasta"
    human_go = cp.load(open("human_go.pik"))
    obo_file = "go.obo"
    homolog_prots = cp.load(open("mohuHits.pik"))
    data = GetHumanProteinData(human_prots, human_go, obo_file, homolog_prots)
    '''
    
    '''
    mouse_cits = cp.load(open("pmid_pmid_mouse"))
    data = GetPaperData(mouse_cits)
    '''
    
    
    mesh_tree = "mesh_tree.bin"
    GetMeshData(mesh_tree)
    
    '''
    uniprot_mouse = "uniprot-mouse.fasta"
    ParseMouseFasta(uniprot_mouse)
    '''
