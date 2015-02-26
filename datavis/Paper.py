class Paper:
    def __init__(self, pmid, organism, citations, meshHeadings, pubTypes):
        self.pmid = pmid
        self.organism = organism
        self.citation = citations
        self.meshHeadings = meshHeadings
        self.pubTypes = pubTypes
    
    