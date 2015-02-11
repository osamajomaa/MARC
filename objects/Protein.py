class Protein:
    def __init__(self, pid, organism, seq, go_terms, papers, homologs):
        self.pid = pid
        self.organism = organism
        self.sequence = seq
        self.go_terms = go_terms
        self.papers = papers
        self.homologs = homologs
    
    