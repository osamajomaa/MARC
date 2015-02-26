class MeshTerm:
    def __init__(self, tid, termName, category, parent, ancestors):
        self.tid = tid
        self.name = termName
        self.category = [category]
        self.parent = [parent]
        self.ancestors = [ancestors]
    
    