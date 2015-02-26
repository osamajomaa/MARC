import cPickle as cp

def getCitedBy(paperList):
    paper_cits = {}
    for paper in paperList:
        paper_cits[paper] = []
    for paper in paperList:
        for ref in paperList[paper].citation:
            if ref in paper_cits:
                paper_cits[ref].append(paper)
    output_handle = open("../data/paper_cits.pik", "w")
    cp.dump(paper_cits, output_handle)
    output_handle.close()


if __name__ == "__main__":
    papers = cp.load(open("../data/papers.pik"))
    getCitedBy(papers)