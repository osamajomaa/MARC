from Bio.Blast import NCBIXML
import cPickle
import os

RESULTS_PATH = "BLAST/Results"

def get_blast_hits():

    files = os.listdir(RESULTS_PATH)
    homologs = {}
    count = 1
    for fname in files:
        result_handle = open(os.path.join(RESULTS_PATH, fname))
        records = NCBIXML.parse(result_handle)
        for rec in records:
            query = rec.query.split()[0]
            hits = []
            if rec.alignments:
                rank = 1
                for alignment in rec.alignments:
                    protein = alignment.title.split('|')[2].split()[1]
                    hits.append((protein, rank))
                    rank += 1
            print "Done with ", count
            homologs[query] = hits
    fHandler = open("mouse_homologs.pik", 'w')
    cPickle.dump(homologs, fHandler)
    fHandler.close()

if __name__ == "__main__":
    get_blast_hits()
    print "Finish!"