from Bio.Blast.Applications import NcbiblastpCommandline
from subprocess import call
import sys
import os

E_VALUE = 10**-8

def build_BLAST_db(fasta_file, blast_db):
    cmd= "makeblastdb -in " + fasta_file + " -dbtype 'prot' -out " + blast_db
    call(cmd, shell=True)
    
def make_BLAST_query(query_file, db_file, out_file):
    cline = NcbiblastpCommandline(query=query_file, db=db_file, outfmt=5, evalue=E_VALUE, out=out_file)
    print (cline)
    cline()

if __name__ == "__main__":
    '''
    fasta_file = "BLAST/human.fasta"
    blast_db = "BLAST/DB/humanDB"
    build_BLAST_db(fasta_file, blast_db)
    print "Finish!!"
    '''
    '''
    query_file = "BLAST/mouse.fasta"
    db_file = "BLAST/DB/humanDB"
    out_file = "results.xml"
    make_BLAST_query(query_file, db_file, out_file)
    print "finish"
    '''
    query_file = sys.argv[1]
    db_file = sys.argv[2]
    out_file = sys.argv[3]
    make_BLAST_query(query_file, db_file, out_file)
    print "finish"
    
    
    '''
    query_file = sys.argv[1]
    db_file = sys.argv[2]
    out_file = sys.argv[3]
    #fasta_file = os.path.join(CURR_PATH, "PCPAnnotationAssessment/Blasting_Data/CD-HIT/human_0.9/nr_query")
    #db_file = os.path.join(CURR_PATH, "PCPAnnotationAssessment/Blasting_Data/human_cit_db/humancitblast.db")
    #build_BLAST_db(fasta_file, db_file)
    #out_file = "dummyresult.xml"
    make_BLAST_query(query_file, db_file, out_file)
    print "Finish!!"
    '''
    