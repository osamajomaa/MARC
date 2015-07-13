## MARC - Mouse humAn Research Classifier

### The Problem
In a study published in PNAS biomedical journal in 2010, a group of researchers studied the correlation in transcriptional responses between mouse and human. The researchers studied three different diseaseses: Burn, Trauma and Endotoxemia and found out that there was poor to no correlation in the genomic responses between mouse and human, let alon within mice themselves. This shows that mouse is still not a perfect lab animal! Therefore it's important to identify research areas suitable to use mosue as a model organism to study human diseases. That's where MARC comes!

### What Is MARC Anyway?
MARC is a project that aims to classify the literature that has been studied on mouse as a model organism for human disease research. It does this through citation networks and a controlled vocabulary which is Medical Subject Heading.

**Citation Network**: It's a directed graph that connects citing papers to cited papers. It models the citation relationship between the citing mouse papers that cite human papers. These papers are then clustered into different classes of research using a unified medical vocabulary. 

**Medical Subject Heading**: It's a unified medical vocabulary that index papers in PubMed. It's organized in a DAG structure and contains 16 top level categories and over 27K terms altogether. MARC used MeSH to group the mouse and human papers in the citation networks into different classes of research

### How Does MARC Work?
MARC first collects the mouse and human paper data from Uniprot-GOA database. It parses the gene association files for mouse and human to extract the protein the paper data. From the PMIDs for the mouse and human papers, it queries another database which is Scopus to get the citation list for each mouse paper. It then stores the paper and citation data in pickled files on disk. To get the MeSH vocabulary data, MARC parses the bin file that contains the entire vocabulary and store it in a special structure on disk. Finally, MARC reads the paper, protein and mesh data from the files on disk and dump the data in MARC Database.

### What Is MARC Composed Of?
MARC is comprised of:<br>
**1. MARC Database:** A MongoDB database that contains data about the mouse and human papers in PubMed, the proteins that they study and the Medical Subject Heading terms that they're annotated to.<br>
**2. MARC Scripts:** Python script files that parse data files and extract data, connect to the database and execute queries, build BLAST database and query it and build graph files.

***
## MARC How-To Guide
In this tutorial I'm going to explain how MARC DB is organized and give example of how you can connect to the database and execute queries. I'm also going to explain the MARC script files, how to use them and what each module is supposed to do.

### MARC DB
##### Description
MARCDB is a NoSQL database that stores data about:
1. Mouse and human papers that describe proteins in the Uniprot-GOA database
2. Mouse and human proteins in the Uniprot-GOA and GenBank database
3. MeSH vocabulary that is used to annotate and describe papers in the MedLine database

#####Structure
In this database we have three collections: Paper, Protein and MeshTerm. Here is a description for each collection:

**Paper**
+ **PMID**: The pubmed id for the paper
+ **PubTypes**: The list of publication types for this paper
+ **MeshHeadings**: The list of the MeSH headings that the paper is annotated to
+ **Citations**: The list of citations of the paper
+ **Organism**: The organism that the paper studies
+ **CitedBy**: The list of papers cited by this paper

**Protein**
+ **PID**: Protein ID,
+ **Organism**: The organism of the protein
+ **Sequence**: The amino acid sequence of the protein
+ **GO Terms**: The Gene Ontology terms that the proteins is annotated to
+ **Papers**: The list of papers that the protein is studied in
+ **Homologs**: The homologs of the protein acquired by BLAST

**MeshTerm**
+ **Name**: The name of the mesh term
+ **TID**: The list of the term IDs
+ **Categories**: The list of categories that the term falls under in the 16 MeSH top-level categories
+ **Parents**: The list direct MeSH parents of the term
+ **AncestorsList**: The list of the term MeSH ancestors 

<br>
##### Restoring a MongodB backup
First you need to restore the binary backup in the MarcDB tar folder:<br>
1. Untar the MarcDB.tar.gz folder<br>
2. Execute this command in the terminal on a computer with MongoDB installed on it:
```sh
$ mongorestore <path to backup directory>
```

##### Sample Queries

1. View all databases in your database server:

  ```sh
    $ show databases
  ```
2. Select a database to execute queries on:

  ```sh
    $ use <db name>
  ```
  ***Note*** All the following commands can only be executed after you have selected a database (command #2).

3. Show all collections(tables) in a database:

  ```sh
    $ show collections
  ```
4. Show all documents(rows) in a collection:

  ```sh
    $ db.<collection name>.find()
  ```
5. Show only one documents from a collection:

  ```sh
    $ db.<collectin name>.findOne()
  ```
  Queries to a Mongo database are in JSON format; {Key:Value}

6. Find Exact Matching. Find the paper with the PMID=9642104:

  ```sh
    $ db.Paper.find({PMID:"9642104"})
  ```
7. Find in array. Find all papers that has the publication type is journal article where the field PubTypes is of type array:

  ```sh
    $ db.Paper.find({PubTypes:"Journal Article"})
  ```
8. Query field of an embedded document. Find all papers that has a "Aggrecans" as a mesh descriptor in the mesh headings array:

  ```sh
    $ db.Paper.find({"MeshHeadings.Descriptor":"Aggrecans"})
  ```
  ***Notice*** the use of quotations around the name of the array and subfield when querying an embedded document

9. Get all ancestors of the "Genetic Code" mesh term in the mesh vocabulary:

  ```sh
    $ db.MeshTerm.find({Name:"Genetic Code"}, {AncestorsList:1})
  ```
  ***Notice*** the 1 is used to limit the fields returned to AncestorsList for all matching documents

10. Get the names of all descendents of the "Genetic Phenomena" mesh term in the mesh vocabulary:

  ```sh
    $ db.MeshTerm.find({"AncestorsList.Ancestors":"Genetic Phenomena"}, {Name:1})
  ```

11. Get the parent(s) codes of the "Genetic Code" mesh term in the mesh vocabulary:

  ```sh
    $ db.MeshTerm.find({Name:"Genetic Code"}, {Parents:1})
  ```

### MARC Scripts
I will list the script modules that I wrote to build MARC in order and what each does.

#####1. database
This module collects the paper, protein and mesh data from different databases, store the data in files on disk and then dump it in MARC DB. It also contain the functionality necessary to connect to the database and execute queries against it. It comprises **MoPapersAndCits.py** which parses the mouse gene association and files to get the mouse proteins and human fasta files to get the human proteins. It also connects to PubMed Entrez API to get the MeSH annotations for each mouse and human paper and build a tree-like structure in a python dictionary and store it on file. **FileDumper.py** reads the paper, protein and mesh data files on disk, and write this data to MARC DB. **DB.py** contains the functions necessary to connect to the database, insert data to each of the three collections (Paper, Protein and MeshTerm) and exectue queries such as searching for a paper by its PMID or getting the list of ancestors for a particulary MeSH term. 

#####2. blast
This module creates a BLAST database, query it and parse the BLAST hit result files. **do_blast.py** creates a BLAST human protein database and uses a mouse fasta protein file to query the database to get the homologs for each mouse protein. **BLAST_Query_Builder.py** builds the mouse protein fasta query file. **blast_parser.py** parses blast XML hit files where it extracts the homologs for the query mouse proteins and store them in a dictionary that maps each mouse protein id to its list of homlogs. **FASTA_Splitter.py** This is an optional file useful to break down fasta files into chunks that can be used to query a BLAST database in a parallel manner using a qsub commands in PBS.

#####3. objects
This modules contains the classes that represent document in each Protein, Paper and MeshTerm collection. These classes used to offlaod and onload data from and to MARC DB.

#####4. datavis
**graphdata.py** includes functions that create graph files in GEXF format. **getCitedBy.py** is just a helper module that include a function which maps each paper to its citations.
