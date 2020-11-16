#"(1) Obtain GeneIDs.py"
#Updated 7/20/20

#This is the first of two programs used in the bioinformatic analysis of the Nedd4 family of
#E3 ligases.

#The purpose of this program is to take in a Tab 3.0 file downloaded from BioGrid
#containing the interactome of a given protein, and write a text file with all the
#unique Gene IDs which can then be uploaded to UniProt (https://www.uniprot.org/uploadlists/)
#in order to obtain a FASTA file that will be used in the second program PxY_Finder


import csv

#Creates the set which all the GeneIDs will be added to
gene_IDs = set()


#Insert the path and file name of the downloaded BioGrid Tab 3.0 file here
BioGridPath = '/Users/michaelpupi/Desktop/Summer Research Fellowship/PY motif project/ITCH PY Motif Analysis/ITCH biogrid file.txt'

#Insert the path and file name of the output GeneIds .txt file here
GeneIDsOutput = '/Users/michaelpupi/Desktop/Summer Research Fellowship/PY motif project/ITCH PY Motif Analysis/ITCH GeneIDs redo.txt'


#Opens the .txt file downloaded from BioGrid containing the protein interactome as a csv
    
with open(BioGridPath, newline = '') as f:
    reader = csv.reader(f, delimiter = '\t')

    #This loop iterates through each line in the csv (skipping the first)
    #and adds the geneIds of the two proteins that are interacting to the set
    # (duplicates are not added because it is a set)
    for iteration in reader:
        if iteration[1] != "Entrez Gene Interactor A":
            gene_IDs.add(iteration[1])
            gene_IDs.add(iteration[2])

#Writes a new .txt file containing all the unique GeneIDs
    #Insert new path and file name for the GeneIDs .txt file here
with open(GeneIDsOutput,"w") as file:
    for ID in gene_IDs:
        file.write(ID + "\n")
