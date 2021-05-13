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
BioGridPath = 'your-biogrid-file-path.txt'

#Insert the path and file name of the output GeneIds .txt file here
GeneIDsOutput = 'your-output-file-name.txt'


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

___________


#"(2) PxYFinder.py"
#Updated 7/20/20

#This is the second of two programs used in the bioinformatic analysis of the Nedd4 family of
#E3 ligases.


#The program begins taking in a list of FASTA sequences downloaded from UniProt as a
#single .fasta file.  Over the course of running this program, the following files will
#be produced:
    # (1) a CSV containing UniProt Entry Name, UniProtKB ID, and FASTA sequence for each interactor
    # (2) a CSV containing UniProt Entry Name, UniProtKB ID, and FASTA sequence for each interactor
        #in addition to 4 columns containing either 0 (no motif) or 1 (motif present) for
        #"PP", "PPY", "PPxY", and "LPxY"
    # (3) a .txt file containing a statistical breakdown summarizing motif presence in the interactome
    # (4) a .txt file containing a list of sequences 10 amino acids upstream and downstream the PPxY
        #and LPxY motifs
    # (5) a .txt file containing Uniprot accession numbers separated by motif for use with PANTHER

import csv

#Insert FASTA file path here
FastaFilePath = "your-fasta-file-path.fasta"

# (1) Insert file path / name for the first CSV to be created (UniProt Entry Name, UniProtKB ID, and FASTA sequence for each interactor)
csv_one_path = "your-output-path-1.csv"

# (2) Insert file path / name for the second CSV to be created (same as above with the addition of "PP", "PPY", "PPxY", and "LPxY")
csv_two_path  = "your-output-path-2.csv"

# (3) Insert file path / name for the .txt file with a statistical breakdown of motifs
summary_path = "your-output-path-summary.txt"

# (4) Insert file path / name for the .txt file containing sequences to be used for WebLogos
weblogo_path  = "your-output-path-weblogo.txt"

# (5) Insert file path / name for the .txt file containing Uniprot accession numbers for PANTHER
accessions_path = "your-output-path-accessions.txt"

#__________________________________________________________________________________________
#FASTA TO CSV
    #This part of the program uses the input FASTA file to create the first csv file (1)
#__________________________________________________________________________________________


with open(FastaFilePath) as f:
    names = []
    accessions = []
    sequences = []

    #This string temporarily holds parts of sequence as they are being put together
    sequence_holder = ""


    #Turns the FASTA file into a list of strings (making the lines easier to navigate back and forth through later)
    file = []
    for line in f:
        file.append(line.strip("\n"))

    #This integer is used to iterate through the lines and keep track of location in the following while loop
    x = 0

    #The while loop creates parallel lists for Entry Name, Acession, and Fasta sequence
    while x < len(file):

        #If a line starts with ">", identifying information is extracted
        if file[x].startswith('>'):
            index_one = 0
            index_two = 0
            #loops through line/string beginning at index of the first "|" (where the accession begins)
            for y in range(4,len(file[x])):
                if file[x][y] == '|':
                    index_one = y
                if file[x][y] == ' ':
                    index_two = y
                    break
            #adds the accession number between the two lines (e.g. |'accession'|)
            accessions.append(file[x][4:index_one])
            #adds the name between the second '|' and the first space ' '
            names.append(file[x][index_one+1:index_two])
            x = x + 1

        #this part of the code processes the actual amino acid sequences,
            #combining different lines into one string for a given protein
        elif file[x].startswith('M'):
            #loop continues until the next identifier is met (ie until the sequence ends)
            while (file[x][0]!= '>'):
                sequence_holder = sequence_holder + file[x]
                x = x + 1

                #if the end of the file is reached, the loop will have to break to avoid
                #an out of bounds exception
                if x>=len(file):
                    sequences.append(sequence_holder)
                    sequence_holder = ""
                    break

            #add the complete sequence to the list of sequences
            sequences.append(sequence_holder)
            #reset the sequence holder
            sequence_holder = ""


#writes parallel lists into csv file
with open(csv_one_path, "w", newline='') as new_file:
    writer = csv.writer(new_file)
    writer.writerow(["name","accession","sequence"])
    for x in range(0,len(accessions)):
        writer.writerow([names[x],accessions[x],sequences[x]])

#__________________________________________________________________________________________
#PY MOTIF ANALYSIS
        #This part of the program uses the CSV created in the first part of the program (1) to
        #look for motifs, adding new columns with either 0 or 1 indicating whether a given
        #motif is present, generating a new CSV (2)
#__________________________________________________________________________________________


#Opens csv created in the previous part of the program
with open(csv_one_path,"r") as csv_file:
    reader = csv.reader(csv_file)
    rows = list(reader)

#each protein is now part of a list of lists, where:
    #the first element is the name
    #the second element is the accession #
    #the third element is the FASTA sequence

#creates the list which will eventually be used to make a new CSV (2) with motif information
data = []

#create parallel lists for different motifs
PP_list = []
PPY_list = []
PPxY_list = []
LPxY_list = []

#loops through each line in the first CSV (1)
for row in rows:
    if row[0] == "name":
        continue

    data.append(row)
    
    #both forward FASTA sequence collected so each can be iterated through
    sequence = row[2]
    

    #using basic "in" command, checks if the forward or reverse sequences contain the "PP" and "PPY" motifs
    if( ("PP" in sequence):
        PP_list.append(1)
    else:
        PP_list.append(0)

    if( ("PPY" in sequence):
        PPY_list.append(1)
    else:
        PPY_list.append(0)


    #for the PPxY and LPxY motifs, since there is variability in the motif
        #this manually goes through the string to see if the motif is present
        #stop at 3 before the very end to prevent and out of bounds error (and because the sequence has ended)
    PPxY = 0
    for x in range(0,len(sequence)-3):
        if sequence[x] == "P" and sequence[x+1] == "P" and sequence[x+3] =="Y":
            PPxY = 1

        
    PPxY_list.append(PPxY)
    
    LPxY = 0
    for x in range(0,len(sequence)-3):
        if sequence[x] == "L" and sequence[x+1] == "P" and sequence[x+3] =="Y":
            LPxY = 1
        
    LPxY_list.append(LPxY)

#writes the list of lists (data) into a new CSV
with open(csv_two_path,"w",newline = '') as new_file:
    writer = csv.writer(new_file)
    writer.writerow(["name","accession","sequence","PP","PPY","PPxY","LPxY"])

    for x in range(0,len(data)):
        name = data[x][0]
        accession = data[x][1]
        FASTA = data[x][2]
        writer.writerow([name,accession,FASTA,PP_list[x],PPY_list[x],PPxY_list[x],LPxY_list[x]])
    
#__________________________________________________________________________________________
#SUMMARY OF MOTIFS (PERCENTAGE / BREAKDOWN)
        #This part of the program analyses the motif data generated in the previous part,
        #outputting the percentage of the interactome that has each motif into a .txt file (3)
#__________________________________________________________________________________________


#Open CSV file with genes, accessions, fasta sequences, and motifs created in the previous section (2)
with open(csv_two_path,"r") as csv_file:
    reader = csv.reader(csv_file)
    rows = list(reader)

#the information that will be extracted from the file
names = []
accessions = []

#these lists will contain either 0 or 1 corresponding to whether or not the motif is present
PP = []
PPY = []
PPxY = []
LPxY = []

#these lists contain the names of all proteins with the corresponding motif
PP_hits = []
PPY_hits = []
PPxY_hits = []
LPxY_hits = []


#extracts useful information from the file into parallel lists
for row in rows:
    if row[0] == 'name':
        continue
    names.append(row[0])
    accessions.append(row[1])
    PP.append(row[3])
    PPY.append(row[4])
    PPxY.append(row[5])
    LPxY.append(row[6])

#counts the number of each motif
for x in range(0,len(names)):
    if PP[x] == "1":
        PP_hits.append(names[x])
    if PPY[x] == "1":
        PPY_hits.append(names[x])
    if PPxY[x] == "1":
        PPxY_hits.append(names[x])
    if LPxY[x] == "1":
        LPxY_hits.append(names[x])

number_of_PP = len(PP_hits)
number_of_PPY = len(PPY_hits)
number_of_PPxY = len(PPxY_hits)
number_of_LPxY = len(LPxY_hits)


#Writes file, calculating motif proportions
with open(summary_path,"w") as output:
    output.write("PP hits: %s\n"%number_of_PP)
    output.write("PPY hits: %s\n"%number_of_PPY)
    output.write("PPxY hits: %s\n"%number_of_PPxY)
    output.write("LPxY hits: %s\n"%number_of_LPxY)

    percentage_of_PP = number_of_PP/len(names)
    percentage_of_PPY = number_of_PPY/len(names)
    percentage_of_PPxY = number_of_PPxY/len(names)
    percentage_of_LPxY = number_of_LPxY/len(names)

    output.write("\n")

    output.write("Fraction of PP: %s\n"%percentage_of_PP)
    output.write("Fraction of PPY: %s\n"%percentage_of_PPY)
    output.write("Fraction of PPxY: %s\n"%percentage_of_PPxY)
    output.write("Fraction of LPxY: %s\n"%percentage_of_LPxY)

    with_motif = list(set(PPxY_hits) | set(LPxY_hits))
    no_motifs = list(set(names)-set(with_motif))
    
    output.write("\n")
    output.write("All proteins with either PPxY or LPxY:\n")
    
    for protein in with_motif:
        output.write("%s\n"%protein)
    
    output.write("\n")
    output.write("All proteins with neither PPxY or LPxY:\n")

    for protein in no_motifs:
        output.write("%s\n"%protein)


#__________________________________________________________________________________________
#OBTAIN SEQUENCES FOR WEBLOGO
    #This part of the program collects all sequences containing the PPxY and LPxY motifs,
    #writing those with at least 10 AAs upstream and downstream into a .txt file (4)
#__________________________________________________________________________________________


#Function that tests for the presence of the PPxY motif beginning at a specific index
    #True if present, false if not
def is_PPxY(start,sequence):
    if start + 3 >= len(sequence):
        return False
    if sequence[start] == "P" and sequence[start+1] == "P" and sequence[start+3] =="Y":
            return True      
    return False

#Function that tests for the presence of the LPxY motif beginning at a specific index
    #True if present, false if not
def is_LPxY(start,sequence):
    if start + 3 >= len(sequence):
        return False
    if sequence[start] == "L" and sequence[start+1] == "P" and sequence[start+3] =="Y":
            return True      
    return False

#These lists will hold fragments of the sequences with either PPxY or LPxY
PPxYs = []
LPxYs = []


with open(csv_one_path,"r") as csv_file:
    reader = csv.reader(csv_file)
    rows = list(reader)

#This for loop collects sequences containing PPxY
for row in rows:
    if row[0] == "name":
        continue 

    #Gets FASTA sequence, forward and reverse to iterate over in both directions
    sequence = row[2]
    rev_seq = sequence[::-1]

    #Loop runs through every index in each sequence, checking for motif
    for x in range(0,len(sequence)):
        if is_PPxY(x,sequence) == False:
            continue
        
        #Out of bounds errors avoided by checking where in the sequence the iterator is 
        if x >= 10 and x+14<len(sequence):
            seq = sequence[x-10:x+14]
        elif x < 10:
            seq = sequence[0:x+14]
        elif x+14 >= len(sequence):
            seq = sequence[x-10:len(sequence)]

        PPxYs.append(seq)

    #Exact same loop as above, but runs for sequence in reverse
    for x in range(0,len(sequence)):
        if is_PPxY(x,rev_seq) == False:
            continue

        if x >= 10 and x+14<len(sequence):
            seq = rev_seq[x-10:x+14]
        elif x < 10:
            seq = rev_seq[0:x+14]
        elif x+14 >= len(sequence):
            seq = rev_seq[x-10:len(sequence)]

        PPxYs.append(seq)
        


#This for loop collects sequences containing LPxY (identical function as above for PPxY)
for row in rows:
    if row[0] == "name":
        continue 

    #Gets FASTA sequence, forward and reverse to iterate over in both directions
    sequence = row[2]
    rev_seq = sequence[::-1]

    #Loop runs through every index in each sequence, checking for motif
    for x in range(0,len(sequence)):
        
        if is_LPxY(x,sequence) == False:
            continue
        
        #Out of bounds errors avoided by checking where in the sequence the iterator is 
        if x >= 10 and x+14<len(sequence):
            seq = sequence[x-10:x+14]
        elif x < 10:
            seq = sequence[0:x+14]
        elif x+14 >= len(sequence):
            seq = sequence[x-10:len(sequence)]

        LPxYs.append(seq)

    #Exact same loop as above, but runs for sequence in reverse
    for x in range(0,len(sequence)):
        if is_LPxY(x,rev_seq) == False:
            continue

        if x >= 10 and x+14<len(sequence):
            seq = rev_seq[x-10:x+14]
        elif x < 10:
            seq = rev_seq[0:x+14]
        elif x+14 >= len(sequence):
            seq = rev_seq[x-10:len(sequence)]

        LPxYs.append(seq)

    
#Writes obtained sequences of proper length (full 10 AA upstream and downstream) into a file
with open(weblogo_path,"w") as output:

    output.write("PPxY sequences")
    output.write("\n")

    skipped_ppxy_count = 0
    for seq in PPxYs:

        #Checks for proper length
        if len(seq)==24:
            output.write(seq)
            output.write("\n")
        else:
            skipped_ppxy_count = skipped_ppxy_count + 1
    output.write("\n")

    skipped_lpxy_count = 0 
    output.write("LPxY sequences")
    output.write("\n")
    for seq in LPxYs:
        #Checks for proper length
        if len(seq)==24:
            output.write(seq)
            output.write("\n")
        else:
            skipped_lpxy_count = skipped_lpxy_count + 1

    output.write("\n")
    output.write("NOTE: Some motifs not included because there weren't 10 AAs upstream or downstream.")
    output.write("\n")
    output.write("\n")
    output.write("%s PPxY sequences not included" % skipped_ppxy_count)
    output.write(" of %s total PPxY sequences." % len(PPxYs))
    output.write("\n")
    output.write("\n")
    output.write("%s LPxY sequences not included" % skipped_lpxy_count)
    output.write(" of %s total LPxY sequences." % len(LPxYs))


#__________________________________________________________________________________________
#OBTAIN ACCESSION NUMBERS FOR PANTHER
    #This section generates a text file with all accession numbers separated by whether a
    #given motif is present or not.
#__________________________________________________________________________________________


with open(csv_two_path,"r") as csv_file:
    reader = csv.reader(csv_file)
    rows = list(reader)

#These lists contain the accessions of proteins with the corresponding motif
PPxY_hits = []
LPxY_hits = []
noPPxY = []
noLPxY = []

#Populates lists with the correct accessions
for row in rows:
    if row[5] == "1":
        PPxY_hits.append(row[1])
    else:
        if row[1]=="accession":
            continue
        noPPxY.append(row[1])
    if row[6] == "1":
        LPxY_hits.append(row[1])
    else:
        if row[1]=="accession":
            continue
        noLPxY.append(row[1])

PPxY_or_LPxY = list(set(PPxY_hits) | set(LPxY_hits) )

neither_PPxY_nor_LPxY = list(set(noPPxY) & set(noLPxY))


with open(accessions_path,"w") as output:
    output.write("PPxY Accessions")
    output.write("\n")
    for num in PPxY_hits:
        output.write(num)
        output.write("\n")
    output.write("\n")

    output.write("LPxY Accessions")
    output.write("\n")
    for num in LPxY_hits:
        output.write(num)
        output.write("\n")
    output.write("\n")

    output.write("PPxY or LPxY Accessions")
    output.write("\n")
    for num in PPxY_or_LPxY:
        output.write(num)
        output.write("\n")
    
    output.write("\n")
    output.write("No PPxY")
    output.write("\n")
    for num in noPPxY:
        output.write(num)
        output.write("\n")
    
    output.write("\n")
    output.write("No LPxY")
    output.write("\n")
    for num in noLPxY:
        output.write(num)
        output.write("\n")
    
    output.write("\n")
    output.write("No PPxY or LPxY")
    output.write("\n")
    for num in neither_PPxY_nor_LPxY:
        output.write(num)
        output.write("\n")

