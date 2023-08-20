#!/usr/local/bin/python3

##################################################################################################
#Script Name    : ICA2_script.py                                                                                          
#Description    : Identify a family of protein sequences from a taxonomic group. Run with "ICA2_function.py" file.                                                                          
#Date           : 2022_11_20                                                                                       
#Author         : B222908                                          
#help document  : B222908-2022.ICA2.manuals.pdf
# version       : 3.1                                     
##################################################################################################

import os,re
import ICA2_function


# run pro_tx_seq() function to let user input protein name and taxonomic group id or name
# check if the input invalid and check the number of sequences retrieved
# run a loop until the input and retrieve is resonable
dataset=ICA2_function.pro_tx_seq()  
# return a list [protein,taxonomy,protaxpar_label]

# Make a dataframe containing contents from fasta file and calculate the sequence length
df=ICA2_function.create_df(dataset[2]) 
# columns: ['protein_accession', 'protein_name', 'genus', 'species',
#           'full_species_name', 'length', 'header', 'sequence', 'raw_seq']

# Check how many species contained in the dataset, output the number of proteins, genus and species
# Output the calculation result of how many proteins per species and how many proteins per genus
ICA2_function.count_species(dataset[2])

# ask user to continue or not
condition=input("Do you want to continue for the alignment? yes/no:\n")
if re.search('n',condition.lower()):
    print("Bye!\n")
    quit()

while True:
    # checks the distribution of sequence length and generate a subset group
    # if the difference of sequences length is too big: (maxLen-minLen)/maxLen > 0.2,
    # ask user to choose an appropriate protein subset for alignment based on the range of sequence length
    sub_label=ICA2_function.seq_len_subset(dataset)

    # determine, and plot, the level of protein sequence conservation
    # basic information about sequences in an input multiple sequence alignment
    ICA2_function.conservation(dataset,sub_label)

    condition=input("How is the conservation of this group? Do you want to pick another range of length? yes/no:\n")
    if re.search('y',condition.lower()):
        continue
    else:
        break


# ask user to continue or not
condition=input("Do you want to continue finding motifs? yes/no:\n")
if re.search('n',condition.lower()):
    print("Bye!\n")
    quit()

# output fasta file for each sequence for motifs scanning
ICA2_function.seq_out(sub_label)

# scan with motifs from PROSITE database , names of motifs
ICA2_function.find_motifs(sub_label)

# remove unwanted files
os.system("rm -f *"+sub_label+"_len_single_seq.fa *"+sub_label+"_len_single.patmatmotifs "+sub_label+"_len.motifs.count count.txt")

# end
print('Program finished. Bye!')

