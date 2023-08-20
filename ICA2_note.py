#!/usr/local/bin/python3

##################################################################################################
#Script Name	: ICA2_script.py                                                                                          
#Description	: Identify a family of protein sequences from a taxonomic group. Run with "ICA2_function.py" file.                                                                          
#Date           : 2022_11_20                                                                                       
#Author       	: B222908                                          
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



esearch -db taxonomy -query "mammals" | efetch -format xml | grep "<TaxId>" | head -1 > tax


prettyplot *align.msf -graph pdf -boxcol -consensus -gtitle "protein for taxonomy" -goutfile prettyplot

patmatmotifs -full -sequence *len_single_seq.fa -outfile *_patmatmotifs

cut -d ' ' -f 3 
os.system("cut -d ' ' -f 3 *count | sort > "+taxonomy+".motifs.txt")

elink -target protein | efilter -organism mammal -source refseq | efetch -format fasta

esearch -db taxonomy -query "mammals" | efetch -format xml  > mammals.xml
esearch -db protein -query "pyruvate dehydrogenase[Protein Name] AND txid4890[Organism]"| grep "Count"


esearch -db protein -query "pyruvate dehydrogenase[Protein Name] AND txid4890[Organism]" | efetch -format acc | wc -l
esearch -db protein -query "ABC transporters[Protein Name] AND mammals[Organism] " | efetch -format xml > test.xmlseq
esearch.fcgi?db=<database>&term=<query>

esearch -db protein -query "pyruvate dehydrogenase AND txid4890[Organism]" | efetch -format fasta  > proteins.fa
clustalo -i txid4890.fa -o txid4890.align.msf --outfmt=msf --threads=100
infoalign txid4890.align.msf -out txid4890.infoalign
sort -k10 txid4890.infoalign > txid4890.infoalign.sort
plotcon -sformat msf txid4890.align.msf -winsize 25 -gsubtitle 'pyruvate dehydrogenase sequences in txid4890' -graph svg -goutfile txid4890

    while True:
        print("Do you want to see more rows or specific genus?")
        test=input("Type in a number for more rows, or type in the initial letter for specific genus:\n")
        test=test.strip().replace(' ','').lower()
        try:
            num=int(test)
            if num>0 and num<=
            print(num,'rows showing:')
            sub_df.head(num)
        except:
            test
plotcon -sformat msf txid40674.align.msf -winsize 25 -gsubtitle 'pyruvate dehydrogenase sequences in txid40674' -graph png -goutfile txid40674


patmatmotifs -full -sequence test.fa1 -outfile txid4890.patmatmotifs

# Create an ambiguous consensus sequence
consambig -sequence txid4890.align.msf -outseq txid4890.aligned.cons
dataset=["ABC transporters",'txid40674']

'esearch -db protein -query "'+protein+'[Protein Name] AND '+taxonomy+'[Organism] NOT PARTIAL'


    # write result to file
    with open(protaxpar_label+'_proteins_species.txt','w') as output:
        output.write('accession\tspecies\n')
        for key,value in dic.items():
            output.write('%s\t%s\n' % (key, value))


def find_motifs(taxonomy,label,df):
    motif='for fn in *.seq.fa \n do patmatmotifs -full -sequence ${fn} -outfile ${fn//.seq.fa/}.patmatmotifs & \n done \n wait'
    os.system(motif)
    # count the number of each motif


    with multiprocessing.Pool() as pool:
        pool.map(moti, proteins)


    # count the number of each motif
    proteins=list(df['pro_number'])
    for protein in proteins:
        protein=protein.replace('|','\|')
        os.system('ls '+protein+'*')
        query_motif='patmatmotifs -full -sequence '+protein+'.seq.fa -outfile '+protein+'.patmatmotifs'
        os.system(query_motif)



            if len(test)==1:
                for i in list(range(0,spe_df.shape[0],1)):
                    if re.search(test.upper(),spe_df.iloc[i]['genus'][0]):
                        ini_test.append(True)
                    else:
                        ini_test.append(False)
                print(spe_df[ini_test].to_string())   

clustalo -i pyruvate_dehydrogenase_ascomycete_fungi_no_par_whole_len.fa  --distmat-out=matrix_file --threads=100

        # add brackets to protein name if protein name contain two words and taxonomy name only one word
        if re.search('\w\s+\w',protein) and taxonomy==taxonomy.replace(' ', '') and not re.search("txid",taxonomy):
            protein='('+protein+')'

    # raw_seq=[]
        # raw_seq.append(full_seq[i][(devide+1):]) # full fasta format

try:
        os.system("display "+protaxpar_label+".seqLenHis.png")
        # sometimes it failed, restart terminal will fix it
    except:
        print("Failed to open the image, try open the file yourself: "+protaxpar_label+".seqLenHis.png\n")
    

        # check if count.txt file exsist, if not, create file
        try: 
            open("count.txt").read()
        except:
            with open("count.txt",'w') as c:
                c.write('')