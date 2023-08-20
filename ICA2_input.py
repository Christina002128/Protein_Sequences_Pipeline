#############################################################################
#Script Name	: ICA2_function.py                                                                                          
#Description	: Contain all functions. Run with "ICA2_script.py" file.                                                                            
#Date           : 2022_11_20                                                                                       
#Author       	: B222908                                          
#help document  : B222908-2022.ICA2.manuals.pdf  
#version        : 3.1                                     
############################################################################

# import model
import os
import pandas as pd
import re
import matplotlib.pyplot as plt
import shutil

# user input protein name and taxinomic group
def pro_tx_seq():
    while True: # a loop till correct input
        protein=input("Please type in protein name:\n").strip()
        taxonomy=input("Please type in a taxonomic name or a taxonomy ID. (taxonomy ID for example: txid4890)\n").lower().strip()
        # test if one of the input empty
        if protein=='' or taxonomy=='':
            print("Input missing. try again.\n")
            continue
        # test if protein name is valid
        if re.search('\W',protein.replace(' ','')):
            print("Invalid protein ID. try again.\n")
            continue
        # test if taxonomy input is valid
        if re.search('txid', taxonomy): # if the input is a taxonomy ID, then
            # test if id is integer
            taxonomy=taxonomy.replace(' ','') # delete spaces
            id=taxonomy.replace("txid",'') # save id number
            try:
                int(id) 
            except:
                print("Invalid taxonomy ID, try again.\n")
                continue  # skip to the next loop, input again
            print('Query for "'+protein+'" protein and "'+taxonomy+'" group.\n')
        else: # if input is only integer
            try:
                int(taxonomy.replace(' ','')) 
                taxonomy='txid'+taxonomy.replace(' ','') # if integer, add 'txid' in the front
                print('Query for "'+protein+'" protein and "'+taxonomy+'" group.\n')
            except: # if not, then it's taxonomic name
                print('Generate query for "'+protein+'" protein and "'+taxonomy+'" taxonomic group.')
        #partial or not partial search
        par=input("Would you like to contain partial proteins or not? yes/no:\n")
        par=par.replace(' ','').lower()
        # count the number of sequences searched from NCBI
        if re.search('n',par): # search not partial proteins
            print("Searching no partial protein sequences from ncbi...")
            ecount='esearch -db protein -query "'+protein+'[protein] AND '+taxonomy+'[organism] NOT PARTIAL" | grep "Count" > count.txt'
            os.system(ecount)
            par_query='NOT PARTIAL'  # for esearch query
            par="no_par"  # for filename label
        else: # search all proteins
            print("Searching all protein sequences from ncbi...")
            ecount='esearch -db protein -query "'+protein+'[protein] AND '+taxonomy+'[organism]" | grep "Count" > count.txt'
            os.system(ecount)
            par_query='' # for esearch query
            par='all' # for filename label
        count=open("count.txt").read()
        count=int(count.replace('<Count>','').replace('</Count>\n','').strip()) 
        # check if the number of sequences too high and ask user if they want to download
        if count>1000: 
            print("Too many sequences("+str(count)+")!  Continue downloading? yes/no:")
            test=input()
            if re.search('n',test.lower()):
                test2=input("Do you want to retype a group of proteins? yes/no:\n")
                if re.search('y',test2.lower()):
                    continue
                else:
                    print("Bye!\n")
                    quit() # quit the script
        # check if the retrieve result is 0
        if count==0 :
            if count==0:
                print("0 result searched.")
                print("Maybe try typing in protein name in singular form instead of plural.\nMaybe try typing in NCBI taxonomy ID instead of taxonomic name.\n")
                continue  # skip to the next loop, input again
        if count==1:
            print("Only 1 result. Cannot do alignment. Try another protein group.\n")
            continue
        else:
            print(str(count)+" sequences downloading...")
            # create label for naming the file
            tax=taxonomy.replace(' ', '_')
            pro=protein.replace(' ', '_').replace('(','').replace(')', '')
            protaxpar_label=pro+'_'+tax+'_'+par
            #download sequences
            fastafetch='esearch -db protein -query "'+protein+'[protein] AND '+taxonomy+'[organism] '+par_query+'" | efetch -format fasta > '+protaxpar_label+'.fa'
            os.system(fastafetch)
            print("protein sequences is inside file: "+protaxpar_label+".fa")
            break
    return [protein,taxonomy,protaxpar_label] # return a list




# Make a dataframe containing contents from fasta file and calculate the sequence length
def create_df(protaxpar_label):
    fas=open(protaxpar_label+".fa").read().strip()
    # get all headers 
    fas_list=fas.split('\n')
    headers=[]
    for line in fas_list:
        if '>' in line:
            headers.append(line)
    # get headers content: accession, protein name, full species name, genus name, species name
    acc=[]
    pro_name=[]
    full_species=[]
    genus=[]
    species=[]
    for pro in headers:
        # for sequences with regular header:  >accession protein_name [species]
        if re.search('\[', pro): 
            search=re.search(r'>(\S*)\s(.*)\[(.*)\]',pro)
            acc.append(search.group(1).strip()) # accession
            pro_name.append(search.group(2).strip()) # protein name
            spec=search.group(3).strip() # [species]
            full_species.append(spec) # full name
            genus.append(spec.split(' ')[0]) # genus name
            sp_list=spec.split(' ')[1:] # species name might include more than one word
            species.append(' '.join(sp_list)) # species name
        # for sequences with strange header, genus and species name columns equal to 'NA'
        else: 
            search=re.search(r'>(\S*)\s(.*$)',pro)
            acc.append(search.group(1).strip())
            pro_name.append(search.group(2).strip())
            full_species.append('empty')
            genus.append('empty')
            species.append('empty')
    # get sequence and sequence length
    full_seq=fas.split('>')[1:] # delete first empty string
    seq=[]
    seqlen=[]
    for i in list(range(0,len(full_seq),1)): 
        # count sequence length for each sequence
        devide=full_seq[i].find("\n") # the first index of '\n' splits the header and sequence
        seq.append(full_seq[i][(devide+1):].replace('\n','')) # only sequence
        seqlen.append(len(seq[i])) # length of sequence
    # make data frame
    df={} 
    df=pd.DataFrame({'protein_accession':acc,"protein_name":pro_name,
    "genus":genus,"species":species,"full_species_name":full_species,
    "length":seqlen,"header":headers,"sequence":seq})
    # sort by genus
    df=df.sort_values('genus',ascending=True) 
    df.to_csv(protaxpar_label+"_summary.csv",sep='\t',index=False) # output file
    return df



# Check how many species contained in the dataset
def count_species(protaxpar_label):
    # load dataframe
    df=pd.read_csv(protaxpar_label+"_summary.csv", sep="\t")
    # count the number of proteins, genus and species. 
    # ('empty' value roughly treated as one species, doesn't affect much. Because the fasta header doesn't contain the species doesn't means it has no organism sources)
    protein_numbers=len(list(df['protein_accession']))
    genus_numbers=len(set(list(df['genus'])))
    species_numbers=len(set(list(df['full_species_name'])))
    empty_numbers=len(df[df['genus']=='empty'])
    # number of protein per species or genus
    sp_percent=round(protein_numbers/species_numbers,2)
    ge_percent=round(protein_numbers/genus_numbers,2)
    # display result
    print("\nResult of retrieved data:")
    print("1. The result contains "+str(protein_numbers)+" sequences, "+str(genus_numbers)+" kinds of genus and "+str(species_numbers)+" kinds of species.")
    print("2. About "+str(sp_percent)+" proteins per species, "+str(ge_percent)+" proteins per genus.")
    # write species summary to csv file
    spe_df=df[['protein_accession',"protein_name",'genus','species']].reset_index(drop=True)
    spe_df.to_csv(protaxpar_label+"_summary_species.csv",sep='\t',index=False)
    print('3. Result inside "'+protaxpar_label+'_summary_species.csv" file.')
    print("(",empty_numbers,"sequences are without species labels, registered as 'empty')")
    # print out first 10 rows
    print("Here are the first 10 rows of the result sorted by genus name(sequence headers without specifying species are signed as 'empty').")
    print(spe_df.head(10))
    # user to choose what genus or number of rows to see
    while True:
        print("\nDo you want to see more rows or specific genus?")
        print("Type in an integer for more rows(positive number for head rows, negative number for tail rows), for example: 8 or -5")
        print("Or type in the initial letter or full name of a specific genus, for example: A or Acinonyx")
        test=input()
        test=test.strip().replace(' ','').lower()
        # input empty, break the loop
        if test=='':
            break
        # input is a number
        try:
            num=int(test)
        except:
            num=0
        if num>0: # heads of rows
            print(num,'head rows showing:')
            print(spe_df.head(num))
        elif num<0: # tails of rows
            num=abs(num)
            print(num,'tail rows showing:')
            print(spe_df.tail(num))
        # input is an initial letter or a word
        elif re.search('^\w*$',test):
            find_genus='^'+test # start with the word
            ini_test=[]
            # find the rows with genus name
            for i in list(range(0,spe_df.shape[0],1)):
                if re.search(find_genus,spe_df.iloc[i]['genus'].lower()):
                    ini_test.append(True)
                else:
                    ini_test.append(False)
            print(spe_df[ini_test].to_string())
        # invalid
        else:
            print("Invalid input, try again.\n")
            continue
        # user choose to stop the loop
        test2=input("Want to choose again? yes/no:\n")
        if re.search(r'^n',test2.replace(' ','').lower()):
            break




# choose a appropriate protein subset for alignment by considering the sequence length
def seq_len_subset(dataset):
    taxonomy=dataset[1]
    protaxpar_label=dataset[2]
    # load dataframe
    df=pd.read_csv(protaxpar_label+"_summary.csv", sep="\t")
    print("\nSummary of sequences length:")
    print(df.describe()) # print out calculations of length
    # display histogram of sequence length
    print("\nGenerating histogram of sequences length. (If it fails to generate plot, please restart a terminal window and run the program again!)")
    y=list(df['length'])
    n, bins, patches = plt.hist(y, 50, color='c')
    plt.grid(True)
    plt.xlabel('Sequence length')
    plt.ylabel('Number of proteins')
    plt.title('histogram of sequences length')
    plt.savefig(protaxpar_label+".seqLenHis.png")
    print("Generate plot succeeded. The plot is in "+protaxpar_label+".seqLenHis.png file.")
    print("Showing histogram.(If it fails to display the plot, try open the file yourself.)")
    print("Close the plot window to continue.")
    plt.show() # sometimes fails, need to reopen the terminal window
    # test if the difference of sequences length is too big: (max-min)/max > 0.2
    percentlen=(df.max()['length']-df.min()['length'])/df.max()['length']
    if percentlen > 0.2:
        test1=input("The difference of length of the sequences is quite big, would you choose a subset based on length range? yes/no:\n")
        if re.search('y',test1.lower()):
            # choose a subset
            test2=input("Would you like to choose sequences with length within range [median - sd, median + sd]? yes/no:\n")
            if re.search('y',test2.lower()):
                mini=int(df['length'].median()-df['length'].std())
                maxi=int(df['length'].median()+df['length'].std())
            else:
                # user input a range of length, test the validity
                while True:
                    print("Please choose a range of length.")
                    mini=input('minimum:\n').replace(' ', '')
                    maxi=input("maximum:\n").replace(' ', '')
                    try:
                        mini=int(mini)
                        maxi=int(maxi)
                    except:
                        print("invalid number. Please type in an integer number. Try again.\n")
                        continue
                    if(maxi>mini):
                        break
                    else:
                        print("minimum number should be less than maximum number. Try again.\n")
                        continue
            # subgroup in length range
            sub_df=df[(df['length']>=mini) & (df['length']<=maxi)].reset_index(drop=True)
            sub_label=protaxpar_label+'_'+str(mini)+'_'+str(maxi)
            # subgroup dataframe
            sub_df.to_csv(sub_label+"_subgroup.csv",sep='\t',index=False)
            # generate subgroup sequences fasta file
            try: # if the file exist, delete it
                os.remove(sub_label+'_len.fa')
                print("Re-generating subset of sequences with length between "+str(mini)+" and "+str(maxi)+".")
            except:
                print("Generating subset of sequences with length between "+str(mini)+" and "+str(maxi)+".")
            with open(sub_label+'_len.fa','w') as out:
                for i in list(range(0,len(sub_df),1)):
                    out.write(sub_df['header'][i]+"\n")
                    out.write(sub_df['sequence'][i]+"\n")
            print(len(sub_df),"sequences are selected.")
        else:
            sub_label=protaxpar_label+'_whole'
            print("Align for all the proteins.\n")
            shutil.copyfile(protaxpar_label+'.fa',sub_label+'_len.fa')
            df.to_csv(sub_label+"_subgroup.csv",sep='\t',index=False)
    else:
        sub_label=protaxpar_label+'_whole'
        print("Align for all the proteins.\n")
        shutil.copyfile(protaxpar_label+'.fa',sub_label+'_len.fa')
        df.to_csv(sub_label+"_subgroup.csv",sep='\t',index=False)
    print("The subset table is in "+sub_label+"_subgroup.csv file.\n")
    return sub_label


# determine, and plot, the level of protein sequence conservation
def conservation(dataset,sub_label):
    protein=dataset[0]
    taxonomy=dataset[1]
    # sequence alignment
    print("Aligning sequences...\n")
    # clustalo, --auto  Set options automatically,  --force  Force file overwriting
    os.system('clustalo -i '+sub_label+'_len.fa -o '+sub_label+'.align.msf --outfmt=msf --auto --force --threads=100')
    # basic information for the alignment
    print("Summarizing alignment results...\n")
    os.system('infoalign '+sub_label+'.align.msf -out '+sub_label+'.infoalign_ing')
    # orgainze the delimiter, change all spaces to tab
    os.system("sed 's/   */\t/g' "+sub_label+".infoalign_ing | sed 's/\t\t*/\t/g' > "+sub_label+".infoalign")
    # dataframe sorted by "Ident" from high to low
    aligninfo=pd.read_csv(sub_label+'.infoalign', sep="\t",na_values=[''])
    aligninfo=aligninfo.sort_values('Ident',ascending=False).reset_index(drop=True)
    # only save useful columns
    aligninfo=aligninfo[['Name', 'SeqLen', 'AlignLen', 'Gaps', 'GapLen','Ident', 'Similar', 'Differ', '% Change']]
    # add species column to the aligninfo dataframe
    spe_df=pd.read_csv(dataset[2]+"_summary_species.csv", sep="\t")
    med_df=pd.DataFrame({'Name':list(spe_df['protein_accession']),'species':list(spe_df['species'])})
    name=[]
    for i in list(med_df['Name']):
        name.append(i.replace('.','_'))
    med_df=pd.DataFrame({'Name':name,'species':list(spe_df['species'])})
    align_df=pd.merge(aligninfo,med_df,on=['Name'],how='left')
    print('Here is the first 20 rows of aligment summary sorted by identity, from high to low.')
    print(align_df.head(20))
    align_df.to_csv(sub_label+".infoalign.csv",sep='\t')
    print("Full result inside \""+sub_label+".infoalign.csv\" file.")
    os.remove(sub_label+'.infoalign_ing')
    os.remove(sub_label+'.infoalign')
    # generate consensus sequence
    os.system('cons '+sub_label+'.align.msf '+sub_label+'_align.cons')
    aligncons=open(sub_label+'_align.cons').read()
    print("Consensus sequence inside \""+sub_label+"_align.cons\" file.\nHere is consensus sequence:")
    print(aligncons)
    # plot the the level of conservation
    plot="plotcon -sformat msf "+sub_label+".align.msf -winsize 20 -gsubtitle '"+protein+" in "+taxonomy+" group' -graph png -goutfile "+sub_label+"_plotcon"
    os.system(plot)
    print("conservation image file generated.\n") # * in case running single function mutiple times
    # display the plot
    os.system('display '+sub_label+'_plotcon*.png') # sometimes "display" just doesn't work
    print("If your terminal can't open image, try open the image yourself or restart a terminal window.\n")


# output fasta file for each sequence
# only
def seq_out(sub_label):
    fas=open(sub_label+'_len.fa').read().strip()
    full_seq=fas.split('>')[1:] # delete first empty string
    filename=[]
    full_fa=[]
    for i in list(range(0,len(full_seq),1)): 
        # output fasta file for each sequence
        filename.append(full_seq[i].split(' ')[0])  # set sequence name as file name
        full_fa.append('>'+full_seq[i])   # add separator('>') back
        with open(filename[i]+'_'+sub_label+'_len_single_seq.fa',"w") as outfile:
            outfile.write(full_fa[i]) # write in single sequence

# scan each sequence with motifs from PROSITE database
def find_motifs(sub_label):
    # load dataframe
    df=pd.read_csv(sub_label+"_subgroup.csv", sep="\t")
    # scan motif for each sequence. Parallel for multiple files in bash command is much quicker
    print("Finding motifs...")
    motif='for fn in *'+sub_label+'_len_single_seq.fa \n do patmatmotifs -full -sequence ${fn} -outfile ${fn//_seq.fa/}.patmatmotifs & \n done \n wait'
    os.system(motif)
    # if motifs are found or not
    if not re.search('len_single.patmatmotifs',','.join(os.listdir())):
        print('no motifs are found.')
    else:
    # get all the motif names from all files
        print("Summarizing motifs...")
        count_motifs="grep 'Motif =' *"+sub_label+"_len_single.patmatmotifs | cut -d ' ' -f 3 | sort > "+sub_label+"_len.motifs.count"
        os.system(count_motifs)
        # count the number of each motif in total
        motifs_n=open(sub_label+"_len.motifs.count").read()
        motifs_n=motifs_n.strip().split('\n')
        motifs_unique=set(motifs_n) 
        moti_count=[]
        moti_name=[]
        for uni in motifs_unique:
            moti_name.append(uni)
            moti_count.append(motifs_n.count(uni))
        # dataframe for motif counts, sort by the amount 
        motif_df=pd.DataFrame({"motif_name":moti_name,"amount":moti_count})
        motif_df=motif_df.sort_values('amount',ascending=False).reset_index(drop=True)
        print("Result of patmatmotifs is in "+sub_label+"_patmatmotifs_summary.csv file.")
        print(motif_df)
        # summary of motifs result and output
        motif_df.to_csv(sub_label+'_patmatmotifs_summary.csv',sep='\t',index=False)

