#!/usr/bin/env python
#__________________________________________________________________DNA BARCODE DATA ANALYSIS SCRIPT________________________________________________________________
#
#followup script of the Barcode_Script_1.4 to analyze the data of the output txt files.
#it generates 3 output-files: 1) 100 top hits of DNA-sequences, 2) 100 top hits of peptide sequences & 3) 100 top hits, sorted in the different libraries
#written by Sabrina Weis, internship March/April 2017
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------
import os
import re
import time

genetic_code= {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

def Codons2Protein(ORF, genetic_code):
    protein=''
    for i in range(0, len(ORF), 3):
        codon=ORF[i:i+3]
        protein += genetic_code[codon]
    return protein

execfile("//netfile1.bioquant.uni-heidelberg.de/RedirDesktops/bq_jhaar/Desktop/NGS Script/Barcode_Script_Janina.conf")

#open the attribution (txt) file and store filenames and the corresponding samplenames in a dictionary called "sample"
with open (filename_sample_file) as temp:
    sample=dict(line.strip().split() for line in temp if line.strip())

#all filenames of the given folder are extracted and those which ends with txt (outputfiles from the Barcode-Script) are stored in the list inputfiles
inputfiles=[]
objects=os.listdir(my_dir)
for i in range(0,len(objects)):
    filename=objects[i]
    if filename.endswith('txt'):
        if filename in sample.keys():
            inputfiles.append(filename)

#the variable "results" is defined as a dictionary containing a dictionary for each file stored in the list filenames
results={}
for l in range(0,len(inputfiles)):
    results[inputfiles[l]]={}

#the files in the list inputfiles are opened in succession and each line (containing a variant and its frequency) is storead in the dictionary "results"
for n in range(0,len(inputfiles)):
    with open(my_dir+inputfiles[n], 'r') as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    for l in content[6:-1]:
            splitlines=l.split('\t')
            results[inputfiles[n]][splitlines[0]]=float(splitlines[1])
info={}
for n in range(0,len(inputfiles)):
    with open(my_dir+inputfiles[n], 'r') as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    for m in content[3:5]:
            splitlines=m.split(': ')
            info[splitlines[0]]=float(splitlines[1])

#each filename that occurs in "results" is replaced by the corresponding samplename (defined in the attribution file)
for n in sample.keys():
    for i in results.keys():
        if n in i:
            results[sample[n]]=results[i]
            del results[i]

#sum of all reads is calculated for each sample to enable the determination of the proportions
results_sum_counts={}
for n in results:
    results_sum_counts[n]=sum(results[n].values())

#all DNA-sequences are translated into peptide-sequences and are stored in a new dictionary "peptides"
peptides={}
for n in results:
    peptides[n]={}
    for i in results[n].keys():
        seq=Codons2Protein(i, genetic_code)
        #if the peptide-sequence already exists, the frequencies of the DNA-sequences that lead to the same peptide-seqence are added
        if seq in peptides[n]:
            peptides[n][seq]+=results[n][i]
        else:
            peptides[n][seq]=results[n][i]

#the peptide-sequences are sorted in the the 4 difderent libraries SGXXXXXX, SGNXXRXXX, RGNXXRXXX, RGXXXXXXX using regular expressions
library1={}
for n in peptides:
    library1[n]={}
    for i in peptides[n].keys():
        match1=re.search('SG[A-M,O-Z,*][A-Z,*]{2}[A-Z,*]{4}',i)
        match2=re.search('SG[A-Z,*]{3}[A-Q,S-Z,*][A-Z,*]{3}',i)
        if match1 or match2:
            library1[n][i]=peptides[n][i]

library2={}
for n in peptides:
    library2[n]={}
    for i in peptides[n].keys():
        match=re.search('SGN[A-Z,*]{2}R[A-Z,*]{3}',i)
        if match:
            library2[n][i]=peptides[n][i]

library3={}
for n in peptides:
    library3[n]={}
    for i in peptides[n].keys():
        match=re.search('RGN[A-Z,*]{2}R[A-Z,*]{3}',i)
        if match:
            library3[n][i]=peptides[n][i]

library4={}
for n in peptides:
    library4[n]={}
    for i in peptides[n].keys():
        match1=re.search('RG[A-M,O-Z,*][A-Z,*]{2}[A-Z,*]{4}',i)
        match2=re.search('RG[A-Z,*]{3}[A-Q,S-Z,*][A-Z,*]{3}',i)
        if match1 or match2:
            library4[n][i]=peptides[n][i]

#get date and time of the analysis (used in the filename of the outputfiles)
timestr = time.strftime("%Y-%m-%d_%H-%M-%S")

#create 3 different outputfiles for each sample
for n in results:
    #"DNA-seq" outputfile containing the 100 most present sequences. Most frequenctly occuring DNA-seq at the top
    outputfile_DNA=my_dir+n+'_DNA-seq_'+timestr+'.txt'
    d=open(outputfile_DNA, 'w')
    d.write('\n'+'Total number of reads:'+' '+str(info['Total number of reads']))
    d.write('\n'+'Reads recovered:'+' '+str(info['Reads recovered']))
    d.write('\n'+'Number of different DNA-sequences: '+str(len(results[n])))
    d.write('\n\n'+'ranking\t'+'DNA-seq\t'+'#\t'+'proportion\t')
    results_sorted_keys=sorted(results[n], key=results[n].get, reverse=True)
    c=1
    for i in results_sorted_keys[:10000]:
        d.write('\n'+str(c)+'\t'+i+'\t'+str(results[n][i])+'\t'+str(results[n][i]/results_sum_counts[n]))
        c+=1
    d.close()
    #"Peptide-seq outputfile" containing the 100 most present sequences. Most frequenctly occuring DNA-seq at the top
    outputfile_Pep=my_dir+n+'_Pedtide-seq_'+timestr+'.txt'
    p=open(outputfile_Pep, 'w')
    p.write('\n'+'Total number of reads:'+' '+str(info['Total number of reads']))
    p.write('\n'+'Reads recovered:'+' '+str(info['Reads recovered']))
    p.write('\n'+'Number of different peptide sequences: '+str(len(peptides[n])))
    p.write('\n\n'+'ranking\t'+'peptide-seq\t'+'#\t'+'proportion\t')
    peptides_sorted_keys=sorted(peptides[n], key=peptides[n].get, reverse=True)
    c=1
    for i in peptides_sorted_keys[:10000]:
        p.write('\n'+str(c)+'\t'+i+'\t'+str(peptides[n][i])+'\t'+str(peptides[n][i]/results_sum_counts[n]))
        c+=1
    p.close()
    #"Libraries outputfile" containing the 100 most present peptide sequences of each library
    outputfile_Lib=my_dir+n+'_Libraries_'+timestr+'.txt'
    l=open(outputfile_Lib, 'w')
    l.write('\n'+'Total number of reads:'+' '+str(info['Total number of reads']))
    l.write('\n'+'Reads recovered:'+' '+str(info['Reads recovered']))
    l.write('\n'+'Number of different peptide sequences: '+str(len(peptides[n])))
    l.write('\n\nLibrary1\t\t\t'+'Library2\t\t\t'+'Library3\t\t\t'+'Library4\t\t\t')
    l.write('\n'+('Peptide-seq\t#\tproportion\t')*4)
    library1_sorted_keys=sorted(library1[n], key=library1[n].get, reverse=True)
    library2_sorted_keys=sorted(library2[n], key=library2[n].get, reverse=True)
    library3_sorted_keys=sorted(library3[n], key=library3[n].get, reverse=True)
    library4_sorted_keys=sorted(library4[n], key=library4[n].get, reverse=True)
    for c in range(0,9999):
        if len(library1_sorted_keys)>c:
            l.write('\n'+library1_sorted_keys[c]+'\t'+str(library1[n][library1_sorted_keys[c]])+'\t'+str(library1[n][library1_sorted_keys[c]]/results_sum_counts[n]))
        else:
            l.write('\n\t\t')
        if len(library2_sorted_keys)>c:
            l.write('\t'+library2_sorted_keys[c]+'\t'+str(library2[n][library2_sorted_keys[c]])+'\t'+str(library2[n][library2_sorted_keys[c]]/results_sum_counts[n]))
        else:
            l.write('\t\t\t')
        if len(library3_sorted_keys)>c:
            l.write('\t'+library3_sorted_keys[c]+'\t'+str(library3[n][library3_sorted_keys[c]])+'\t'+str(library3[n][library3_sorted_keys[c]]/results_sum_counts[n]))
        else:
            l.write('\t\t\t')
        if len(library4_sorted_keys)>c:
            l.write('\t'+library4_sorted_keys[c]+'\t'+str(library4[n][library4_sorted_keys[c]])+'\t'+str(library4[n][library4_sorted_keys[c]]/results_sum_counts[n]))
        else:
            l.write('\t')
    l.close()

