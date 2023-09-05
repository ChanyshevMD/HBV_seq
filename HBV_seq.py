#!/usr/bin/env python3
from Bio import SeqIO
from Bio import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
import fastq as fq
from Levenshtein import distance as lev
import random

FileNames = []

import os
for root, dirs, files in os.walk("."):  
    for filename in files:
        if '.fastq' in filename and 'R1' in filename: FileNames.append(filename)

#FileNames = FileNames[1:2]                             #process group of files
#FileNames = ['HBV-161_S59_L001_R1_001_merged.fastq']   #process given files
Amplicon_query = ''
SingleMode = 0
if len(FileNames) == 1:
    SingleMode = 1
    print('input amplicon number from 1 to 20 or press enter to whole assembly')
    Amplicon_query = input()
    if Amplicon_query != '': Amplicon_query = int(Amplicon_query)

FileOut = open('HBV_statistics.txt','w')      #file for statistics output
FileOut.write('Sample_Name,Genotype,Subtype,Seq,Coverage\n')
  
FastaOut = open('Fasta_out.fa','w')           #file for HBV consensus output 

Alleles_fr_all = []
Alleles_rv_all = []
for seq_record in SeqIO.parse("HBV_seq.fasta", "fasta"):
    Name = seq_record.description
    Alleles_fr_all.append([Name,seq_record.seq,seq_record.id])
    Alleles_rv_all.append([Name,seq_record.seq.reverse_complement()])

Alleles_fr = Alleles_fr_all[:8]
Alleles_rv = Alleles_rv_all[:8]
Alleles_fr_exact = Alleles_fr_all[8:49]

Start_fr = [  79, 263, 487, 709, 862,1006,1192,1299,1541,1679,1798,1967,2167,2321,2466,2685,2811,2896,3023,3257] #fr primers positions
Start_rv = [2922,2747,2546,2301,2157,2001,1863,1698,1566,1399,1285,1049, 849, 757, 547, 394, 306,  90,3245,3077] #rv primers positions
#Amplicon      1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20

Amplicon_length = [256, #HBV_01
                   247, #HBV_02
                   224, #HBV_03
                   247, #HBV_04
                   238, #HBV_05
                   250, #HBV_06
                   202, #HBV_07
                   260, #HBV_08
                   150, #HBV_09
                   185, #HBV_10
                   138, #HBV_11
                   241, #HBV_12
                   235, #HBV_13
                   173, #HBV_14
                   244, #HBV_15
                   178, #HBV_16
                   107, #HBV_17
                   238, #HBV_18
                   246, #HBV_19
                   170] #HBV_20

Forward = [['CTGGTGGCTCCAGTTCAG'],                                                              #HBV_1_Forward
           ['CACAGAGTCTAGACTCGTGGTG'],                                                          #HBV_2_Forward
           ['GCCCGTTTGTCCTCTAATTC'],                                                            #HBV_3_Forward
           ['CCATTTGTTCAGTGGTTCGTAG'],                                                          #HBV_4_Forward
           ['CACAAAACAAAAAGATGGGG'],                                                            #HBV_5_Forward
           ['GAAAGTATGTCAAAGAATTGTGGG'],                                                        #HBV_6_Forward
           ['CTGTGCCAAGTGTTTGCTGA'],                                                            #HBV-F7-2
           ['TCCTAGCCGCTTGTTTTG'],                                                              #HBV_8_Forward
           ['GGGCGCACCTCTCTTTAC'],                                                              #HBV_9_Forward
           ['CATAAGAGGACTCTTGGACTTTCAG'],                                                       #HBV_10_Forward
           ['GGAGGCTGTAGGCATAAATTG'],                                                           #HBV_11_Forward
           ['ATTGACCCGTATAAAGAATTTGG'],                                                         #HBV_12_Forward
           ['TGGGTGGGAAGTAATTTGG'],                                                             #HBV_13_Forward
           ['CGGAGTGTGGATTCGCAC'],                                                              #HBV-F14-2
           ['CAATCGCCGCGTCGCAGA'],                                                              #HBV-F15-2
           ['AAAAGAAGATTGCAATTGATTATGCC','AAAGAAGATTAAACTTAGTCATGCC'],                          #HBV-F16-3/4
           ['CTATTTACACACTCTATGGAAGGC','CATTATTTACATACTCTTTGGAAGGC','ATTTACACACCCTATGGAAGGC'],  #HBV-F17-2/3/4
           ['TGGAACAAGATCTACAGCATGGG'],                                                         #HBV-F18-1
           ['CGCAAATCCAGATTGGGAC','AGCAGTCCCGACTGGGAC'],                                        #HBV-F19-1/2
           ['CCTCAGGCCATGCAGTGG']]                                                              #HBV-F20-1
           
Reverse = [['GAGGACAACAGGTTGGTGAG'],                                                            #HBV_1_Reverse
           ['GCAGGTTCTGCATGGTCC','AGGTTCGGCAGGGTCC'],                                           #HBV-R2-2/3
           ['CCAGACAGTGGGGGAAAG'],                                                              #HBV_3_Reverse
           ['ATCAATAGGCCTGTTTACAGGAAG'],                                                        #HBV_4_Reverse
           ['GAAAGGCCTTGTAAGTTGGCGA','GAAAGGCCTTATACGTTGGCGA'],                                 #HBV-R5-2/3
           ['CAGTATGGATCGGCAGAGG'],                                                             #HBV_6_Reverse
           ['GCGCAGGATCCAGTTGGC','GCGAAGAATCCAGTTGGC'],                                         #HBV-R7-4
           ['GTCCGGCAGATGAGAAGG'],                                                              #HBV_8_Reverse
           ['AAGTATGCCTCAAGGTCGG'],                                                             #HBV_9_Reverse
           ['GCTTGGAGGCTTGAACAGTAG','TCACAGCTTGGAGGCTTGAA'],                                    #HBV_10_Reverse/R10-2
           ['CAAAAACGAGAGTAACTCCACAG','GGCAAAAAAGAGAACAATTCCACAG'],                             #HBV-R11-2/3
           ['TCCTGAACTTTAGGCCCATATTAGT','GCCTGATTTTTAAACCCATATTAGT'],                           #HBV-R12-2/3
           ['CGAGGGAGTTCTTCTTCTAGGG'],                                                          #HBV_13_Reverse
           ['CCCACCTTATGAGTCCAAGG'],                                                            #HBV-R14-2
           ['CCCTTATCCAATGGTAAATATTT','TTTCTCTAGGGGCAAATATTT','GATACCTTTGTCTAATGGTAAATATTT'],   #HBV-R15-2/3/4 cutted
           ['AGATCTTGTTCCCAAGAATATGG','TCGTGTTCCCGAGAATAAGG'],                                  #HBV_16_Reverse + HBV-R16-2
           ['TTCGGGAAAGAATCCCAGAGG','CGGGAACACATCCCAGAGG'],                                     #HBV-R17-2/3
           ['TATTGGTGGAGGCAGGAGG','CGATTGGTAGAAGCAGGAGG'],                                      #HBV-R18-1/2
           ['ATCTTGCAGAGTTTGGTGGAA','ATCTAGCAGAGCTTGATGAAA'],                                   #HBV-R19-1/2
           ['AACACGAGCAGGGGTCCTAG']]                                                            #HBV-R20-4

References1 = [['GAACAGTAAACCCTGTTCTGACTACTGCCTCTCCCTTATC'],                                                                                                                                  #HBV_1_Forward
               ['GACTTCTCTCAATTTTCTAGGGGGAGCTCCCGTGTGTCTT'],                                                                                                                                  #HBV_2_Forward
               ['CAGGATCTACGACCACCAGCACGGGACCATGCAAAACCTG','CAGGAACATCAACTACCAGCACGGGACCATGCAAGACCTG'],                                                                                       #HBV_3_Forward
               ['GGCTTTCCCCCACTGTCTGGCTTTCAGTTATATGGATGAT'],                                                                                                                                  #HBV_4_Forward
               ['CTATTCCCTTAACTTCATGGGATATGTAATTGGAAGTTGG','TTACTCTCTAAATTTTATGGGTTATGTCATTGGATGTTAT'],                                                                                       #HBV_5_Forward
               ['TCTTTTGGGCTTTGCTGCCCCTTTTACACAATGTGGATAT','TCTGTTGGGTTTCGCTGCTCCTTTCACCCAATGTGGTTAC'],                                                                                       #HBV_6_Forward
               ['CGCAACCCCCACTGGCTGGGGCTTGGTCATGGGCCATCAG','CGCAACCCCCACTGGCTGGGGCTTGGCCATAGGCCATCAG','CGCAACCCCCACTGGTTGGGGCTTGGCCATTGGCCATCAG'],                                            #HBV-F7-2
               ['CTCGCAGCAGGTCTGGAGCGAAACTTATCGGGACGGATAA','CTCGCAGCCGGTCTGGAGCGGACATTATCGGCACTGACAA'],                                                                                       #HBV_8_Forward
               ['GCGGACTCCCCGTCTGTGCCTTCTCATCTGCCGGACCGTG','ACTGCCTACCCGTCTGTGCCTTCTCATCTGACGGACCGTG'],                                                                                       #HBV_9_Forward
               ['CAATGTCAACGACCGACCTTGAGGCATACTTCAAAGACTG','TTATGTCAACAACCGGGGTGGAGAAATACTTCAAGGACTG','GACGGTCAATGACCTGGATCGAAGACTACATCAAAGACTG'],                                            #HBV_10_Forward
               ['GTCTGCGCACCAGCACCATGCAACTTTTTCACCTCTGCCT'],                                                                                                                                  #HBV_11_Forward
               ['AGCTACTGTGGAGTTACTCTCGTTTTTGCCTTCTGACTTC','AGCTTCTGTGGAATTGTTCTCTTTTTTGCCTTCTGACTTC'],                                                                                       #HBV_12_Forward
               ['AAGATCCAGCATCCAGGGACCTAGTAGTCAGTTATGTCAA','AGGATCCTGCGGCTAGAGATCTAGTAGTTAATTATGTCAA'],                                                                                       #HBV_13_Forward
               ['TCCTCCAGCTTATAGACCACCAAATGCCCCTATCCTATCA','ACCCCCCGCATATAGGCCACAAAATGCCCCTATCTTATCA'],                                                                                       #HBV-F14-2
               ['AGATCTCAATCTCGGGAATCTCAATGTTAGTATTCCTTGG','AGATCTGCATCTCCAGCTTCCCAATGTTAGTATTCCTTGG'],                                                                                       #HBV-F15-2
               ['TGCCAGGTTTTATCCAAAGGTTACCAAATATTTACCATTG','TGCTAGGTTCTATCCCAATGTTACTAAATATTTGCCCTTA','TGCTAGATTTTTTCCCATCTCTACGAAATATTTGCCCCTA'],                                            #HBV-F16-3/4
               ['TGGTATTCTATATAAGCGGGAAACCACACGTAGCGCATCA','GGGTATCTTATATAAAAGAGAGACAACACGTAGCGCCTCA','GGGCATCTTATATAAAAGAGAAACTACACGGTACGCCTCA'],                                            #HBV-F17-2/3/4                                                                                     #HBV-F17-2/3/4
               ['GCAGAATCTTTCCACCAGCAATCCTCTGGGATTCTTTCCC','AGGTTGGTCATCAAAACCTCGCAAAGGCATGGGGACGAAT','AGCACCTCTCTCAACGACAAGAAGGGGCATGGGACAGAAT','GCTTTCTTGGACGGTCCCTCTCGAATGGGGGAAGAATCAT'], #HBV-F18-1
               ['TTCAATCCCAACAAGGACACCTGGCCAGACGCCAACAAGG','TTCAACAAAAACAAGGACACTTGGCCAATGGCAAACAAGG','TTCAACCCCATCAAGGACCACTGGCCACAAGCCAATCAGG','TTCAACCCCAACAAGGACCCTTGGCCAGAGGCCTGGCAGG'], #HBV-F19-1/2
               ['AATTCCACAACCTTCCACCAAACTCTGCAAGATCCCAGAG','AACTCAACTCACTTCCACCAGGCTCTGTTGGATCCGAGGG','AACTCTACAGCATTCCACCAAGCTCTACAAAATCCCAAAG']]                                            #HBV-F20-1
               
References2 = [['TGATTGGAGGTTGGGGACTGCGAATTTTGGCCAAGACACA'],                                                                                                                                  #HBV_1_Reverse
               ['CGTACTGGTGGTTGATGATCCTGGAATTAGAGGACAAACG','CGTGCTGGTGGTTGTAGATCCTGGAAGTAGAGGACACACG','CGTGCTGGTGGTTGAAGATCCTGGAATTAGAGGACAAACG'],                                            #HBV-R2-2/3
               ['CCCTACGAACCACTGAACAAATGGCACTAGTAAACTGAGC'],                                                                                                                                  #HBV_3_Reverse
               ['TTTTCTAAAACATTCTTTGATTTTTTGTATGATGTGTTCT','TTTCCTAAAACATTGTTTGATTTTTAACACAATATGTTCC','TTTTCTGAAACATTCTTTGATTTTCTGTGTGATGTGATCT'],                                            #HBV_4_Reverse
               ['GAAAGTGAAAGCCTGTTTAGCTTGTATACATGCATACAAA','AAAAGTAAAAGCCTGCTTATATTGAATACATGCATATAAG','AAAGACAAAAGCCTGTTTCGCAGTGATACAGTTATACAGA'],                                            #HBV-R5-2/3
               ['AGCCACAAAGGTTCCACGCATGCGCTGATGGCCTATGGCC'],                                                                                                                                  #HBV_6_Reverse
               ['AGCACAGCCTAGCAGCCATGGAAACGATGTATATTTGCGG','AGTACAACCTAGCAGCCATGGGAAAGATGTATATTTCCGA','AGCACAGCCTAGCAGCCATGGGAAGGAGGTGTACTTCCGA'],                                            #HBV-R7-4
               ['CACAGACGGGGAGTCCGCGTAAAGAGAGGTGCGCCCCGTG','CACAGACGGGTAGGCAGTGGAATGAGAGGTGCGCTCCGTG'],                                                                                       #HBV_8_Reverse
               ['TCGTTGACATTGCTGAGAGTCCAAGAGTCCTCTTATGTAA','TCGTTGACCGGGGCGAAAGTCCAAGAGTCCTTTTATGTAA','TTGTTGACATAACAAACAGTCCAAGAGTCCACTTATATAA'],                                            #HBV_9_Reverse
               ['GACATGAACAAGAGATGATTAGGCAGAGGTGAAAAAGTTG','CAGTAGGACATGAACAAGAGATGATTAGGCAGAGGTGAAA'],                                                                                       #HBV_10_Reverse/R10-2
               ['TAGCTCCAAATTCTTTATAAGGGTCGATGTCCATGCCCCA','TAGCTCCAAATTCTTTATAAGGGTCAATGTCTAAGCCAAA','AAGCTCCAAATTCTTTATACGGGTCAATGTCCATGCCCCA'],                                            #HBV-R11-2/3
               ['GTTGACATAACTGACTACTAGGTCTCTAGACGCTGGATCT','ATTAACATAATTGACTACTAGATCCCTGGATGCTGGATCT','GTTGACATAATTAACTACTAGATCTCTAGCCGCAGGATCC'],                                            #HBV-R12-2/3
               ['GACCTGCCTCGGTCCCGTCGTCTAACAACAGTAGTTTCCG','GACCTGCCTCTTCGTCTAACAACAGTAGTTTCCGGAAGTG'],                                                                                       #HBV_13_Reverse
               ['AATACTAACATTGAGATTCCCGAGATTGAGATCTTCTGCG','AATACTAACATTGGGAAGCTGGAGATGCAGATCTTCTGCG'],                                                                                       #HBV-R14-2
               ['GGTAACCTTTGGATAAAACCTGGCAGGCATAATCAATTGC','AGTGTGGGTAGGATAGAATCTAGCAGGCATAATTAATTTC','AGTAACATTGGGATAGAACCTAGCAGGCATAATCAATTTT','CGTAGAGATGGGAAAAAATCTAGCAGGCATGACTAAGTTT'], #HBV-R15-2/3/4 cutted
               ['TGACCCACAAAATGAGGCGCTACGTGTTGTTTCTCTCTTA','TGACCCACAGAATGAGGCGTACCGTGTAGTTTCTCTTTTA'],                                                                                       #HBV_16_Reverse + HBV-R16-2
               ['ATTGGGAACAGAAAGATTCGTCCCCATGCCTTTGCGAGGT','ATTGGGCACAGAGAGATTCTGTCCCATGCCCCTTCTTGTC','ATTGGTGGTGGAATGATTCTTCCCCCATTCGAGAGGGACC','ATTGCTGGTGGAAAGATTCTGCCCCATGCTGTAGCTCTTG'], #HBV-R17-2/3 
               ['CGGATTTGCTGGCAAAGTTTGTAGTATGCCCTGAGCCTGA','AGGAATTGTTGACACTGTGGTCAATATGCCCTGAGCCTGA','CGGATCTGCTGGCAATGTTGTAAGCACACCCTGTGCCTGT','TGGATCTGGTGGCGAGGTTGTCAGAATGCCCTGTGCCTGA'], #HBV-R18-1/2
               ['GGTTGTGGAATTCCACTGCATGGCCTGAGGATGAGTGTTT','CGTGGTGGAGTTCCACTGCATGGCCTGAGGATGAGTATCC','GTGAGTTGAGTTCCACTGCATTGCCTGTGGATGTGTGTCT','CTGTGTTGAGTTCCACTGCATGGCCTGTGGATGTGTGTCC'], #HBV-R19-1/2
               ['GAATCCTGATGTGATGTTCTCCATGTTCAGCGCAGGGTCC','GAGTCCTGATGTGATGTTGTCCATGTTCATAGCAGGGCCC','GAGTCCTGATGCGATGTTCTCCATGTTCGGTACAGGGTCC']]                                            #HBV-R20-4
 
def CHECK(seq1,seq2):                 #function to check similarity of two sequences; one wrong nucleotide is acceptable
    if lev(seq1,seq2) > 1: return 0
    else: return 1   

def CONSENSUS(List_of_reads):         #function to make one consensus seq from reads of one amplicon
    #set same length
    Lengths = []
    for i in range(len(List_of_reads)):
        Lengths.append(len(List_of_reads[i].body))
    for i in range(len(List_of_reads)):
        List_of_reads[i].body = List_of_reads[i].body + '-'*(max(Lengths)-len(List_of_reads[i].body)) #make one length for reads
        List_of_reads[i].qstr = List_of_reads[i].qstr + '-'*(max(Lengths)-len(List_of_reads[i].qstr)) #make one length for reads
    
    N = len(List_of_reads)
    consensus = ''     #resulting string sequence
    conservative = []  #share of minor nucleotides in each position
    phred_mean = []    #mean phred for each position
    
    for i in range(len(List_of_reads[0].body)): #i = letter number
        Temp = []
        Temp_phred = []
        for j in range(len(List_of_reads)):                                                  #record number            
            Temp_phred.append(ord(List_of_reads[j].qstr[i])-33)                              
            if ord(List_of_reads[j].qstr[i])-33 >= 20: Temp.append(List_of_reads[j].body[i]) #nucleotide is recorded only if phred >= 20
        Counts = [Temp.count('A'),Temp.count('C'),Temp.count('G'),Temp.count('T')]           #counting nucleotides
        Counts = np.array(Counts)/len(List_of_reads)                                         #normalization

        if max(Counts) == Counts[0]: consensus += 'A'
        elif max(Counts) == Counts[1]: consensus += 'C'
        elif max(Counts) == Counts[2]: consensus += 'G'
        elif max(Counts) == Counts[3]: consensus += 'T'
        conservative.append(sum(Counts)-max(Counts))   #record share of minor nucleotides

    #cutting by phred
    for i in range(len(phred_mean)):
        if phred_mean[i] < 30 and np.mean(phred_mean[i:])<25: #if mean phred for some position < 30 and mean phred onward < 25 - cut consensus here
            consensus = consensus[:i]
            conservative = conservative[:i]

    conservative = np.array(conservative)
    return [consensus,conservative]

for FileName in FileNames:
    print('---------------------------------')
    print(FileName)
    print('---------------------------------')
    AMPLICONS = []
    AMPLICONS_NUMBERS = []
    AMPLICONS_VARIABLES = []

    AMPLICONS_WITH_VARIABLE_PLACES = []
    GENOTYPES_BY_AMPLICONS = []
    GENOTYPES_PERCENT = []
    MATRIX_PERCENT = np.zeros((9, 20))

    A_RANGE = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19] #in normal condition handle all amplicons
    if Amplicon_query != '': A_RANGE = [Amplicon_query-1]         #if one amplicon is input
    
    for a in A_RANGE: #amplicon number
        R1 = 1
        R2 = 1
        k = 0
        RES_fr = []
        RES_rv = []
        Record_fr = []
        Record_rv = []

        F1 = fq.read(FileName)
        F2 = fq.read(FileName.replace('R1','R2'))

        F1 = list(F1)
        F2 = list(F2)

        AllReads = []
        for i in range(len(F1)):
            AllReads.append([F1[i],F2[i]])
        random.shuffle(AllReads)

        F1 = []
        F2 = []
        for i in range(len(AllReads)):
            F1.append(AllReads[i][0])
            F2.append(AllReads[i][1])

        AllReads = []

        for i in range(len(F1)):
            if len(RES_fr) == 1000: break   #for time saving consider no more than 1000 random reads of one amplicon
            R1 = F1[i]
            R2 = F2[i]
            Found = False
            for fr in Forward[a]:
                if CHECK(fr,R1.body[:len(fr)]) == 1 and Found == False:              #check, if forward read starts with correct primer
                     for rv in Reverse[a]:
                         if CHECK(rv,R2.body[:len(rv)]) == 1 and Found == False:     #check, if reverse read starts with correct primer
                             for rf1 in References1[a]:
                                 for rf2 in References2[a]:
                                     L1 = lev(str(R1.body[len(fr):len(fr)+40]), rf1) #check, if 40 nt after forward primer are similar with HBV
                                     L2 = lev(str(R2.body[len(rv):len(rv)+40]), rf2) #check, if 40 nt after reverse primer are similar with HBV
                                     if L1 <= 12 and L2 <= 12 and Found == False:    
                                         Found = True
                                         R1.body=R1.body[len(fr):] #cut primers
                                         R2.body=R2.body[len(rv):] #cut primers
                                         R1.qstr=R1.qstr[len(fr):] #cut primers
                                         R2.qstr=R2.qstr[len(rv):] #cut primers

                                         if len(R1.body) > Amplicon_length[a]: R1.body=R1.body[:Amplicon_length[a]] #cut by length
                                         if len(R2.body) > Amplicon_length[a]: R2.body=R2.body[:Amplicon_length[a]] #cut by length
                                         if len(R1.qstr) > Amplicon_length[a]: R1.qstr=R1.qstr[:Amplicon_length[a]] #cut by length
                                         if len(R1.qstr) > Amplicon_length[a]: R2.qstr=R2.qstr[:Amplicon_length[a]] #cut by length
                                         
                                         RES_fr.append(R1)
                                         RES_rv.append(R2)

        F1 = []
        F2 = []

        N = len(RES_fr)
        if N < 10:                  #not enough reads
            AMPLICONS.append('')
            GENOTYPES_BY_AMPLICONS.append(8)
            AMPLICONS_NUMBERS.append(N)
            print(str(a+1)+' amplicon not enough reads')
            continue
 
        Genotypes_counts = [0 for i in range(8+1)]
        Genotypes_fr = [[] for i in range(8+1)]
        Genotypes_rv = [[] for i in range(8+1)]
              
        for i in range(len(RES_fr)):
            Scores = []
            for j in range(8):
                ref = str(Alleles_fr[j][1][Start_fr[a]:Start_fr[a]+len(RES_fr[i])+400])
                ref = ref.replace('-','')
                if len(ref) < len(RES_fr[i]):
                    ref = Alleles_fr[j][1]+Alleles_fr[j][1]
                    ref = str(ref[Start_fr[a]:Start_fr[a]+len(RES_fr[i])+400])
                    ref = ref.replace('-','')
                ref = ref[:len(RES_fr[i])]                                
                Scores.append(lev(RES_fr[i].body,ref))
                
            for j in [3,0,1,2,4,5,6,7]:
                if Scores[j] == min(Scores):
                    if Scores[j]/len(RES_fr[i]) < 0.2: #similarity check
                        Genotypes_counts[j] += 1
                        Genotypes_fr[j].append(RES_fr[i])
                        Genotypes_rv[j].append(RES_rv[i])
                        break
                    else:
                        Genotypes_counts[8] += 1
                        Genotypes_fr[8].append(RES_fr[i])
                        Genotypes_rv[8].append(RES_rv[i])
                        break

        for i in range(9):
            MATRIX_PERCENT[i][a] = Genotypes_counts[i]/N

        for i in range(8):
            if Genotypes_counts[i] == max(Genotypes_counts[:8]):
                RES_fr = Genotypes_fr[i]
                RES_rv = Genotypes_rv[i]
                GENOTYPES_BY_AMPLICONS.append(i)
                AMPLICONS_NUMBERS.append(max(Genotypes_counts[:8]))
                break

        consensus_fr,conservative_fr = CONSENSUS(RES_fr)
        consensus_rv,conservative_rv = CONSENSUS(RES_rv)
        consensus_rv = Seq.Seq(consensus_rv).reverse_complement()

        #checking for insdels
        for j in range(len(conservative_fr)-20):
            if conservative_fr[j] > 0.25 and conservative_fr[j+1] > 0.25 and conservative_fr[j+2] > 0.25:
                if np.mean(conservative_fr[j:]) > 0.25:
                    print(str(a+1)+' amplicon insertion/deletion warning')

                    INDEL_fr = [[RES_fr[0]]]
                    INDEL_rv = [[RES_rv[0]]]
                    criteria = int(0.865*len(RES_fr[0]))
                    for j in range(1,len(RES_fr)):
                        flag = 0
                        for k in range(len(INDEL_fr)):
                            if flag == 0:
                                Align = pairwise2.align.localms(RES_fr[j].body, INDEL_fr[k][0].body, 1, 0, -5, -5)
                                if Align[0].score > criteria:
                                    INDEL_fr[k].append(RES_fr[j])
                                    INDEL_rv[k].append(RES_rv[j])
                                    flag = 1
                        if flag == 0:
                            INDEL_fr.append([RES_fr[j]])
                            INDEL_rv.append([RES_rv[j]])

                    INDEL_NUMBERS = []
                    for j in range(len(INDEL_fr)): INDEL_NUMBERS.append(len(INDEL_fr[j]))
                    INDEL_NUMBERS.sort(reverse = True)
                    for j in range(len(INDEL_fr)):
                        if len(INDEL_fr[j]) == INDEL_NUMBERS[0]:
                            RES_fr_new = INDEL_fr[j]
                            RES_rv_new = INDEL_rv[j]

                    consensus_fr,conservative_fr = CONSENSUS(RES_fr_new)
                    consensus_rv,conservative_rv = CONSENSUS(RES_rv_new)
                    consensus_rv = Seq.Seq(consensus_rv).reverse_complement()
                    
                    break
        
        ref = Alleles_fr[i][1].replace('-','')        
        
        Align_fr = pairwise2.align.localms(ref,consensus_fr, 1, 0, -5, -1)
        Align_rv = pairwise2.align.localms(ref,consensus_rv, 1, 0, -5, -1)
        if Align_fr[0].seqA == Align_rv[0].seqA and a != 18:
            consensus = Align_fr[0].seqB[Align_fr[0].start:Align_fr[0].end] + Align_rv[0].seqB[Align_fr[0].end:Align_rv[0].end]
        elif Align_fr[0].seqA == Align_rv[0].seqA[:len(Align_fr[0].seqA)] and a == 18:
            consensus = Align_fr[0].seqB[Align_fr[0].start:Align_fr[0].end] + Align_rv[0].seqB[Align_fr[0].end:Align_rv[0].end]
        else:
            print(str(a+1)+' amplicon alignment warning')
            consensus = Align_fr[0].seqB[Align_fr[0].start:Align_fr[0].end] + Align_rv[0].seqB[Align_fr[0].end:Align_rv[0].end]

        consensus = consensus.replace('-','')
        AMPLICONS.append(consensus)

    if Amplicon_query != '':
        print('consensus made from fr reads')
        print(consensus_fr)
        print('consensus made from rv reads')
        print(consensus_rv)
        print('final amplicon consensus')
        print(consensus)
        print('print random 20 reads? input f or r for forward or reverse reads or enter to skip')
        Choice_reads = input()
        if Choice_reads == 'f':
            for i in range(min(20,len(RES_fr))): print(RES_fr[i].body)
        if Choice_reads == 'r':
            for i in range(min(20,len(RES_rv))): print(RES_rv[i].body)
        
        continue #skip genome assemble if one amplicon is query

    #ASSEBLING GENOME
    AMPLICONS_TRUE_GENOTYPE = []
    max_item = lambda s: max(t := {i: s.count(i) for i in s}, key=t.get)
    True_genotype = max_item(GENOTYPES_BY_AMPLICONS)

    if True_genotype == 8: True_genotype = 3

    ref = Alleles_fr[True_genotype][1].replace('-','')  #reference genome of true genotype

    Result = []                                         #list of amplicons with cut overlaps

    AMPLICONS = [AMPLICONS[19]]+AMPLICONS[:19]                         #20th amplicon becomes 1st
    AMPLICONS_NUMBERS = [AMPLICONS_NUMBERS[19]]+AMPLICONS_NUMBERS[:19]

    End = 0
    i = 0
    Result.append(AMPLICONS[i])
    if AMPLICONS[i] != '':
        Align = pairwise2.align.localms(ref, AMPLICONS[i], 1, 0, -1, -1) #align 1st amplicon to reference, record position of end
        End = Align[0].end

    for i in range(1,20):
        if AMPLICONS[i] == '':
            Result.append(AMPLICONS[i]) 
        if AMPLICONS[i] != '':
            Align = pairwise2.align.localms(ref, AMPLICONS[i], 1, 0, -1, -1) #align i amplicon to reference
            Result.append(Align[0].seqB[Align[0].start:Align[0].end])        
            #cutting ends 
            if Align[0].start < End and Result[i-1] != '':
                if AMPLICONS_NUMBERS[i-1]>AMPLICONS_NUMBERS[i] and AMPLICONS_NUMBERS[i-1]>100 and AMPLICONS_NUMBERS[i]<10:
                    Result[i] = Result[i][(End-Align[0].start):] #in case that left amplicon is OK and right has little number of reads, overlap is taken by left
                else:
                    Result[i-1] = Result[i-1][:-(End-Align[0].start)]   #in normal case overlap region is taken by right amplicon and left is cut

            End = Align[0].end

    R = ''                   #resulting consensus string
    for r in Result: R += r

    print('Consensus sequence')
    print(R)
    print('Coverage')
    print('%.2f'%(100*len(R.replace('-',''))/len(ref)))
    Genotype_letters = ['A','B','C','D','E','F','G','H']

    #exact allele
    Genotypes_exact = []
    Genotype_exact = ''

    scores = []
    for j in range(40):
        scores.append(lev(R,str(Alleles_fr_exact[j][1])))
    for j in range(40):
        if scores[j] == min(scores):
            Genotypes_exact.append(Alleles_fr_exact[j][2])
    Genotypes_exact = list(set(Genotypes_exact))
    for G in Genotypes_exact: Genotype_exact += G
    print('Genotype')
    print(Genotype_exact)

    S = FileName+','+Genotype_letters[True_genotype]+','+Genotype_exact
    S += ','+R+','+str(100*len(R.replace('-',''))/len(ref))+','
   
    S = S[:-1]

    FileOut.write(S+'\n')
    FastaOut.write('>'+FileName+'\n')
    FastaOut.write(R+'\n')

FileOut.close()
FastaOut.close()
