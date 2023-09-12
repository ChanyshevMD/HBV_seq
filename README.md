# HBV sequencing panel
Hepatitis B virus (HBV) remains a pressing global public health concern. The clinical course of the disease, particularly its tendency towards chronicity and response to therapy, is significantly influenced by the HBV genotype and specific mutations.
The HBV panel developed for simple and robust sequencing of hepatitis B virus holds immense potential for utilization in research, epidemiological monitoring, and the advancement of personalized medicine approaches.
Present script is written to process data obtained by the NGS amplification panel for sequencing the hepatitis B virus.

Being placed in one folder with the raw NGS data, the script assembles HBV genomes and analyzes data. The script can be run from any Python IDE or from the command line in any operation system; it does not require the installation of any other programs other than Python packages BioPython, Levenshtein, Numpy, fastq. The use of the script does not require any knowledge in the field of programming and bioinformatics. At the same time, the script can be easily modified to work with another fasta file and other primers. 

## Usage
Check you have Python and modules Biopython, numpy, fastq, Levenshtein installed.
Download script HBV_seq.py and attached files HBV_seq.fasta, HBV-training_R1_001.fastq.gz, HBV-training_R2_001.fastq.gz into one empty folder.
Run script, when it prints "input amplicon number from 1 to 20 or press enter to whole assembly" press enter.

The output of the script should look like follows:

>HBV-training_R1_001.fastq.gz
>
>13 amplicon not enough reads

>Consensus sequence
>AATTCCACAACCTTCCATCAAACTCTGCAAGATCCC...AGAGTGAAAGGCCTGTATTTCCCTGCTGGTGGCTCCAGTTCAGGAACAGTAAACCCTGTTCCGACTACTGCCTCTCACATATCGTCAACCTTCTCGAGGATTGGGGACCCTGCGCTGAACATGGAGAACATCACATCAGGACTCCTAGGACCCCTGCTCGTGTTACAGGCGGGGGTTTTCTTGTTGACAAGAATCCTCACAATACCGCAGAGTCTAGACTCGTGGTGGACTTCTCTCAATTTTCTAGGGGGAACTACCGTGTGTCTTGGCCAAAATTCGCAGTCCCCAACCTCCAATCACTCACCAACCTCCTGTCCTCCAACTTGTCCTGGTTATCGTTGGATGTGTCTGCGGCGTTTTATCATCTTCCTCTTCATCCTGCTGCTATGCCTCATCTTCTTGTTGGTTCTTCTGGACTATCAAGGTATGTTGCCCGTTTGTCCTCTAATTCCAGGATCCTCAACCACCAGCACGGGACCATGCAGGACCTGCACGACTCCTGCTCGAGGAAACTCTACGTATCCCTCTTGTTGCTGTACCAAACCTTCGGACGGAAATTGCACCTGTATTCCCATCCCATCATCCTGGGCTTTCGGAAAATTCCTATGGGAGTGGGCCTCAGCCCGTTTCTCCTGGCTCAGTTTACTAGTGCCATTTGTTCAGTGGTTCGTAGGGCATTCCCCCACTGTTTGGCTTTCAGTTATATGGATGATGTGGTATTGGGGGCCAAGTCTGTACAGCACCTTGAGTCCCTTTTTACCGCTGTTACCAATTTTCTTTTGTCTTTGGGTATACATTTAAACCCTAACAAAACAAAGAGATGGGGCTATTCCCTAAATTTTATGGGTTATGTCATTGGATGTCATGGATCATTGCCACAAGAACACATCAGACAAAAAATCAAAGACTGTTTTAGAAAACTTCCTGTTAACAGGCCTATTGATTGGAAAGTGTGTCAAAGAATTGTGGGCCTTTTGGGTTTTGCTGCACCTTTTACACAATGTGGTTATCCTGCTTTAATGCCCTTGTATGCATGTATTCAATCTAAGCAGGCTTTCACTTTCTCGCCAACTTACAAGGCCTTTCTGTGTAAACAATACCTGAACCTTTACCCCGTTGCCCGGCAACGGCCAGGTCTGTGCCAAGTGTTTGCTGACGCAACCCCCACTGGCTGGGGCTTGGTCATAGGCCATCAGCGCATGCGTGGAACCTTTATGGCTCCTCTGCCGATCCATACTGCGGAACTCCTAGCCGCTTGTTTTGCTCGCAGCAGGTCTGGAGCAAACATTCTCGGGACTGATAACTCTGTTGTCCTCTCCCGCAAATATACATCGTTTCCATGGCTGCTAGGCTGTGCTGCCAACTGGATCCTGCGCGGGACGTCCTTTGTTTACGTCCCCTCGGCGCTGAATCCCGCGGACGACCCTTCTCGGGGTCGCTTGGGACTCTCTCGTCCCCTTCTCCGTCTGCCGTACCGACCGACGACGGGGCGCACCTCTCTTTACGCGGCCTCCCCGTCTGTGCCTTCTCATCTGCCGGACCGTGTGCACTTCGCTTCACCTCTGCACGTCGCATGGAGACCACCGTGAACGCCCACCAATTCTTGCCCAAGGTCTTACATAAGAGGACTCTTGGACTCTCTGCAATGTCAACGACCGACCTTGAGACATACTTCAAAGACTGTTTGTTTAAAGACTGGGAGGAGTTGGGGGAGGAGCTTAGATTAAAGGTCTTTGTACTAGGAGGCTGTAGGCATAAATTGGTCTGCGCACCAGCACCATGCAACTTTTTCACCTCTGCCTAATCATCTCTTGTTCATGTCCTACTGTTCAAGCCTCCAAGCTGTGCCTTGGGTGGCTTTAGGGCATGGACATCGACCCTTATAAAGAATTTGGAGCTTCCGTGGAGTTACTCTCGTTTTTGCCTTCTGACTTCTTTCCTTCAGTACGAGATCTTCTAGATACCGCCTCAGCTCTGTATCGGGATGCCTTAGAGTCTCCTGAACATTGTACACCTCACCATACTGCACTCAGGCAAGCTATTATTTGCTGGGGGGAATTAATGACTCTAGCTACCTGGGTGGGTGGTAATTTAGACGATCAAAGATCTAGGGAACTAGTAGTCGGTTATGTCAACTCCTCCAGCTTATAGACCACCAAATGCCCCTATCTTATCAACACTTCCGGAGACTACTGTTGTTAGATGCCGAGGCAGGTCCCCTAGAAGAAGAACTCCCTCGCCTCGCAGACGAAGATCTCAATCGCCGCGTCGCAGAAGATCTCAATCTAGGGAACCTCAATGTTAGTATTCCTTGGACTCATAAGGTGGGAAACTTTACGGGGCTTTATTCTTCTACTGTACCTGTCTTTAATCCTCATTGGAAAACACCCTCTTTTCCTAATATACATTTACACCAAGACATTATCAAAAAATGTGAACAGTTTGTAGGCCCACTTACAGTTAACGAAAAAAGAAGATTGCAATTGATTATGCCTGCTAGGTTTTATCCAAATGTCACCAAATATTTGCCATTGGATAAGGGTATTAAACCTTATTATCCAGAACATCTAGTTAATCATTACTTCCTAACTAGACACTATTTACACACTCTATGGAAGGCGGGTGTATTATATAAGAGAGAAACAACACATAGCGCCTCATTTTGTGGGTCACCATATTCTTGGGAACAAGAGCTACAGCATGGGGCAGAATCTCTCCACCAGCAATCCTCTGGGATTCTTTCCCGACCACCAGTTGGATCCAGCCTTCAGAGCAAACACTGCAAATCCAGATTGGGACTTCAATCCCAACAAGGACACTTGGCCGGACGCCAACAAGGTAGGAGCGGGAGCATTCGGGCTGGGTTACACCCCACCGCACGGAGGCCTTTTGGGGTGGAGCCCTCAGGCTCAGGGCATAATACAAACTTTGCCAGCAAATCCGCCTCCTGCCTCCACCAATCGCCAGTCAGGAAGGCAGCCTACCCCGCTGTCGCCTCCTTTGCGAGACACTCATCCTCAGGCCATGCAGTGG

>Coverage
>96.45

>Genotype
>D3

If you get this output, then the script is working and you can process real sequencing data.
In case of any problems with the script or any questions, please write to chanishq@gmail.com

The script is able to work in several modes. 
If the folder contains the data of several samples, the script processes each one in turn, prints name of fastq file, possible errors during data processing, resulting consensus sequence, coverage, and HBV genotype.
Also the script creates fasta file "Fasta_out.fa" containing consensus sequences and names. In addition, the script creates file "HBV_statistics.txt", containig data from console: names, genotypes, subtypes, consensus sequences, and coverages.

If the folder contains the data of one sample (for example HBV-training_R1_001.fastq.gz and HBV-training_R2_001.fastq.gz) the script asks if there is an amplicon of interest.
If answer is no (press Enter) the script assembles one sample as usual.
If answer is yes (input amplicon number) the script processes only amplicon of interest and prints consensus of forward reads, consensus of reverse reads and resulting consensus.
Also script asks to print 20 random forward or reverse reads so user may look at raw data.

## Samples_001:003 data
Sequencing data of the three samples described in the paper can be downloaded from the following links:

**1. HBV_seq data**

Data, obtained by present HBV sequence panel. Suitable for the script HBV_seq.py

[Samples_001-003_HBV_seq_data.rar](https://drive.google.com/file/d/1aKnKUC0SyxoT3cgZUKv8wRjUjgUJ1cbW/view)

**2. Nextera data**

Data, obtained by Nextera XT DNA Library Preparation Kit. Not suitable for the script HBV_seq.py, but may be processed by standard bioinformatics approach, as described in paper.
Attached as a proof of results obtained by HBV_seq.

[Samples_001-003_Nextera_data.rar](https://drive.google.com/file/d/1iPW5YDHStMVgH1oQIXTGBIeN79ZSuaH5/view)

**3. Sanger data**

Data of Sanger sequencing. Attached as a proof of results obtained by HBV_seq.

[Samples_001-003_Sanger_data.rar](https://drive.google.com/file/d/1F1ODrUiy5qQv-eIl7AU5Wzk0qDfS69vx/view)


