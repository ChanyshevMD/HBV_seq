# HBV sequencing panel
Present script is written to process data obtained by the NGS amplification panel for sequencing the hepatitis B virus.
The details of library preparation are described in the paper in Virologica Sinica.
Hepatitis B virus (HBV) remains a pressing global public health concern. The clinical course of the disease, particularly its tendency towards chronicity and response to therapy, is significantly influenced by the HBV genotype and specific mutations.
The HBV panel developed for simple and robust sequencing of hepatitis B virus holds immense potential for utilization in research, epidemiological monitoring, and the advancement of personalized medicine approaches.

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
>96.4487743557511

>Genotype
>D3

If you get this output, then the script is working and you can process real sequencing data.
In case of any problems with the script or any questions, please write to chanishq@gmail.com
