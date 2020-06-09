# Quantify-sequence-homo-heterozygosity

This code calculates the relative proportions of homo- and heterzygotes for a large data set of sequence alignments.
First parsimony-informative sites are identified, their nucleotide states queried, and then the highest- to lowest-frequency nucleotide states are ordered. The nucleotide state of either of the haplotypes for each individuals at each site is queried and compared with the reference nucleotide (most frequent) and alternate nucleotides (less frequent). Then the proportions of heterozygotes(coded 1) and homozygotes(coded 0,2) are averaged across all sequences
