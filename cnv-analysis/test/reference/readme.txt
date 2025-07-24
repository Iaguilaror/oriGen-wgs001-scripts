## telomere_positions.txt
Downloaded from: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.37_GRCh38.p11/README.txt
zgrep on the file *_genomic_gaps.txt.gz

# centromeres.tsv
# centromeres downloaded by: https://www.biostars.org/p/435003/
Hi, thank you for the answer.

Eventually I got what I needed with few adjustments.

Here my solution to obtain centromeric coordinates for hg38:

Go to the Table Browser: http://genome.ucsc.edu/cgi-bin/hgTables
Choose the Mapping and Sequencing group
Select the "Chromosome Band (Ideogram)" track
Select filter, and enter "cen" in the gieStain field
Press "submit" and then "get output"
Each chromosome will have two entries which overlap.

They can be simply merged into a single entry.

I hope this could be helpful for others.

ADD REPLY • link4.2 years ago by Simo ▴ 50
0
As of 2021 the filter you want for this approach is "acen" rather than "cen", as noted by Simo above.


