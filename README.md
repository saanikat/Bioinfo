# Bioinformatics Codes
Here are some basic codes and pipelines used extensively in bioinformatics. Given below is a brief description for each of the provided codes. 
1. SNP_Calling.sh: This code takes FASTQ reads and a reference genome as an input. It aligns the FASTQ reads to the reference genome to create an alignment file. It then processes the alignment file followed by variant calling. This pipeline employs GATK (Genome Analysis Toolkit- https://gatk.broadinstitute.org/hc/en-us) for this. 
2. Gene_prediction.py: Gene Prediction is a crucial step in Computational Genomics. This code takes the contigs .fasta files as input which can be generated performing Genome Assembly. The code uses 3 different gene prediction tools and compares their CPU Usage, RAM, and duration. The tools compared are:
   1. Prodigal:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119
   2. FragGeneScan: https://pubmed.ncbi.nlm.nih.gov/20805240/
   3. Balrog: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008727
3. gene_prediction_comparison.py: ORForise is a tool used for comparison of Gene Prediction results from two different tools. This code takes the following input:
gene prediction results files from Prodigal and FragGeneScan
It gives the following output: csv file containing the comparison of the tools based on various metrics.
Here is the link to GitHub of ORForise: https://github.com/NickJD/ORForise
5. needleman_wunsch.py: The Needleman-Wunsch Alogirthm is a popular algorithm used to align protein or nucleotide sequences using Global Alignment. This code has been developed without importing numpy and is an easy to understand approach for biologists new to coding!
6. smith_watermann.py: Smith-Watermann Algorithm is a popular algorithm used to align protein or nucleotide sequences using Local Alignment. This code has been developed without importing numpy and is an easy to understand approach for biologists new to coding!
