'''
ORForise is a tool used for comparison of Gene Prediction results from two different tools. This code takes the following input:
gene prediction results files from Prodigal and FragGeneScan
It gives the following output:
csv file containing the comparison of the tools based on various metrics. 

Here is the link to GitHub of ORForise: https://github.com/NickJD/ORForise


As this code required to be implemented locally, I hardcoded the file paths. 

'''

import subprocess
import glob
import pathlib

def prodigal():
    prodigal_files=glob.glob('../prodigal/prodigal_output/*') #location of prodigal output files- include the path where your Prodigal output files are
    dna_files=glob.glob('../tests/dna_files/*') #location of test contigs- edit to the path where your test contig .fasta files are
    contig_ids=[x.split("/")[-1].split("_")[0] for x in dna_files] #extracting the contig ids
    prodigal_dir = pathlib.Path(prodigal_files[0]).parent
    dna_file_dir=pathlib.Path(dna_files[0]).parent
    subprocess.run(['mkdir','-p','/home/orf/prodigal_orf'],stdout=2, stderr=subprocess.STDOUT) #edit the path to the directory where you want to store your output 
    for id in contig_ids:
        print("On contig- %s"%id)
        process=subprocess.Popen(['Aggregate-Compare','-dna',dna_file_dir /f'{id}_contigs.fasta','-ref','/home/orf/genomic.gff','-t', 'Prodigal','-tp',prodigal_dir /f'{id}_contigs.gff','-o','/home/orf/prodigal_orf/'f'{id}_prodigal_orf.csv']) #edit paths


def fraggenescan():
    fraggenescan_files=glob.glob('../fraggenescan/full_contigs/*') #edit to the path of your FragGeneScan results
    dna_files=glob.glob('/home/dna_files/*') #edit to the path of your dna test contig .fasta files
    contig_ids=[x.split("/")[-1].split("_")[0] for x in dna_files] #extracting the contig ids
    fraggenescan_dir = pathlib.Path(fraggenescan_files[0]).parent
    dna_file_dir=pathlib.Path(dna_files[0]).parent
    subprocess.run(['mkdir','-p','/home/orf/fraggenescan_orf'], stdout=2, stderr=subprocess.STDOUT) #edit the path to where you wish to store your output
    for id in contig_ids:
        print("On contig- %s" % id)
        process=subprocess.Popen(['Aggregate-Compare','-dna',dna_file_dir /f'{id}_contigs.fasta','-ref','/home/orf/genomic.gff','-t', 'FragGeneScan','-tp',fraggenescan_dir /f'{id}.gff','-o','/home/orf/fraggenescan_orf/'f'{id}_fraggenescan_orf.csv']) #edit paths 

if __name__ == "__main__":
    prodigal()
    fraggenescan()
