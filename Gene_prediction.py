'''
Gene Prediction is an important step in computational genomics. This code takes the contigs .fasta files as input which can be generated performing Genome Assembly. The code uses 3 different gene prediction tools and compares their CPU Usage, RAM, and duration. 

input: contigs .fasta files
output: .gff and .faa files

The tools compared are:
1. Prodigal:https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119
2. FragGeneScan: https://pubmed.ncbi.nlm.nih.gov/20805240/
3. Balrog: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008727

Alternatively, if you simply want to generate the gene prediction results, you can comment out any of the tools you do not want results from.
'''

#!/usr/bin/env python3

import subprocess
import pathlib
from timeit import default_timer as timer
import psutil
import numpy as np
import threading
import os
import humanfriendly
import glob

# Change this to the directory where the assembly results are, without the trailing "/"
contig_dir = 'assembly/final_results'

pathlib.Path('prediction/').mkdir(exist_ok=True)
os.chdir('prediction/')
input_filenames = [x.split('/')[-1] for x in glob.glob(f'../{contig_dir}/*.fasta')]
contig_ids = [x.split("/")[-1].split("_")[0] for x in input_filenames] #extracting the contig ids

class TrackCPU(threading.Thread):
    '''
    Utility for tracking CPU use by creating a new thread and storing CPU use of the current process in an array.
    Returns the average usage across all CPU cores.
    '''
    def run(self):
        self.running = True
        currentProcess = psutil.Process()
        self.cpu_usage = []
        while self.running:
            self.cpu_usage.append(currentProcess.cpu_percent())

    def stop(self):
        self.running = False
        return np.mean(self.cpu_usage)/ psutil.cpu_count()
    
def prodigal():
    pathlib.Path('prodigal_output/').mkdir(exist_ok=True)
    os.chdir('prodigal_output/')
    cpu=TrackCPU()
    cpu.start()
    time=[]
    memory=[]
    start=timer()    
    with open('loginfo', 'w') as h:
        for i,id in enumerate(contig_ids):
            p=subprocess.Popen(['prodigal', '-f', 'gff', '-i', f'../../{contig_dir}/{id}_contigs.fasta', '-o', f'{id}.gff', '-a', f'{id}.faa', '-d', f'{id}.fna'], stdout=h, stderr=h)
            #subprocess.run(['psrecord', f'{p.pid}', '--plot', f'plot_{i}.png'])
            process = os.wait4(p.pid, os.WUNTRACED | os.WCONTINUED)
            time.append(process[2][0] + process[2][1])
            memory.append(process[2][3] + process[2][4])
    end=timer()
    cpu_usage=cpu.stop()
    print(time)
    print(memory)
    os.chdir('../')
    print(f"\nProdigal runtime: {round((end - start)/60, 2)} minutes \nCPU Usage {round(cpu_usage,2)} % \nTime (User+System): {humanfriendly.format_timespan(np.mean(time))} \nMemory (Shared+Unshared): {humanfriendly.format_size(np.mean(memory))}")
  
def fragGeneScan():
    pathlib.Path('fgs_output/').mkdir(exist_ok=True)
    cpu=TrackCPU()
    cpu.start()
    time=[]
    memory=[]
    start=timer()    
    with open('fgs_output/loginfo', 'w') as h:
        for id in contig_ids:   
            genome = f'-genome=../{contig_dir}/{id}_contigs.fasta'
            out = f'-out=fgs_output/{id}'
            p = subprocess.Popen(['run_FragGeneScan.pl', genome, out, '-complete=1', '-train=complete'], stdout=h, stderr=subprocess.STDOUT)
            process = os.wait4(p.pid, os.WUNTRACED | os.WCONTINUED)
            time.append(process[2][0] + process[2][1])
            memory.append(process[2][3] + process[2][4])
    end=timer()
    cpu_usage=cpu.stop()
    print(f"\nFragGeneScan runtime: {round((end - start)/60, 2)} minutes \nCPU Usage {round(cpu_usage,2)} % \nTime (User+System): {humanfriendly.format_timespan(np.mean(time))} \nMemory (Shared+Unshared): {humanfriendly.format_size(np.mean(memory))}")

def balrog():
    pathlib.Path('balrog_output').mkdir(exist_ok=True)
    os.chdir('balrog_output/')
  
    # Separate each assembly's contigs into separate FASTA files
    for file in input_filenames:
        temp_dir = f"{file.split(('.'))[0]}_split_files"
        subprocess.run(f'cp ../../{contig_dir}/{file} . ; splitfasta {file}', shell=True)
        subprocess.run(f'mv {temp_dir}/* . ; rm -r {temp_dir} {file}', shell=True)
    
    cpu=TrackCPU()
    cpu.start()
    time=[]
    memory=[]
    start=timer()  

    # Run Balrog on each individual contig
    with open('loginfo', 'w') as h:
        for file in glob.glob('*fasta'):
            p = subprocess.Popen(f'balrog -i {file} --mmseqs -o {file.split(".")[0]}.gff', stdout=h, stderr=subprocess.STDOUT, shell=True)
            process = os.wait4(p.pid, os.WUNTRACED | os.WCONTINUED)
            time.append(process[2][0] + process[2][1])
            memory.append(process[2][3] + process[2][4])
            os.remove(file)
    end=timer()
    cpu_usage=cpu.stop()
    print(f"\nBalrog runtime: {round((end - start)/60, 2)} minutes \nCPU Usage {round(cpu_usage,2)} % \nTime (User+System): {humanfriendly.format_timespan(np.mean(time))} \nMemory (Shared+Unshared): {humanfriendly.format_size(np.mean(memory))}")

    # Combine contig GFFs for each assembly and remove temporary files
    for id in contig_ids:
        temp = [x for x in glob.glob('*gff') if id in x]
        contigs = sorted(temp, key=lambda x: x.split('_')[-1].split('.')[0])
        g = open(f'{id}.gff', 'w')
        for x in contigs:
            with open(x) as h:
                g.write(h.read())
            os.remove(x)
        g.close()

    os.chdir('../')

if __name__ == "__main__":
    prodigal()
    fragGeneScan()
    balrog()
