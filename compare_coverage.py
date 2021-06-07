#!/usr/bin/python
import sys
import os
import pandas as pd
import random
import gzip
import glob
import tempfile

def filter_fq(in_fq,out_fq,target):
    with open(out_fq,'w') as outpf:
        with gzip.open(in_fq,'r') as inpf:
             run = False
             for line in inpf:
                 line = line.decode('utf-8')
                 if line[0] == '@':
                     if target in line:
                         run = True
                     else:
                         run = False
                 if run:
                     outpf.write(line)

def get_coverage(bam_dir,name):

    tmp_file = tempfile.NamedTemporaryFile().name
    os.system(f'samtools coverage {bam_dir}/{name}.sort.bam > {tmp_file}')
    df = pd.read_csv(tmp_file,sep='\t')
    df['src'] = name
    os.system(f'rm {tmp_file}')
    return df

if __name__ == '__main__':

    uid = sys.argv[1]
    target = sys.argv[2]

    filter_fq(f'./data/{uid}_1.fastq.gz',f'tmp/{uid}_1.fastq',target)
    filter_fq(f'./data/{uid}_2.fastq.gz',f'tmp/{uid}_2.fastq',target)

    #map
    ref = glob.glob(f'../refs/136Fstrains/separate_dbs/{target}*.1.bt2')[0].replace('.1.bt2','')
    os.system(f'bowtie2 -q --quiet -1 tmp/{uid}_1.fastq -2 tmp/{uid}_2.fastq -x {ref} -S tmp/{uid}.sam -p 10 --very-sensitive --no-unal')
    os.system(f'samtools view -bS tmp/{uid}.sam -F 256 -q 1 > tmp/{uid}.bam')
    os.system(f'samtools sort  tmp/{uid}.sam -o tmp/{uid}.sort.bam')
    os.system(f'samtools index  tmp/{uid}.sort.bam')

    df = get_coverage('tmp',uid)
    print(((df['coverage']*df['endpos'])/100).sum()/df['endpos'].sum(),(df.meandepth*df.covbases).sum()/df.covbases.sum()   )

 

