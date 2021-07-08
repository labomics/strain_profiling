#!/usr/bin/python
import sys
import os
import pandas as pd
import random
import gzip
import glob
import tempfile
import argparse

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

def get_coverage(bam_dir,name,samtools='samtools'):

    tmp_file = tempfile.NamedTemporaryFile().name
    os.system(f'{samtools} coverage {bam_dir}/{name}.sort.bam > {tmp_file}')
    df = pd.read_csv(tmp_file,sep='\t')
    df['src'] = name
    os.system(f'rm {tmp_file}')
    return df

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--bowtie',type=str,help='The bowtie2 bin path',default='bowtie2')
    parser.add_argument('--samtools',type=str,help='The samtools bin path',default='samtools')
    parser.add_argument('--data_path',type=str,help='The data file dir',default='test/data/')
    parser.add_argument('--tmp_path',type=str,help='The tmp file dir',default='test/outputs/tmp/')
    parser.add_argument('--uid',type=str,help='The file uid',default='')
    parser.add_argument('--target',type=str,help='The target name',default='')
    parser.add_argument('--target_ref_separate_db_path',type=str,help='References built from specific target strains',default='test/refs/targets/separate_dbs/')
    parser.add_argument('--threads',type=int,help='Threads number for multiprocessing',default=4)
    
    args = parser.parse_args()

    uid = args.uid
    target = args.target
    sep_db = args.target_ref_separate_db_path
    threads = args.threads
    tmp_path = args.tmp_path
    data_path = args.data_path
    bowtie = args.bowtie
    samtools = args.samtools

    filter_fq(f'{data_path}/{uid}_1.fastq.gz',f'{tmp_path}/{uid}_1.fastq',target)
    filter_fq(f'{data_path}/{uid}_2.fastq.gz',f'{tmp_path}/{uid}_2.fastq',target)

    #map
    ref = glob.glob(f'{sep_db}/{target}*.1.bt2')[0].replace('.1.bt2','')
    os.system(f'{bowtie} -q --quiet -1 {tmp_path}/{uid}_1.fastq -2 {tmp_path}/{uid}_2.fastq -x {ref} -S {tmp_path}/{uid}.sam -p {threads} --very-sensitive --no-unal')
    os.system(f'{samtools} view -bS {tmp_path}/{uid}.sam -F 256 -q 1 > {tmp_path}/{uid}.bam')
    os.system(f'{samtools} sort  {tmp_path}/{uid}.sam -o {tmp_path}/{uid}.sort.bam')
    os.system(f'{samtools} index  {tmp_path}/{uid}.sort.bam')

    df = get_coverage(tmp_path,uid,samtools)
    print(((df['coverage']*df['endpos'])/100).sum()/df['endpos'].sum(),(df.meandepth*df.covbases).sum()/df.covbases.sum()   )

 

