#!/usr/bin/python

import sys
import os
import pandas as pd
import glob
import argparse


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--bowtie',type=str,help='The bowtie2 bin path',default='bowtie2')
    parser.add_argument('--data_path',type=str,help='The fastq data path',default='test/data/')
    parser.add_argument('--output_path',type=str,help='The output path',default='test/outputs/')
    parser.add_argument('--background_ref_db_all',type=str,help='References built from all background genomes',default='test/refs/backgrounds/merged')

    args = parser.parse_args()

    bowtie = args.bowtie
    data_path = args.data_path
    output_path = args.output_path
    background_ref = args.ref_db_all


    if not os.path.exists(data_path):
        print('data path not exists')
        exit(0)
    if not os.path.exists(data_path+'/unmapped'):
        os.system(f'mkdir {data_path}/unmapped')
    if not os.path.exists(output_path):
        os.system(f'mkdir {output_path}')
    if not os.path.exists(output_path+'/tmp'):
        os.system(f'mkdir {output_path}/tmp')

    for f in glob.glob(f'{data_path}/*_1P.fq.gz'):
        smp = os.path.basename(f).replace('_1P.fq.gz','')
        print(smp)
        data1 =  f
        data2 = f.replace('_1P','_2P')
        if os.path.exists(data1) and os.path.exists(data2):
            print(f'{smp} bowtie on whole genome')
            cmd = f'time {bowtie} -x {background_ref} -1 {data1} -2 {data2}  -S {output_path}/tmp/{smp}.sam --no-unal -p {threads} --un-conc-gz {data_path}/unmapped/{smp}.unmapped.gz'
            os.system(cmd)
            os.system(f'rm {output_path}/tmp/{smp}.sam')

