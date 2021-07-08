#!/usr/bin/python

import random
import numpy as np
import os
import glob
import uuid
import pandas as pd
import tempfile
import argparse

def compare_coverage(bowtie,samtools,data_path,tmp_path,uid,target,sep_db,threads,util_path='utils'):
    tmp_file = tempfile.NamedTemporaryFile().name
    os.system(f'python {util_path}/compare_coverage.py --bowtie {bowtie} --samtools {samtools} --data_path {data_path} --tmp_path {tmp_path} --uid {uid} --target  {target} --target_ref_separate_db_path {sep_db} --threads {threads} > {tmp_file}')
    with open(f'{tmp_file}','r') as inpf:
        for line in inpf:
            line = line.strip()
            cov,dep = line.split(' ')
            break
    os.system(f'rm {tmp_file}')
    return [cov,dep]

def mut_seq(seq):
    mut_pos = random.sample(range(len(seq)),5)
    i = -1
    for pos in mut_pos:
        i+= 1
        r = np.random.normal(loc=0,scale=3)
        if abs(r)>i and abs(r)<i+1:
            seq = seq[:pos] + random.choice(['A','T','G','C']) + seq[(pos)+1:]
    return seq

def seq2dict(fa):
    _dict = {}
    with open(fa,'r') as inpf:
        seq = ''
        seqname = ''
        for line in inpf:
            line = line.strip()
            if line[0] == '>':
                if seqname != '':
                    _dict[seqname] = seq
                    seq = ''
                seqname = line[1:].split(' ')[0]
            else:
                seq += line
        if len(seq) > 0:
            _dict[seqname] = seq
    return _dict

def write_fq(seqs,outfile,sign):#sign: read1 or read2, avaliable: [1,2]
    if sign not in [1,2]:
        print('sign err to write fastq')
        return
    with open(outfile,'w') as outpf:
        for name,seq in seqs:
            outpf.write(f'@{name}/{sign}\n')
            outpf.write(f'{seq}\n')
            outpf.write('+\n')
            outpf.write(f"{'I'*len(seq)}\n")

def gzip(f):
    os.system(f'gzip {f}')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--bowtie',type=str,help='The bowtie2 bin path',default='bowtie2')
    parser.add_argument('--samtools',type=str,help='The samtools bin path',default='samtools')
    parser.add_argument('--background_ref_fna_path',type=str,help='The path of background genome references',default='test/refs/backgrounds/fnas/')
    parser.add_argument('--background_ref_db_path',type=str,help='The path of background genome references',default='test/refs/backgrounds/merged')

    parser.add_argument('--target_ref_fna_path',type=str,help='The path of target genome references',default='test/refs/targets/fnas/')
    parser.add_argument('--target_ref_db_path',type=str,help='The path of target genome references',default='test/refs/targets/merged')
    parser.add_argument('--target_ref_separate_db_path',type=str,help='References built from specific target strains',default='test/refs/targets/separate_dbs/')

    parser.add_argument('--background_max_num',type=int,help='Upper limit for random background genomes selection',default=100)
    parser.add_argument('--background_min_num',type=int,help='Lower limit for random background genomes selection',default=1)
    parser.add_argument('--target_max_num',type=int,help='Upper limit for random target genomes selection',default=10)
    parser.add_argument('--target_min_num',type=int,help='Lower limit for random target genomes selection',default=1)
    parser.add_argument('--data_dir',type=str,help='Dir for generated data files',default='test/data/')
    parser.add_argument('--output_dir',type=str,help='Output dir for results files',default='test/outputs/')
    parser.add_argument('--threads',type=int,help='Threads number for multiprocessing',default=4)


    args = parser.parse_args()

    curr_path = os.path.dirname(os.path.abspath(__file__))

    uid =  str(uuid.uuid4())
    bowtie = args.bowtie
    samtools = args.samtools
    background_db = args.background_ref_db_path
    d_background = args.background_ref_fna_path
    d_target = args.target_ref_fna_path
    data_path = args.data_dir
    output_path = args.output_dir
    threads = args.threads


    if not os.path.exists(output_path):
        os.mkdir(output_path)
    if not os.path.exists(data_path):
        os.mkdir(data_path)
    if not os.path.exists(f'{output_path}/info'):
        os.mkdir(f'{output_path}/info')
    if not os.path.exists(f'{output_path}/tmp'):
        os.mkdir(f'{output_path}/tmp')
    if not os.path.exists(f'{output_path}/stats'):
        os.mkdir(f'{output_path}/stats')
    if not os.path.exists(f'{output_path}/res'):
        os.mkdir(f'{output_path}/res')
    if not os.path.exists(f'{data_path}/unmapped'):
        os.mkdir(f'{data_path}/unmapped')

    
    if args.target_min_num > args.target_max_num:
        print('target reference number error!')
        exit(0)
    if args.background_min_num > args.background_max_num:
        print('background reference number error!')
        exit(0)

    target_no = random.choice(range(args.target_min_num,args.target_max_num))
    background_no = random.choice(range(args.background_min_num,args.background_max_num))
    total_no = target_no + background_no

    print(f'target_no: {target_no}, background_no: {background_no}')

    pool = []
    frac = []
    for i in range(total_no):
        r = round(abs(np.random.normal(loc=0,scale=5)))
        if r==0:
            r+=1
        pool += [i]*r
        frac += [r]

    target_fs = random.sample(glob.glob(f'{d_target}/*.fna'),target_no)
    background_fs = random.sample(glob.glob(f'{d_background}/*.fna'),background_no)

    basenames = [os.path.basename(x).replace('.fna','') for x in target_fs+background_fs]
    genomes = [seq2dict(x) for x in target_fs+background_fs]
    count = (random.choice(range(10)) + 1)*100000
    print('count: ',count)

    r1_lst = []
    r2_lst = []
    for i in range(count):

        g_i = random.choice(pool)
        genome_dict = genomes[g_i]
        contigs = list(genome_dict.keys())
        contig = random.choice(contigs)
        seq = genome_dict[contig]

        start = random.choice(range(len(seq)))

        if start + 100 > len(seq):
            start = max(0,len(seq) - 100 - 1)

        end = start + 100 + random.choice(range(50,200)) + 100

        if end >= len(seq):
            end = len(seq)-1

        r1 = seq[start:start+100]
        r2 = seq[max(0,end-100):end]
        if len(r1)<100 or len(r2)<100:
            continue

        r1_lst.append((f'{basenames[g_i]}-{i}',r1))
        r2_lst.append((f'{basenames[g_i]}-{i}',r2))


    write_fq(r1_lst,f'{data_path}/{uid}_1.fastq',1)
    write_fq(r2_lst,f'{data_path}/{uid}_2.fastq',2)
    gzip(f'{data_path}/{uid}_1.fastq')
    gzip(f'{data_path}/{uid}_2.fastq')

    with open(f'{output_path}/info/{uid}.info','w') as outpf:
        outpf.write(f'#{count}\n')
        outpf.write(f'#{target_no}\t{background_no}\n')
        frac_str = '\t'.join([str(x) for x in frac])
        outpf.write(f"#{frac_str}\n")
        i = -1
        for bn in basenames:
            i += 1
            if i + 1 <= target_no:
                outpf.write(f'{bn}\t{frac[i]}\ttarget\n')
            else:
                outpf.write(f'{bn}\t{frac[i]}\tbackground\n')

    os.system(f'time {bowtie}  -x {background_db} -1 {data_path}/{uid}_1.fastq.gz -2 {data_path}/{uid}_2.fastq.gz -S {output_path}/tmp/{uid}.sam --no-unal -p 20 --un-conc-gz {data_path}/unmapped/{uid}_%.fastq.gz')
    os.system(f'rm {output_path}/tmp/{uid}.sam')
    os.system(f'python {curr_path}/profiling.py --bowtie {bowtie} --samtools {samtools} --fq1 {data_path}/unmapped/{uid}_1.fastq.gz --outbase {output_path} --ref_fna_dir {d_target} --ref_db_all {args.target_ref_db_path} --ref_db_separate_path {args.target_ref_separate_db_path} --threads {args.threads}')
    
    df = pd.read_csv(f'{output_path}/stats/stat.{uid}.csv')

    cov_dep_lst = []
    i = -1
    res = pd.DataFrame(columns=['name','real_cov','real_dep','frac','rank','predicted_cov','predicted_dep','best_mean_nm'])
    target_bnames = []
    print('frac: ',frac)
    for f in target_fs:
        i+=1
        bname = os.path.basename(f).replace('.fna','')
        target_bnames.append(bname)
        item = []
        item += [bname] + compare_coverage(bowtie,samtools,data_path,output_path+'/tmp',uid,bname,args.target_ref_separate_db_path,threads,curr_path+'/utils') + [frac[i]/sum(frac)]  
        row = df[df['src'] == f'{uid}-{bname}']
        rank = row.index[0]
        predicted_cov = list(row['sum-cov'])[0]
        predicted_dep = list(row['mean_dep'])[0]
        nm = list(row['nm-best'])[0].split('-')[2]  
        item += [rank,predicted_cov,predicted_dep,nm]
        cov_dep_lst.append(item)
        res.loc[len(res)] = item
    for f in glob.glob(f'{d_target}/*.fna'):
        bname = os.path.basename(f).replace('.fna','')
        if bname in target_bnames:
            continue
        item = []
        item += [bname] + [0.0, 0.0, 0.0] 
        row = df[df['src'] == f'{uid}-{bname}']
        rank = row.index[0]
        predicted_cov = list(row['sum-cov'])[0]
        predicted_dep = list(row['mean_dep'])[0]
        nm = list(row['nm-best'])[0].split('-')[2]  
        item += [rank,predicted_cov,predicted_dep,nm]
        cov_dep_lst.append(item)
        res.loc[len(res)] = item




    res.to_csv(f'{output_path}/res/{uid}.csv')
