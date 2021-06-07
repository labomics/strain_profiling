#!/usr/bin/python

import random
import numpy as np
import os
import glob
import uuid
import pandas as pd
import tempfile

def compare_coverage(uid,target):
    tmp_file = tempfile.NamedTemporaryFile().name
    os.system(f'python compare_coverage.py {uid} {target} > {tmp_file}')
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

    uid =  str(uuid.uuid4())
    d_background = '../refs/oldAllRef_without_GCA000166035/fnas/'
    d_target = '../refs/136Fstrains/fnas/'


    target_no = random.choice(range(1,11))
    background_no = random.choice(range(1,101))
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


    write_fq(r1_lst,f'data/{uid}_1.fastq',1)
    write_fq(r2_lst,f'data/{uid}_2.fastq',2)
    gzip(f'data/{uid}_1.fastq')
    gzip(f'data/{uid}_2.fastq')

    with open(f'info/{uid}.info','w') as outpf:
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

    os.system(f'time bowtie2  -x ../refs/oldAllRef_without_GCA000166035/merged -1 data/{uid}_1.fastq.gz -2 data/{uid}_2.fastq.gz -S tmp/{uid}.sam --no-unal -p 20 --un-conc-gz unmapped/{uid}_%.fastq.gz')
    os.system(f'rm tmp/{uid}.sam')
    os.system(f'python profiling.py unmapped/{uid}_1.fastq.gz .')
    
    df = pd.read_csv(f'stats/stat.{uid}.csv')

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
        item += [bname] + compare_coverage(uid,bname) + [frac[i]/sum(frac)]  
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




    res.to_csv(f'res/{uid}.csv')
    os.system(f'rm output*/{uid}*')
    os.system(f'rm output*/{uid}*')
