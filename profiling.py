#!/usr/bin/python
import sys
import os
import glob
from multiprocessing import Pool
import tempfile
import pandas as pd
import pysam
import numpy as np
import re
import time
import datetime

def mapping(fq1,fq2,name,ref,output,minq=0,nthreads=10):
    t1 = time.time()
    if (os.path.exists(f'{output}/{name}.bam') or os.path.exists(f'{output}/{name}.sort.bam')):
        return
    os.system(f'bowtie2 -q --quiet -1 {fq1} -2 {fq2} -x {ref}  -S {output}/{name}.sam -p {nthreads} --very-sensitive --no-unal')
    t2 = time.time()
    os.system(f'samtools view -bS {output}/{name}.sam  -F 256 -q {minq} > {output}/{name}.bam')
    os.system(f'rm {output}/{name}.sam')

    os.system(f'samtools sort {output}/{name}.bam -o {output}/{name}.sort.bam')
    os.system(f'samtools index {output}/{name}.sort.bam')
    os.system(f'rm {output}/{name}.bam')
    t3 = time.time()
    print(f'{name}, mapping time:{(t2 - t1)/1000}, other time:{(t3-t2)/1000}') 
    return

def get_coverage(bam_dir,name):

    if not os.path.exists(f'{bam_dir}/{name}.sort.bam'):
        os.system(f'samtools sort {bam_dir}/{name}.bam -o {bam_dir}/{name}.sort.bam')
        os.system(f'samtools index {bam_dir}/{name}.sort.bam')

    tmp_file = tempfile.NamedTemporaryFile().name
    os.system(f'samtools coverage {bam_dir}/{name}.sort.bam > {tmp_file}')
    df = pd.read_csv(tmp_file,sep='\t')
    df['src'] = name
    os.system(f'rm {tmp_file}')
    return df

def get_intsect(reg1,reg2):
    l1,r1 = reg1
    l2,r2 = reg2

    if max(l1,l2) < min(r1,r2):
        return max(l1,l2),min(r1,r2)
    return -1,-1

def get_partial_nm(sub_reg,reg,md):
    md = md.replace('^','')
    sub_l,sub_r = sub_reg
    l,r = reg

    seq = ''
    for item in re.findall('\d+|[A-Z]',md):
        if item == '0':
            continue
        if item.isnumeric():
            seq += '*'*int(item)
        else:
            seq += item

    return sum([x!='*' for x in seq[sub_l-l:sub_r-l+1]])

def merge(items):

    if len(items) <= 1:
        return

    tail = items[-1]
    t_key,t_start,t_end,t_cigar,t_nm,t_md = tail

    pre_tail = items[-2]
    p_key,p_start,p_end,p_cigar,p_nm,p_md = pre_tail

    if t_start >= p_end:
        return 

    isl,isr = get_intsect([p_start,p_end],[t_start,t_end])
    assert isl!=-1,'intersect exception'

    t_partial_nm = get_partial_nm([isl,isr],[t_start,t_end],t_md)
    p_partial_nm = get_partial_nm([isl,isr],[p_start,p_end],p_md)

    if t_partial_nm != p_partial_nm:
        if t_nm > p_nm:
            items.pop()
            merge(items)
    return 


if __name__ == '__main__':

    fq1 = sys.argv[1]

    if len(sys.argv) == 3:
        out_base = sys.argv[2]
    else:
        out_base = '.'

    print('out_base:',out_base)

    fq2 = fq1.replace('_1.fastq.gz','_2.fastq.gz').replace('.unmapped.1.gz','.unmapped.2.gz')

    ref_dir = '/data1/chenyaowen/workspace/metaSNP_finalversion/refs/136Fstrains/'
    overall_ref = f'{ref_dir}/merged'
    separate_ref = f'{ref_dir}/separate_dbs/'

    name = os.path.basename(fq1).replace('_1.fastq.gz','').replace('.unmapped.1.gz','')
    print(name)

    if os.path.exists(f'{out_base}/stats/stat.{name}.csv'):
        exit(0)
        
    output_overall = f'{out_base}/output_overall'
    output_sep = f'{out_base}/output_sep'

    if not os.path.exists(output_overall):
        os.system(f'mkdir {output_overall}')
    if not os.path.exists(output_sep):
        os.system(f'mkdir {output_sep}')

    if not os.path.exists(f'{out_base}/stats'):
        os.system(f'mkdir {out_base}/stats')
    if not os.path.exists(f'{out_base}/coverage'):
        os.system(f'mkdir {out_base}/coverage')


    nthreads = 10
    print('#step1, mapping to overall_ref')
    st_time = time.time()
    print('time: ',datetime.datetime.now().ctime())
    mapping(fq1,fq2,name,overall_ref,output_overall,minq=0,nthreads=nthreads)#no need to keep unique mapped reads
    print('time: ',st_time - time.time())
    print('time: ',datetime.datetime.now().ctime())

    print('#step2, mapping to separate_ref')
    print('time: ',datetime.datetime.now().ctime())
    st_time = time.time()
    refs = []
    with Pool(nthreads) as pool:
        for f in glob.glob(f'{ref_dir}/fnas/*.fna'):
            ref = os.path.basename(f).replace('.fna','')
            refs.append(ref)
            pool.apply_async(mapping,(fq1,fq2,f'{name}-{ref}',f'{separate_ref}/{ref}',output_sep,0,8))#to use 0 minq, ignore intra-genome similarity

        pool.close()
        pool.join()

    print('time: ',st_time - time.time())
    print('time: ',datetime.datetime.now().ctime())

    print('#step4, reads type assignment and nm distribution')
    print('#step4, keep the best reads for same region')

    st_time = time.time()
    print('time: ',datetime.datetime.now().ctime())
    save = pysam.set_verbosity(0)


    overall_reads2contig = {}
    bamfile = pysam.AlignmentFile(f'{out_base}/output_overall/{name}.sort.bam','rb')
    pysam.set_verbosity(save)

    for aln in bamfile:
        if aln.is_read1:
            key = aln.query_name + '_1'
        else:
            key = aln.query_name + '_2'
        nm = int(dict(aln.tags)['NM'])
        if nm/aln.query_length <= 0.05:
            overall_reads2contig[key] = (aln.reference_name,nm,aln.mapping_quality)
    bamfile.close()

    def classfy_reads(ref):
        total,best,consist,falsepositive,equal = [0]*5
        nms = []
        best_nms = []
        consist_nms = []
        seqs_lst = {}
        sep_bam = f'{output_sep}/{name}-{ref}.sort.bam'
        if not os.path.exists(sep_bam):
            return None
        bamfile = pysam.AlignmentFile(sep_bam,"rb")
        pysam.set_verbosity(save)
        output_bam_f = f'{output_sep}/{name}-{ref}.best.bam'
        if True:
            output_bam = pysam.AlignmentFile(output_bam_f, "wb", template=bamfile)

            aln_dict = {}


            for aln in bamfile:

                target = aln.reference_name
                if target == '' or target == '*':
                    continue

                tags = dict(aln.tags)
                ref_start = aln.reference_start
                ref_end = aln.reference_end
                cigar = aln.cigarstring
                nm = int(tags['NM'])
                md = tags['MD']

                if nm/aln.query_length > 0.05:
                    continue

                nms.append(nm)
                total += 1

                if aln.is_read1:
                    key = aln.query_name + '_1'
                else:
                    key = aln.query_name + '_2'
                
                is_fp = False
                is_equal = False
                if key not in overall_reads2contig:
                    continue #abnormal
                if overall_reads2contig[key][0] == target:
                    if overall_reads2contig[key][2] <= 10:
                        equal += 1
                        is_equal = True
                    else:
                        consist += 1
                        consist_nms.append(nm)
                else: # map to else-where
                    #1. same similarity, this could be called equal..
                    if overall_reads2contig[key][1] >= nm:
                        equal += 1
                        is_equal = True
                    else:
                        falsepositive += 1
                        is_fp = True
                if is_fp:
                    continue
                record = [key,ref_start,ref_end,cigar,nm,md]
                if target not in seqs_lst:
                    seqs_lst[target] = []
                seqs_lst[target].append(record)
                try:
                    merge(seqs_lst[target])
                except Exception as e:
                    print('err')
                    print(e)
                if seqs_lst[target][-1] == record:
                    aln_dict[key] = aln

            with open(f'{output_sep}/{name}-{ref}.bestmapped.reads.txt','w') as outpf:
                for contig in seqs_lst:
                    for item in seqs_lst[contig]:
                        best_nms.append(item[-2])
                        best += 1
                        output_bam.write(aln_dict[item[0]])
                        outpf.write(f'{item[0]}\n')

            output_bam.close()
            bamfile.close()

        best_df = get_coverage(output_sep,f'{name}-{ref}.best')
        best_df['width'] = (best_df['coverage']/100)*best_df['endpos']
        sum_cov = (best_df['width'].sum()/best_df['endpos'].sum())*100
        mean_dep = (best_df.meandepth*best_df.covbases).sum()/best_df.covbases.sum()

        nms = np.array(nms)
        best_nms = np.array(best_nms)
        consist_nms = np.array(consist_nms)

        if len(nms) == 0:
            nm_total = np.nan
        else:
            nm_total = f'{nms.max()}-{nms.min()}-{nms.mean()}-{np.median(nms)}-{np.std(nms)}'

        if len(best_nms) == 0:
            nm_best = np.nan
        else:
            nm_best = f'{best_nms.max()}-{best_nms.min()}-{best_nms.mean()}-{np.median(best_nms)}-{np.std(best_nms)}'

        if len(consist_nms) == 0:
            nm_consist = np.nan
        else:
            nm_consist = f'{consist_nms.max()}-{consist_nms.min()}-{consist_nms.mean()}-{np.median(consist_nms)}-{np.std(consist_nms)}'


        return [f'{name}-{ref}',total,consist,best,falsepositive,equal,sum_cov,mean_dep,nm_total,nm_best,nm_consist]

    with Pool(nthreads) as pool:
        arr = pool.map(classfy_reads,refs)

        new_arr = []
        for item in arr:
            if item is not None:
                new_arr.append(item)
        arr = new_arr

        stat_df = pd.DataFrame(arr)
        stat_df.columns = ['src','total','consist','best','falsepositive','equal','sum-cov','mean_dep','nm-total','nm-best','nm-consist']
        stat_df = stat_df.sort_values(by='sum-cov',ascending=False).reset_index(drop=True)
        stat_df.to_csv(f'{out_base}/stats/stat.{name}.csv',index=None)

        best_reads = []
        for bf in glob.glob(f'{output_sep}/{name}-*.bestmapped.reads.txt'):
            with open(bf,'r') as inpf:
                for line in inpf:
                    best_reads.append(line.strip())
        best_reads = len(set(best_reads))

        with open(f'{out_base}/stats/best_number.{name}.txt','w') as outpf:
            outpf.write(f'{best_reads}\n')

    print('time: ',st_time - time.time())
    print('time: ',datetime.datetime.now().ctime())

    #os.system(f'rm {output_overall}/{name}.*')
    #os.system(f'rm {output_sep}/{name}*genomic.sort.bam')
    #os.system(f'rm {output_sep}/{name}*genomic.sort.bam.bai')
    #os.system(f'rm {output_sep}/{name}*.best.bam')
    #os.system(f'rm {output_sep}/{name}*bestmapped.reads.txt')


