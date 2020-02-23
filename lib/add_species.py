# -*- coding: utf-8 -*-
# 2019.12.24
# @zlk
# 在已有结果中添加一个物种
from collections import OrderedDict
import argparse
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import gene_fa
import multiprocessing
import tandem5
import glob
import flanking_gene_total9
import synteny6


def pick_file(path, end_str):
    files1 = os.listdir(path)
    for f1 in files1:
        if f1.endswith(end_str):
            file_path = path + f1
    return file_path


def add_species(fa,gff,threads_num,output_file,file_num,tools):
    fa_file_list = fa.split(',')
    gff_file_list = gff.split(',')
    # 提取蛋白序列
    pool = multiprocessing.Pool(threads_num)
    for fa_index in range(len(fa_file_list)):
        result = pool.apply_async(gene_fa.gene_fa, args=(
        gff_file_list[fa_index], fa_file_list[fa_index], fa_index + file_num, output_file))
    result.get()
    pool.close()
    pool.join()
    tandem5.self_blast(len(fa_file_list)+file_num, output_file, threads_num, tools, list(range(file_num)))
    # 挑选tandem，并去除在基因序列文件中去除tandem
    tandem5.run_pick_tandem(output_file, threads_num, len(fa_file_list)+file_num)
    # 两两物种进行blast
    file_list=[]
    for file_index in range(len(fa_file_list)+file_num):
        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        for f in files:
            if re.search(r'file\d*_file\d', f):
                file_list.append(f)
    flanking_gene_total9.resp_blast(len(fa_file_list)+file_num, output_file, tools, threads_num,file_list)
    # 对两两的blast结果进行共线性分析
    flanking_gene_total9.run_flanking(threads_num, output_file, len(fa_file_list)+file_num)
    # 进行多基因组共线性分析
    synteny6.synteny(len(fa_file_list)+file_num)