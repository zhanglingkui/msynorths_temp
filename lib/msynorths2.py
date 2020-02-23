# -*- coding: utf-8 -*-
# 2019.08.08
# @zlk
# 进行多个物种的共线性分析
# 将最后共线性文件添加模块度文件
from collections import OrderedDict
import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import gene_fa
import multiprocessing
import tandem5
import glob
import flanking_gene_total9
import synteny6

def msynorths(fa,gff,threads_num,output_file,tools):
    fa_file_list = fa.split(',')
    gff_file_list = gff.split(',')
    # 提取蛋白序列
    gene_fa.run_gene_fa(fa_file_list, gff_file_list, output_file,threads_num)
    # 自身blast，用作去除tandem
    tandem5.self_blast(len(fa_file_list), output_file, threads_num,tools,[])

    # 挑选tandem，并去除在基因序列文件中去除tandem
    tandem5.run_pick_tandem(output_file, threads_num, len(fa_file_list))
    # 两两物种进行blast
    flanking_gene_total9.resp_blast(len(fa_file_list), output_file, tools, threads_num,[])
    # 对两两的blast结果进行共线性分析
    flanking_gene_total9.run_flanking(threads_num, output_file, len(fa_file_list))
    # 进行多基因组共线性分析
    synteny6.synteny(len(fa_file_list))
