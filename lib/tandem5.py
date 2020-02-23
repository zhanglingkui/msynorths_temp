# -*- coding: utf-8 -*-
# 2019.10.25
# @zlk
# 挑选出tandem，并将其在蛋白质序列中删除，写成*_notandem.fa文件
# 原有文件没有将删除的tandem 形成一个文件做后续分析
# 原有挑选tandem的原则是一串tandem累加后计算容忍基因数量，现改为两两之间计算。

from collections import OrderedDict
from Bio import SeqIO
import multiprocessing
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import sys
import os


def pick_tandem(gff_file, prot_file, blast_file):
    input_gff_file = open(gff_file, 'r')
    input_blast_file = open(blast_file, 'r')
    input_prot_file = prot_file
    out_prot_file = open(input_prot_file + '_notandem', 'w')
    out_tandem_array = open(input_prot_file + 'tandem_array', 'w')

    input_nohom_gene_num = 4
    input_identity = 30
    input_coverage = 0.2
    # 获得gene的位置顺序，写入字典
    gene_index_length_dic = {}
    index = 1
    for line in input_gff_file:
        line_list = line.strip().split('\t')
        gene_index_length_dic[line_list[0]] = []
        # 将基因的顺序添加进来
        gene_index_length_dic[line_list[0]].append(index)
        index += 1
        # 将cds长度添加进来
        gene_index_length_dic[line_list[0]].append(line_list[-1].strip())
    # 将同一个基因的同源基因写入同一个字典
    gene_blast_dic = OrderedDict()
    for line in input_blast_file:
        line_list = line.split('\t')
        # 设置identity最低值
        if float(line_list[2]) < input_identity:
            continue
        # 设置coverage

        elif int(line_list[3].strip()) / int(gene_index_length_dic[line_list[0]][1]) < input_coverage:
            continue
        elif line_list[0] in gene_blast_dic.keys():
            # 将每一个基因的同源基因ID和他的index合成一个列表写进字典
            gene_index_list = [line_list[1], gene_index_length_dic[line_list[1]][0]]
            if gene_index_list in gene_blast_dic[line_list[0]]:
                continue
            else:

                gene_blast_dic[line_list[0]].append(gene_index_list)
        else:
            gene_index_list = [line_list[1], gene_index_length_dic[line_list[1]][0]]
            gene_blast_dic[line_list[0]] = []
            gene_blast_dic[line_list[0]].append(gene_index_list)

    del_gene_dict = {}
    for key, value in gene_blast_dic.items():
        this_gene_index = int(gene_index_length_dic[key][0])
        # 将list按第二列的index排序
        value_sort = sorted(value, key=(lambda s: s[1]))
        # 遍历这个value，看其是否有tandem
        for blast_list in value_sort:
            # index小于等于query gene的直接跳过
            if blast_list[1] <= this_gene_index:
                continue
            # 当同源基因内包含的非同源基因数大于等于五就停止
            # 设置非同源基因容忍度
            elif (blast_list[1] - this_gene_index) < input_nohom_gene_num:
                del_gene_dict[key] = blast_list[0]
                break
    # 将一对一对的tandem根据相同的ID进行连接起来
    used_list = []
    for key, value in del_gene_dict.items():
        if key in used_list:
            continue
        value = [value]
        for i in value:

            if i in del_gene_dict.keys():
                value.append(del_gene_dict[i])
                used_list.append(i)
        del_gene_dict[key] = value
    # 遍历tandem字典，将其写入tandem_array文件，并且将需要删除的ID写入列表
    del_only_gene_list = []
    for key, value in del_gene_dict.items():
        if key in used_list:
            continue
        out_tandem_array.write(key + '\t')
        for i in value:
            out_tandem_array.write(i + '\t')
            del_only_gene_list.append(i)
        out_tandem_array.write('\n')

    prot_dict = SeqIO.to_dict(SeqIO.parse(input_prot_file, "fasta"))
    for key in prot_dict.keys():
        if key in del_only_gene_list:
            continue
        else:
            out_prot_file.write('>' + key + '\n' + str(prot_dict[key].seq) + '\n')
def self_blast(loop_num,output_file,threads_num,tools,do_list):
    for file_index in range(loop_num):
        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        if file_index in do_list:
            continue
        for f in files:
            if f.endswith('_prot'):
                prot_file = work_path + f
        if tools == 'diamond':
            commond11 = 'mkdir ' + work_path + 'db/'
            os.system(commond11)
            commond1 = 'diamond makedb --in ' + prot_file + ' -d ' + work_path + 'db/ref'
            os.system(commond1)
            commond2 = 'diamond blastp -d ' + work_path + 'db/ref -q ' + prot_file + ' --quiet --sensitive -e 1e-5  -p ' + str(
                threads_num) + ' -o ' + work_path + 'self_blast'
            os.system(commond2)
        else:
            # 建库
            commond1 = 'makeblastdb -in ' + prot_file + ' -dbtype prot -parse_seqids -out ' + work_path + 'db/ref'
            os.system(commond1)
            # 比对
            commond2 = 'blastp -db ' + work_path + 'db/ref -query ' + prot_file + '  -num_threads ' + str(
                threads_num) + ' -outfmt 6 -evalue 1e-5 -out ' + work_path + 'self_blast'
            os.system(commond2)
def run_pick_tandem(output_file,threads_num,loop_num):
    pool = multiprocessing.Pool(threads_num)
    for file_index in range(loop_num):
        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        for f in files:
            if f.endswith('_prot'):
                prot_file = work_path + f
            if f.endswith('_longest.gff'):
                gff_file = work_path + f
        result = pool.apply_async(pick_tandem, (gff_file, prot_file, work_path + 'self_blast'))
    result.get()
    pool.close()
    pool.join()
if __name__ == "__main__":
    gff_file = sys.argv[1]
    prot_file = sys.argv[2]
    blast_file = sys.argv[3]
    pick_tandem(gff_file, prot_file, blast_file)
