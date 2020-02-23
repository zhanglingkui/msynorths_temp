# -*- coding: utf-8 -*-
# 2019.11.19
# @zlk
# 2019.11.29 添加一个排序，将基因按照染色体，和物理位置进行排序。
'''
输入gff文件，fasta文件，输出文件，得到一个新的文件类型*_new.gff，这里面包含有基因和mrna的id 以及cds序列的位置、个数、长度
然后得到一个*_longest.gff，这个是经过去除可变剪辑的文件，然后根据这个文件去fa文件里提取相应的序列，然后进行翻译
'''
import argparse
from collections import OrderedDict
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re
import sys
import multiprocessing

# 将整个文件封装成一个函数
def gene_fa(gff_file, fa_file, file_index,output_file):
    file_name = output_file+'/file' + str(file_index) + '/'
    commond_mkdir='mkdir '+ file_name
    os.system(commond_mkdir)
    gff_title = gff_file.split('/')[-1]
    fa_title = fa_file.split('/')[-1]
    newFile = open(file_name + gff_title + '_new.gff', 'w')
    longest_mrna_file = open(file_name + gff_title + '_longest.gff', 'w')
    output_pro_file = open(file_name + fa_title + '_prot', 'w')
    output_chr_file = open(file_name + 'chr_length.txt', 'w')

    def get_chr_length(file_name,out_chr_file):
        fa_dict = SeqIO.to_dict(SeqIO.parse(file_name, "fasta", IUPAC.unambiguous_dna))
        for i in fa_dict.keys():
            chr_length = i + '\t' + str(len(fa_dict[i]))+'\n'
            out_chr_file.write(chr_length)

    def sort_file_func(path):
        input_file = open(path, 'r')
        chr_dict = OrderedDict()
        for line in input_file:
            line_list = line.strip().split('\t')
            line3_list = line_list[3].split(':')
            if line3_list[0] in chr_dict.keys():
                chr_dict[line3_list[0]].append([line3_list[1], line])
            else:
                chr_dict[line3_list[0]] = []
                chr_dict[line3_list[0]].append([line3_list[1], line])
        input_file.close()
        output_file = open(path, 'w')
        for key, value in chr_dict.items():
            value.sort(key=lambda x: int(x[0]))
            for line in value:
                output_file.write(line[1])
        output_file.close()

    # 在gff文件第九列提取出其ID
    def pick_id(list):
        list9 = list.split(';')
        # 遍历第九列
        for m in range(len(list9)):
            if list9[m].startswith('ID'):
                list_id = list9[m].strip()
            else:
                continue
            return list_id

    def pick_parent(list):
        list9 = list.split(';')
        # 遍历第九列
        for m in range(len(list9)):
            if list9[m].startswith('Parent'):
                list_parent = list9[m].strip()
            else:
                continue
            return list_parent

    # 判断调整第九列的
    def adjustCol9(list):
        # 将第一行的第九列保存
        list_1 = list[0][8]

        # 遍历genelist
        for i in range(len(list)):
            # 将每一行的第九列以；切分成列表
            # list_9 = (list[i][8]).split(';')
            # 对第一行操作将其只剩下ID列
            if list[0][2] == 'gene' and i == 0:
                list[i][8] = pick_id(list_1)
                list_mrna = list[i][8]
            elif list[0][2] == "mRNA" and i == 0:
                list[i][8] = pick_id(list_1)
                list_mrna = list[i][8]
            elif list[i][2] == 'mRNA' and i > 0:
                list_mrna = list[i][8]
                if list[0][2] == 'mRNA':
                    list[i][8] = pick_id(list_mrna)
                else:
                    list[i][8] = pick_id(list_mrna) + ';Parent' + pick_id(list_1)[2:]
            else:

                list[i][8] = 'Parent' + pick_id(list_mrna)[2:]

    # 将列表写入

    def newGffWrite(filename, list):
        for part in list:
            filename.write(str(part) + '\t')
        filename.write('\n')

    # 生成新的gff文件类型

    def new_gff(list):
        # 生成新的gff文件
        cds_num = 4
        cds_num2 = 1
        cds_total_length = 0
        listNew = ['gene1', 'mrna1', 'strand', 'mrna_position', 'num(cds)', 'length']
        for j in range(len(list)):
            # 如果基因列表的第一行的第三列是gene就在新列表第一列写入'gene'id
            if list[j][2] == 'gene' and j == 0:
                listNew[0] = list[0][8]
                listNew[2] = list[0][6]
                if file_type == 1:
                    listNew[1] = '.'
                else:
                    continue
            # 如果基因列表的第一行的第二列为mrna
            elif list[j][2] == 'mRNA' and j == 0:
                listNew[0] = list[0][8]
                listNew[1] = list[0][8]
                listNew[2] = list[0][6]
                listNew[3] = list[j][0] + ':' + list[j][3] + ':' + list[j][4]
            # 如果第二行为mrna就将它的ID写入第二列
            elif list[j][2] == 'mRNA' and j == 1:

                mrna_id = list[1][8]
                listNew[1] = pick_id(mrna_id)
                listNew[3] = list[j][0] + ':' + list[j][3] + ':' + list[j][4]
            elif list[j][2] == 'CDS':
                # 将cds序列写入列表
                cds_location = list[j][0] + ':' + list[j][3] + ':' + list[j][4]
                listNew.insert(cds_num, cds_location)
                cds_num += 1
                # 写入cds数量
                listNew[-2] = cds_num2
                # 写入cds长度
                cds_length = int(list[j][4]) - int(list[j][3]) + 1
                cds_total_length += int(cds_length)
                listNew[-1] = cds_total_length
                cds_num2 += 1
            elif list[j][2] == 'five_prime_UTR' or list[j][2] == 'three_prime_UTR' or list[j][2] == 'UTR' \
                    or list[j][2] == '3_UTR' or list[j][2] == '5_UTR' or list[j][2] == 'UTR_5' or list[j][2] == 'UTR_3':
                # 将cds序列写入列表
                cds_location = list[j][0] + ':' + list[j][3] + ':' + list[j][4]
                listNew.insert(cds_num, cds_location)
                cds_num += 1
            else:
                print(list[j])
                listNew[2] = list[0][6]
        newGffWrite(newFile, listNew)

    # 翻译rna序列
    def rna_translate(seq):
        protein_table = {'UUU': 'F', 'CUU': 'L', 'AUU': 'I', 'GUU': 'V',
                         'UUC': 'F', 'CUC': 'L', 'AUC': 'I', 'GUC': 'V',
                         'UUA': 'L', 'CUA': 'L', 'AUA': 'I', 'GUA': 'V',
                         'UUG': 'L', 'CUG': 'L', 'AUG': 'M', 'GUG': 'V',
                         'UCU': 'S', 'CCU': 'P', 'ACU': 'T', 'GCU': 'A',
                         'UCC': 'S', 'CCC': 'P', 'ACC': 'T', 'GCC': 'A',
                         'UCA': 'S', 'CCA': 'P', 'ACA': 'T', 'GCA': 'A',
                         'UCG': 'S', 'CCG': 'P', 'ACG': 'T', 'GCG': 'A',
                         'UAU': 'Y', 'CAU': 'H', 'AAU': 'N', 'GAU': 'D',
                         'UAC': 'Y', 'CAC': 'H', 'AAC': 'N', 'GAC': 'D',
                         'UAA': 'Stop', 'CAA': 'Q', 'AAA': 'K', 'GAA': 'E',
                         'UAG': 'Stop', 'CAG': 'Q', 'AAG': 'K', 'GAG': 'E',
                         'UGU': 'C', 'CGU': 'R', 'AGU': 'S', 'GGU': 'G',
                         'UGC': 'C', 'CGC': 'R', 'AGC': 'S', 'GGC': 'G',
                         'UGA': 'Stop', 'CGA': 'R', 'AGA': 'R', 'GGA': 'G',
                         'UGG': 'W', 'CGG': 'R', 'AGG': 'R', 'GGG': 'G'}
        if len(seq) % 3 == 2:
            seq += 'N'
        elif len(seq) % 3 == 1:
            seq += 'NN'
        prot = ''
        for i in range(0, len(seq) + 1, 3):

            protein = seq[i:i + 3]

            if protein in protein_table:
                if protein_table[protein] == 'Stop':
                    prot += '*'
                else:
                    prot += protein_table[protein]
            else:
                prot += ''
        return prot

    # 将文件按照第一行和第三行排序（默认第一行为染色体位置）
    commond_2 = r"sort   -k1,1 -k4n,4   -t $'\t' " + gff_file
    sort_file = os.popen(commond_2, 'r')
    # 去掉无用行
    file_list = []
    for line in sort_file:
        lineList = line.split('\t')

        if line.startswith('#') or line.strip() == '':
            continue
        # maker注释的文件会将序列信息放在文件下部
        elif line.startswith('>Contig'):
            break
        # 判断是否为无用行
        elif ('miRNA' in lineList[8]) or ('ncRNA' in lineList[8]) or ('snoRNA' in lineList[8]) or (
                'snRNA' in lineList[8]) or ('tRNA' in lineList[8]) or ('rRNA' in lineList[8]) or (
                'other_RNA' in lineList[8]) or ('transposable_element_gene' in lineList[8]):
            continue
        elif lineList[2] == 'gene' or lineList[2] == 'mRNA' or lineList[2] == 'CDS':
            file_list.append(lineList)
        elif lineList[2] == 'transcript':
            lineList[2] = 'mRNA'
            file_list.append(lineList)
        else:
            continue

    # 判断文件类型
    s = 0
    r = 0
    file_type = ''
    for i in range(len(file_list)):
        if i == 1000:
            break
        elif file_list[i][2] == 'gene':
            s = 1
        elif file_list[i][2] == 'mRNA':
            r = 2
        elif s + r == 3:
            file_type = 'gene+mrna'
            break
    if s + r == 1:
        file_type = 'gene_only'
    elif s + r == 2:
        file_type = 'mrna_only'

    # 根据文件类型进行相应的操作
    gene_dic = OrderedDict()
    mrna_dic = OrderedDict()
    gene_mrna_dic = {}
    gene_list = []
    # 当文件类型为gene加mrna时
    if file_type == 'gene+mrna':
        for i in range(len(file_list)):
            list_i = file_list[i]
            # 将基因按照ID存进字典
            if list_i[2] == 'gene':
                gene_id = pick_id(list_i[8])[3:]
                gene_list.append(gene_id)
                gene_dic[gene_id] = list_i
            elif list_i[2] == 'mRNA':

                mrna_id = pick_id(list_i[8])[3:]
                if 'Parent' in list_i[8]:
                    mrna_parent = pick_parent(list_i[8])[7:]
                else:
                    mrna_parent = pick_id(list_i[8])[3:]
                # 判断
                gene_mrna_dic[mrna_id] = mrna_parent
                # 判断基因字典如果已经有这个键就加入否则就新建一个
                if mrna_id in mrna_dic.keys():
                    # 出现mrna就将其插入第一位
                    mrna_dic[mrna_id].insert(0, list_i)
                else:
                    mrna_dic[mrna_id] = []
                    mrna_dic[mrna_id].append(list_i)

            else:
                # print(file_list[i])
                # seq_id=pick_id(file_list[i][8])[3:]
                seq_parent = pick_parent(file_list[i][8])[7:]
                # 有一个gff文件中parent，写了两部分，以‘,’分隔，一部分是mrnaID 一部分是蛋白质ID
                if ',' in seq_parent:
                    seq_parent = seq_parent.split(',')[0]
                # 如果已经存有这个键就加入否则就新建一个
                if seq_parent in mrna_dic.keys():
                    mrna_dic[seq_parent].append(list_i)
                else:
                    mrna_dic[seq_parent] = []
                    mrna_dic[seq_parent].append(list_i)
        # 将gene加进来
        for mrna_key, value in mrna_dic.items():
            # 根据字典得到基因id
            # 判断这个mrna是否有对应的基因（有可能是转座子基因被跳过）
            if gene_mrna_dic[mrna_key] in gene_list:
                add_gene_id = gene_mrna_dic[mrna_key]
                mrna_dic[mrna_key].insert(0, gene_dic[add_gene_id])
            else:
                continue
            # 利用mrnaID找到这个字典，把对应的基因字典添加进来，并且插到第一位
            # 找不到对应的字典将其删除

            # 对嵌套的基因列表进行排序
            value.sort(key=lambda x: int(x[3]))
            # 调整第九列
            adjustCol9(value)
            new_gff(value)

    # 当文件中只有基因时
    elif file_type == 'gene_only':
        for i in range(len(file_list)):
            list_i = file_list[i]
            if list_i[2] == 'gene':

                gene_id = pick_id(list_i[8])[3:]
                if gene_id in gene_dic.keys():
                    gene_dic[gene_id].insert(0, list_i)
                else:
                    gene_dic[gene_id] = []
                    gene_dic[gene_id].append(list_i)
            else:
                seq_parent = pick_parent(list_i[8])[7:]
                if seq_parent in gene_dic.keys():

                    gene_dic[seq_parent].append(list_i)
                else:
                    gene_dic[seq_parent] = []
                    gene_dic[seq_parent].append(list_i)

        for keys, value in gene_dic.items():
            value.sort(key=lambda x: int(x[3]))
            # 调整第九列
            adjustCol9(value)
            new_gff(value)

    # 当文件中只有mrna时
    elif file_type == 'mrna_only':
        for i in range(len(file_list)):
            list_i = file_list[i]
            if list_i[2] == 'mRNA':
                gene_id = pick_id(list_i[8])[3:]
                if gene_id in gene_dic.keys():
                    gene_dic[gene_id].insert(0, list_i)
                else:
                    gene_dic[gene_id] = []
                    gene_dic[gene_id].append(list_i)
            else:
                seq_parent = pick_parent(list_i[8])[7:]
                if seq_parent in gene_dic.keys():
                    gene_dic[seq_parent].append(list_i)
                else:
                    gene_dic[seq_parent] = []
                    gene_dic[seq_parent].append(list_i)
        for keys, value in gene_dic.items():
            value.sort(key=lambda x: int(x[3]))
            # 调整第九列
            adjustCol9(value)
            new_gff(value)

    else:
        print ('错误文件')

    newFile.close()
    input_position_file = open(file_name + gff_title + '_new.gff', 'r')
    # 去除可变剪切
    gene_mrna_dict = OrderedDict()
    for line in input_position_file:
        line = line.strip()
        line_list = line.split('\t')
        # 留下最长的可变剪切行
        if line_list[0] in gene_mrna_dict.keys():
            # 如果遇到比已有字典里的长，就替换
            if line_list[-1] > gene_mrna_dict[line_list[0]][-1]:
                gene_mrna_dict[line_list[0]] = line_list
            else:
                continue
        else:
            gene_mrna_dict[line_list[0]] = []
            gene_mrna_dict[line_list[0]] = line_list
    input_position_file.close()
    # 将字典里处理过可变剪辑的行写入文件
    for value in gene_mrna_dict.values():
        for part in value:
            longest_mrna_file.write(part + '\t')
        longest_mrna_file.write('\n')

    longest_mrna_file.close()
    commond_4 = 'rm -rf ' + file_name + gff_title + '_new.gff'
    os.system(commond_4)
    # 对gff 文件进行排序
    sort_file_func(file_name + gff_title + '_longest.gff')
    # output_file = open(file_name + gff_title + '_longest.gff', 'w')
    # for key in sorted(total_line_dict.keys()):
    #     for line in total_line_dict[key]:
    #         output_file.write(line[1])
    # for line in spec_list:
    #     output_file.write(line[1])
    # output_file.close()
    longest_mrna_file = open(file_name + gff_title + '_longest.gff', 'r')
    # 将cds序列提取翻译写入文件
    fa_dict = SeqIO.to_dict(SeqIO.parse(fa_file, "fasta", IUPAC.unambiguous_dna))

    for line2 in longest_mrna_file:
        line2 = line2.strip()
        line2_list = line2.split('\t')
        # 去除空行
        if len(line2) == 0:
            continue
        cds_seq = Seq('', IUPAC.unambiguous_dna)
        for position in line2_list[4:-2]:
            position_list = position.split(':')
            chr_pos = position_list[0]
            # 位置坐标从1开始
            start = int(position_list[1]) - 1
            stop = int(position_list[2])
            cds_seq += fa_dict[chr_pos].seq[start:stop]
        # my_seq = Seq(str(cds_seq), IUPAC.unambiguous_dna)
        if line2_list[2] == '+':
            protein_seq = rna_translate(cds_seq.transcribe())
        else:
            protein_seq = rna_translate(cds_seq.reverse_complement().transcribe())
        output_pro_file.write('>' + line2_list[0] + '\n' + str(protein_seq) + '\n')
    output_pro_file.close()
    longest_mrna_file.close()
    get_chr_length(fa_file,output_chr_file)
def run_gene_fa(fa_file_list, gff_file_list, output_file,threads_num):
    pool = multiprocessing.Pool(threads_num)
    for fa_index in range(len(fa_file_list)):
        result = pool.apply_async(gene_fa,
                                  args=(gff_file_list[fa_index], fa_file_list[fa_index], fa_index, output_file))
    result.get()
    pool.close()
    pool.join()
if __name__ == '__main__':
    gff_file=sys.argv[1]
    fa_file=sys.argv[2]
    gene_fa(gff_file,fa_file, 100)