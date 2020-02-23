# -*- coding: utf-8 -*-
# 2019.12.16
# @zlk
# 根据两个基因组的blastp结果，计算其flanking基因比值，得出结果文件，后续将对其过滤筛选(前几版本程序均是将同源基因写入字典进行配对，后发现没有必要)。
# 之前计算侧翼基因用的是所有的blastp的同源基因对，现在使用best hit的同源基因对, 现在要对blast结果进行进一步过滤要与tandem过滤条件相同
# 原有版本需要将潜在的共线性基因写入文件，再读取判断。现在直接将其写入字典
# 原本需要进行的两次共线性flanking判断是先后进行，现在决定对其进行多进程处理，加快运行速度。
# 11.22 原来没有考虑将不同染色体的进行单独分析，可能会导致染色体边缘会出现问题，还有就是一些乱序的scaffold可能会出现问题
# 原有版本需要跑两个小时才有结果太慢了，需要改进以染色体为单位，进行存储顺序坐标
# 原有版本在挑选blast前五个结果时，只考虑了一种情况，现在正反都要做。
from collections import OrderedDict
import os
import glob
import multiprocessing
from multiprocessing import Process, Manager, Pool
import sys
import math
import synteny6


def pick_file(path, end_str):
    files1 = os.listdir(path)
    for f1 in files1:
        if f1.endswith(end_str):
            file_path = path + f1
    return file_path
def flanking_gene(gff1_title, gff2_title, blast_file, outfile):
    input_blast = open(blast_file, 'r')
    input_identity = 0
    input_coverage = 0
    order_num = 0.5
    flanking_control_num = 4
    input_evalue = math.exp(-10)



    def filter_blast(blast_file1):
        front_dict = {}
        reverse_dict = {}
        for line2 in blast_file1:
            line_list2 = line2.strip().split('\t')
            if line_list2[0] in front_dict.keys():
                if len(front_dict[line_list2[0]]) < 7:
                    front_dict[line_list2[0]].append([line_list2[-2], line2])
            else:
                front_dict[line_list2[0]] = []
                front_dict[line_list2[0]].append([line_list2[-2], line2])
            if line_list2[1] in reverse_dict.keys():
                if len(reverse_dict[line_list2[1]]) < 7:
                    reverse_dict[line_list2[1]].append([line_list2[-2], line2])
                else:
                    reverse_dict[line_list2[1]] = sorted(reverse_dict[line_list2[1]], key=lambda p: float(p[0]),
                                                        reverse=False)
                    if line_list2[-2] < reverse_dict[line_list2[1]][-1][0]:
                        reverse_dict[line_list2[1]][-1] = [line_list2[-2], line2]
            else:
                reverse_dict[line_list2[1]] = []
                reverse_dict[line_list2[1]].append([line_list2[-2], line2])
        reverse_list = []
        for key1, value1 in reverse_dict.items():
            for line3 in value1:
                reverse_list.append(line3[1])
        front_list=[]
        for key2, value2 in front_dict.items():
            for line4 in value2:
                front_list.append(line4[1])
        if len(reverse_list)>len(front_list):
            output_list=reverse_list
        else:
            output_list=front_list
        return output_list

    def write_list(w_list):
        output_str = ''
        for i in w_list:
            output_str += (str(i) + ',')
        return output_str

    # 将checklist里出现一次以上的元素输出一个列表
    def do_check_list(do_list):
        do_dict = {}
        return_list = []
        for i in do_list:
            if i in do_dict.keys():
                return_list.append(i)
            else:
                do_dict[i] = 1
        return return_list

    def compare_20to100(query_start, query_stop, query_chr, ref_start, ref_stop, ref_chr, homgene_dict):
        flanking_num = 0
        flanking_index_list = []
        for query20_gene in query_gene_list_dict[query_chr][query_start:query_stop]:
            if query20_gene in homgene_dict.keys():
                for ref100_gene in ref_gene_list_dict[ref_chr][ref_start:ref_stop]:
                    if homgene_dict[query20_gene][0] == ref100_gene:
                        flanking_num += 1
                        flanking_index_list.append(ref_gene_index_dict[ref_chr][ref100_gene][1])
                        break
            else:
                continue
        return flanking_num, flanking_index_list

    def compare2_20to100(query_start, query_stop, query_chr,ref_start, ref_stop, ref_chr,homgene_dict):
        flanking_num = 0
        flanking_index_list = []
        for ref100_gene in ref_gene_list_dict[ref_chr][ref_start:ref_stop]:
            if ref100_gene in homgene_dict.keys():
                for query20_gene in query_gene_list_dict[query_chr][query_start:query_stop]:
                    if homgene_dict[ref100_gene][0] == query20_gene:
                        flanking_num += 1
                        flanking_index_list.append(query_gene_index_dict[query_chr][query20_gene][1])
                        break
            else:
                continue
        return flanking_num, flanking_index_list

    def big_syn(input_file, syn_dict, query_compare_num, ref_compare_num, compare_def, output_gene_dict,
                output_syn_list):
        input_file1 = open(input_file, 'r')
        input_file1_list=filter_blast(input_file1)
        for line4 in input_file1_list:

            line_list=line4.strip().split()
            if float(line_list[2]) < input_identity:
                continue
            elif float(line_list[10]) > input_evalue:
                continue
            query_gene = line_list[0]
            query_gene_chr = query_gene_position_dict[query_gene].split(':')[0]
            query_gene_index = int(query_gene_index_dict[query_gene_chr][query_gene][1])
            total_query_gene_num = len(query_gene_index_dict[query_gene_chr])
            ref_gene = line_list[1]
            ref_gene_chr = ref_gene_position_dict[ref_gene].split(':')[0]
            ref_gene_index = int(ref_gene_index_dict[ref_gene_chr][ref_gene][1])
            total_ref_gene_num = len(ref_gene_index_dict[ref_gene_chr])
            # 设置coverage
            if int(line_list[3].strip()) / int(query_gene_index_dict[query_gene_chr][line_list[0]][0]) < input_coverage:
                continue
            # 判断query gene前面是否有20个基因
            if query_gene_index > query_compare_num:
                pre_gene_num = query_compare_num
                query_gene_pre_index = query_gene_index - query_compare_num
            else:
                pre_gene_num = query_gene_index
                query_gene_pre_index = 1
            # 判断下ref 基因是否背后有100个基因
            if total_ref_gene_num - ref_gene_index < ref_compare_num:
                behind_gene_num2 = total_ref_gene_num - ref_gene_index
                ref_gene_behind_index = total_ref_gene_num
            else:
                behind_gene_num2 = ref_compare_num
                ref_gene_behind_index = ref_gene_index + ref_compare_num
            # 判断下ref gene前面是否有100个基因
            if ref_gene_index > ref_compare_num:
                pre_gene_num2 = ref_compare_num
                ref_gene_pre_index = ref_gene_index - ref_compare_num
            else:
                pre_gene_num2 = ref_gene_index
                ref_gene_pre_index = 1
            # 判断下query gene 是否背后有20个基因
            if total_query_gene_num - query_gene_index < query_compare_num:
                behind_gene_num = total_query_gene_num - query_gene_index
                query_gene_behind_index = total_query_gene_num
            else:
                behind_gene_num = query_compare_num
                query_gene_behind_index = query_gene_index + query_compare_num

            # ========================
            flanking_pre_pre_num, pre_pre_index_list = compare_def(query_gene_pre_index - 1, query_gene_index - 1,query_gene_chr,
                                                                   ref_gene_pre_index - 1,
                                                                   ref_gene_index - 1,ref_gene_chr, syn_dict)

            # ========================
            flanking_behind_behind_num, behind_behind_index_list = compare_def(query_gene_index,
                                                                               query_gene_behind_index,query_gene_chr,
                                                                               ref_gene_index,
                                                                               ref_gene_behind_index,ref_gene_chr, syn_dict)

            # ========================
            flanking_pre_behind_num, pre_behind_index_list = compare_def(query_gene_pre_index - 1, query_gene_index - 1,query_gene_chr,
                                                                         ref_gene_index,
                                                                         ref_gene_behind_index,ref_gene_chr, syn_dict)

            # ========================
            flanking_behind_pre_num, behind_pre_index_list = compare_def(query_gene_index, query_gene_behind_index,query_gene_chr,
                                                                         ref_gene_pre_index - 1,
                                                                         ref_gene_index - 1,ref_gene_chr, syn_dict)

            if flanking_pre_pre_num >= flanking_control_num or flanking_behind_behind_num >= flanking_control_num:
                pre_ratio = str(flanking_pre_pre_num) + '/' + str(pre_gene_num)
                pre_ratio2 = str(flanking_pre_pre_num) + '/' + str(pre_gene_num2)
                behind_ratio = str(flanking_behind_behind_num) + '/' + str(behind_gene_num)
                behind_ratio2 = str(flanking_behind_behind_num) + '/' + str(behind_gene_num2)
                index_list = list(set(pre_pre_index_list + behind_behind_index_list))
                # 需要判断下输出list的排列顺序，以ref为第一前两列
                if output_gene_dict == syn_gene_dict:
                    gene_id = line_list[0]
                    syn_list = [line_list[0], query_gene_position_dict[line_list[0]], line_list[1],
                                ref_gene_position_dict[line_list[1]], line_list[10], pre_ratio, behind_ratio, '+',
                                write_list(index_list)]
                else:
                    gene_id = line_list[1]
                    syn_list = [line_list[0], query_gene_position_dict[line_list[0]], line_list[1],
                                ref_gene_position_dict[line_list[1]], line_list[10], pre_ratio2, behind_ratio2, '+',
                                write_list(index_list)]

                if gene_id in output_gene_dict.keys():
                    output_gene_dict[gene_id].append(syn_list)
                else:
                    output_gene_dict[gene_id] = []
                    output_gene_dict[gene_id].append(syn_list)
            if flanking_pre_behind_num >= flanking_control_num or flanking_behind_pre_num >= flanking_control_num:
                pre_ratio = str(flanking_pre_behind_num) + '/' + str(pre_gene_num)
                pre_ratio2 = str(flanking_pre_behind_num) + '/' + str(pre_gene_num2)
                behind_ratio = str(flanking_behind_pre_num) + '/' + str(behind_gene_num)
                behind_ratio2 = str(flanking_behind_pre_num) + '/' + str(behind_gene_num2)
                index_list = list(set(pre_behind_index_list + behind_pre_index_list))
                # 需要判断下输出list的排列顺序，以ref为第一前两列

                if output_gene_dict == syn_gene_dict:
                    gene_id = line_list[0]
                    syn_list = [line_list[0], query_gene_position_dict[line_list[0]], line_list[1],
                                ref_gene_position_dict[line_list[1]], line_list[10], pre_ratio, behind_ratio, '-',
                                write_list(index_list)]
                else:
                    gene_id = line_list[1]
                    syn_list = [line_list[0], query_gene_position_dict[line_list[0]], line_list[1],
                                ref_gene_position_dict[line_list[1]], line_list[10], pre_ratio2, behind_ratio2, '-',
                                write_list(index_list)]

                if gene_id in output_gene_dict.keys():
                    output_gene_dict[gene_id].append(syn_list)
                else:
                    output_gene_dict[gene_id] = []
                    output_gene_dict[gene_id].append(syn_list)
        input_file1.close()
        evalue = 1
        newlist = []
        for key, value in output_gene_dict.items():
            big_num = 0
            if len(value) == 1:
                output_syn_list.append(value[0])
                continue
            else:
                for list1 in value:
                    pre_num = int(list1[-4].split('/')[0])
                    behind_num = int(list1[-3].split('/')[0])
                    total_num = pre_num + behind_num

                    if total_num > big_num:
                        big_num = total_num
                        newlist = list1
                        evalue = float(list1[-5])
                    elif total_num < big_num:
                        continue
                    else:
                        if float(list1[-5]) <= evalue:
                            big_num = total_num
                            newlist = list1
                            evalue = float(list1[-5])
                        else:
                            continue
                output_syn_list.append(newlist)

    # 此函数用于获取，基因的位置坐标字典
    def get_position(input_gff, gene_index_dict, gene_list_dict, gene_position_dict):
        input_file = open(input_gff, 'r')
        index1 = 0
        chr_name = 'zzzz'
        for line in input_file:
            line_list = line.strip().split('\t')
            line3_list = line_list[3].strip().split(':')
            # 将物理位置存起来
            gene_position_dict[line_list[0]] = line_list[3]
            # 用于判断是否还在同一条染色体
            if line3_list[0] != chr_name:
                chr_name = line3_list[0]
                index1 = 0
            if line3_list[0] in gene_list_dict.keys():
                gene_list_dict[line3_list[0]].append(line_list[0])
            else:
                gene_list_dict[line3_list[0]] = []
                gene_list_dict[line3_list[0]].append(line_list[0])
            # 返回基因长度和index
            if line3_list[0] in gene_index_dict.keys():
                gene_index_dict[line3_list[0]][line_list[0]] = [line_list[-1], index1]
            else:
                gene_index_dict[line3_list[0]] = {}
                gene_index_dict[line3_list[0]][line_list[0]] = [line_list[-1], index1]
            index1 += 1
        input_file.close()

    # 拿到ref基因组的顺序字典，基因列表字典，物理位置字典
    ref_gene_index_dict = {}
    ref_gene_list_dict = {}
    ref_gene_position_dict = {}
    get_position(gff1_title, ref_gene_index_dict, ref_gene_list_dict, ref_gene_position_dict)
    # 拿到query基因组的顺序字典，基因列表字典，物理位置字典
    query_gene_index_dict = {}
    query_gene_list_dict = {}
    query_gene_position_dict = {}
    get_position(gff2_title, query_gene_index_dict, query_gene_list_dict, query_gene_position_dict)
 # 将同源基因的best hit记录下来
    homgene_dict1 = {}
    homgene_dict2 = {}
    for line in input_blast:
        line_list = line.split('\t')
        if line_list[0] in homgene_dict1.keys():
            if float(line_list[2]) > homgene_dict1[line_list[0]][1]:
                homgene_dict1[line_list[0]] = [line_list[1], float(line_list[2])]
        else:
            homgene_dict1[line_list[0]] = [line_list[1], float(line_list[2])]
        if line_list[1] in homgene_dict2.keys():
            if float(line_list[2]) > homgene_dict2[line_list[1]][1]:
                homgene_dict2[line_list[1]] = [line_list[0], float(line_list[2])]
        else:
            homgene_dict2[line_list[1]] = [line_list[0], float(line_list[2])]
    input_blast.close()
    # 遍历这个blast，判断每个同源基因对的共线性关系，过了门槛的将其加入到，共线性基因对的字典里。
    # 这两个字典分别以query 和ref 为键，判断下哪个字典的共线性基因多，就以哪个为最后的结果

    syn_gene_dict = OrderedDict()
    syn_gene_dict2 = OrderedDict()
    syn_gene_list = []
    syn_gene_list2 = []
    if len(query_gene_position_dict)>len(ref_gene_position_dict):
        query_num=20
        ref_num=100
        big_syn(blast_file, homgene_dict1, query_num, ref_num, compare_20to100, syn_gene_dict, syn_gene_list)
        big_syn(blast_file, homgene_dict2, query_num, ref_num, compare2_20to100, syn_gene_dict2, syn_gene_list2)

    else:
        query_num = 100
        ref_num = 20
        big_syn(blast_file, homgene_dict1, query_num, ref_num, compare_20to100, syn_gene_dict, syn_gene_list)
        big_syn(blast_file, homgene_dict2, query_num, ref_num, compare2_20to100, syn_gene_dict2, syn_gene_list2)
    if len(syn_gene_list) > len(syn_gene_list2):
        output_syn_file = open(outfile, 'w')
        more_gene_list = syn_gene_list
    else:
        output_syn_file = open(outfile, 'w')
        more_gene_list = []
        # 将列表按照ref基因排序
        chr_dict=OrderedDict()
        for syn_line in syn_gene_list2:
            line3_list=syn_line[3].split(':')
            if line3_list[0] in chr_dict.keys():
                chr_dict[line3_list[0]].append([line3_list[1], syn_line])
            else:
                chr_dict[line3_list[0]] = []
                chr_dict[line3_list[0]].append([line3_list[1], syn_line])
        for key, value in chr_dict.items():
            value.sort(key=lambda x: int(x[0]))
            for line in value:
                more_gene_list.append(line[1])


    # 将最后一行写入一个总的列表
    total_line8_list = []
    for line in more_gene_list:
        line8_list = line[-1].split(',')
        total_line8_list.append(line8_list)
    # 最后根据最后一行check一遍，将错误的index写入列表
    index = 0
    total_index = len(total_line8_list) - 1
    wrong_syn_list = []
    for line in total_line8_list:
        right_num = 0
        if index == 0:
            check_num = [1, 2, 3, 4]
        elif index == 1:
            check_num = [0, 2, 3, 4]
        elif total_index - index == 0:
            check_num = [index - 4, index - 3, index - 2, index - 1]
        elif total_index - index == 1:
            check_num = [index - 3, index - 2, index - 1, total_index]
        else:
            check_num = [index - 2, index - 1, index + 1, index + 2]
        check_list_tmp = []
        for i in check_num:
            check_list_tmp += total_line8_list[i]
        check_list = do_check_list(check_list_tmp)
        for i in line:
            if i in check_list:
                right_num += 1
            else:
                continue
        # 如果比对上的flanking基因数量占自己flanking基因数量的比值小于0.8 就认为它有问题
        if right_num / len(line) < order_num:
            wrong_syn_list.append(index)
        index += 1
    # 遍历这个需要删除的列表，切记列表的下标会随着删除而改变
    del_num = 0
    for i in wrong_syn_list:
        del more_gene_list[i - del_num]
        del_num += 1

    for syn_list in more_gene_list:
        for gene_list in syn_list[:-1]:
            output_syn_file.write(gene_list + '\t')
        output_syn_file.write('\n')
def resp_blast(loop_num,output_file,tools, threads_num,file_file_list):
    for file_index in range(loop_num):

        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        for f in files:
            if f.endswith('_prot_notandem'):
                prot_file = work_path + f
        if tools=='diamond':
            commond11 = 'mkdir ' + work_path + 'db2/'
            os.system(commond11)
            commond1='diamond makedb --in '+prot_file+' -d '+ work_path + 'db2/ref'
            os.system(commond1)
        else:
            commond1 = 'makeblastdb -in ' + prot_file + ' -dbtype prot -parse_seqids -out ' + work_path + 'db2/ref '
            os.system(commond1)
    # 第一遍遍历，从第一个文件到倒数第二个文件，每一次都是得到一个query基因组的文件夹
    for file_index in range(loop_num - 1, 0, -1):
        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        for f in files:
            if f.endswith('_notandem'):
                query_fa = work_path + f

        # 第二遍遍历，从当前query文件到最后一个文件，每一次得到的都是库的文件夹
        for db_file_index in range(file_index - 1, -1, -1):
            if ('file' + str(file_index) + '_file' + str(db_file_index)) in file_file_list:
                continue
            db_file = output_file + '/file' + str(db_file_index) + '/db2/ref'
            if tools=='diamond':
                blast_commond='diamond blastp -d  '+ db_file + ' -q '+query_fa+' --more-sensitive -e 1e-5  --quiet -p '+str(threads_num)+' -o '+ work_path + 'file' + str(
                file_index) + '_file' + str(db_file_index)
            else:
                blast_commond = 'blastp -db ' + db_file + ' -query ' + query_fa + ' -num_threads '+str(threads_num)+' -outfmt 6 -evalue 1e-10 -out ' + work_path + 'file' + str(
                file_index) + '_file' + str(db_file_index)
            os.system(blast_commond)
def run_flanking(threads,output_file,file_num):
    pool = multiprocessing.Pool(threads)
    for file_index in range(file_num - 1, 0, -1):
        work_path = output_file + '/file' + str(file_index) + '/'
        # 同一个文件夹下可能得到多个文件，会形成一个列表遍历这个列表进行操作
        work_file_list = glob.glob(work_path + 'file*_file*')
        for work_file_path in work_file_list:
            if work_file_path.endswith('_syn'):
                continue
            path_list = work_file_path.split('/')
            file_list = path_list[-1].split('_')
            file1 = output_file + '/' + file_list[0] + '/'
            file2 = output_file + '/' + file_list[1] + '/'
            result = pool.apply_async(flanking_gene,
                                      (pick_file(file2, '_longest.gff'), pick_file(file1, '_longest.gff'),
                                       work_file_path, work_file_path + '_syn'))
    result.get()
    pool.close()
    pool.join()

# if __name__ == '__main__':
    # gff1_title = sys.argv[2]
    # gff2_title = sys.argv[1]
    # blast_file = sys.argv[3]
    # outfile = sys.argv[4]
    #
    # flanking_gene(gff1_title, gff2_title, blast_file, outfile)
#     def pick_file(path, end_str):
#         files1 = os.listdir(path)
#         for f1 in files1:
#             if f1.endswith(end_str):
#                 file_path = path + f1
#         return file_path
#     now_path = os.getcwd()
#     pool = multiprocessing.Pool(50)
#     for file_index in range(3 - 1, 0, -1):
#         work_path = now_path + '/file' + str(file_index) + '/'
#         # 同一个文件夹下可能得到多个文件，会形成一个列表遍历这个列表进行操作
#         work_file_list = glob.glob(work_path + 'file*_file*')
#         for work_file_path in work_file_list:
#             if work_file_path.endswith('_syn'):
#                 continue
#             path_list = work_file_path.split('/')
#             file_list = path_list[-1].split('_')
#             file1 = now_path + '/' + file_list[0] + '/'
#             file2 = now_path + '/' + file_list[1] + '/'
#             result = pool.apply_async(flanking_gene,(pick_file(file2, '_longest.gff'), pick_file(file1, '_longest.gff'),
#                                        work_file_path, work_file_path + '_syn'))
#     pool.close()
#     pool.join()
# synteny6.synteny(3)
