# -*- coding: utf-8 -*-
# 2019.12.05
# 对基因ID进行修改，作为其后续排列的证据
# 使用社区发现算法，对每行的基因进行模块度计算。
# 对于可能出现的注释错误，就是一个基因可能注释成两个基因这种，在最后的total_syn 这个列表中
# 将这种基因给与一定标注，以便模块度计算可以跳过
# 将对每行共线性基因的分类算法，换成cpm派系过滤算法
# 在社区划分时，会将单独一条线或者两条线的基因不显示在分的社区中，现在将他拉回到现有的社区中
# 将划分好的社区进行分行，同时把基因分开，还有关系文件
import os
import multiprocessing
import glob
from collections import OrderedDict
import gc
import cpm_algo3
import re
def get_set(set_list):
    index = 1
    line_num = 5000
    gene_num_dict = {}
    for line in set_list:
        line_list = line
        if index > line_num:
            break
        for part in line_list:
            if part == 'exemption' or part == '':
                continue
            part_list = part
            if len(part_list) in gene_num_dict.keys():
                gene_num_dict[len(part_list)] += 1
            else:
                gene_num_dict[len(part_list)] = 1
            index += 1
    # average_num=total_num/index
    # average_num+=average_num/3\
    value_100 = 0
    for key, value in gene_num_dict.items():
        value_100 += value
    for value in sorted(gene_num_dict.items(), key=lambda d: d[0], reverse=True):
        if value[1] / value_100 > 0.1:
            re_value=value[0]
            break
    return re_value
def synteny(file_num):
    file_name = OrderedDict()
    names = globals()
    now_path = os.getcwd()

    # 将所有共线性文件的路径写入字典，便于下一步进行遍历
    for file_index in range(file_num-1, -1, -1):
        work_path = now_path + '/file' + str(file_index) + '/'
        work_file_list = glob.glob(work_path + 'file*_file*')
        gff_file_list = glob.glob(work_path + '*_longest.gff')
        names['gene_index_dict' + str(file_index)] = {}

        # 遍历这些共线性文件，以ref为键，添加到字典里
        for syn_file in work_file_list:
            if syn_file.endswith('_syn'):
                s_list = syn_file.split('/')
                sn_list = s_list[-1].split('_')
                if sn_list[0] in file_name.keys():
                    file_name[sn_list[0]].append(sn_list[1])
                else:
                    file_name[sn_list[0]] = []
                    file_name[sn_list[0]].append(sn_list[1])
        # 将每个基因组的基因index写入字典
        for gff_file in gff_file_list:
            gene_index_dict=names.get('gene_index_dict'+str(file_index))
            input_file=open(gff_file,'r')
            index=0
            for line in input_file:
                line_list=line.strip().split('\t')
                gene_index_dict[line_list[0][3:]]=index
                index+=1
            input_file.close()
    # 对字典里的value进行排序，根据其query id。
    for key, value in file_name.items():
        file_name[key] = sorted(value, key=lambda i: i, reverse=True)

    total_syn_list = []
    modularity_list = []
    # 建一个字典将所有双向的共线性基因对都存起来，
    total_syn_dict = {}

    for syn_file in file_name.keys():

        for ref_file in file_name[syn_file]:

            # 挑出文件的序号
            syn_index = re.sub("\D", "", syn_file)
            ref_index = re.sub("\D", "", ref_file)
            syn_file_name = syn_file + '_' + ref_file + '_syn'
            syn_file_path = now_path + '/' + syn_file + '/' + syn_file_name
            open_file = open(syn_file_path, 'r')
            # 将ID加入字典时对其进行了修改，增加了其所属的文件index
            for line in open_file:

                line_list = line.split('\t')
                new_query_id = syn_index + '=' + line_list[0].split('=')[1]
                new_ref_id = ref_index + '=' + line_list[2].split('=')[1]
                # 按照正向将共线性基因对加入字典
                if new_query_id in total_syn_dict.keys():
                    total_syn_dict[new_query_id].append(new_ref_id)
                else:
                    total_syn_dict[new_query_id] = []
                    total_syn_dict[new_query_id].append(new_ref_id)
                # 按照反向将共线性基因对加入字典
                if new_ref_id in total_syn_dict.keys():
                    total_syn_dict[new_ref_id].append(new_query_id)
                else:
                    total_syn_dict[new_ref_id] = []
                    total_syn_dict[new_ref_id].append(new_query_id)
            open_file.close()

    # 遍历字典，将所有的共线性基因拉进来
    total_syn_dict_tem = total_syn_dict.copy()
    for key, value in total_syn_dict.items():
        if key in total_syn_dict_tem.keys():
            add_num = 1
            total_add_list = []
            # 将第一次遍历的共线性关系对加入模块列表
            for value_index in value:
                add_modu_list = [0, add_num]
                add_num += 1
                total_add_list.append(add_modu_list)
            add_list = [key] + value
            del total_syn_dict_tem[key]
            for id in add_list:
                if id in total_syn_dict_tem.keys():
                    for new_id in total_syn_dict[id]:
                        # 看下这个基因是否已经被拉进来过了
                        if new_id in total_syn_dict_tem.keys():
                            pull_id_index = add_list.index(id)
                            if new_id in add_list:
                                # 获取ID位置
                                new_id_index = add_list.index(new_id)
                            else:
                                add_list.append(new_id)
                                # 查看当前新添加的ID为第几个
                                new_id_index = len(add_list) - 1
                            add_modu_list = [pull_id_index, new_id_index]
                            if (add_modu_list in total_add_list) or ([new_id_index, pull_id_index] in total_add_list):
                                continue
                            else:
                                total_add_list.append(add_modu_list)
                    del total_syn_dict_tem[id]

            modularity_list.append(total_add_list)
            total_syn_list.append(add_list)

    # 现在需要一个标准集，
    set_num=get_set(modularity_list)

    index=-1
    total_write_list=[]
    for line in modularity_list:
        index+=1
        part_list=cpm_algo3.make_graph(line,3,set_num)
        write_list=list(part_list)
        total_write_list.append(write_list)
    output_file_gene=open('total_syn_split','w')
    output_file_modu=open('modu_file_split','w')

    line_index=0
    for write_list in total_write_list:
        # 拿到基因id和序号的对应字典
        id_index_dict={}
        id_index=0
        for gene_id in total_syn_list[line_index]:
            id_index_dict[id_index]=gene_id
            id_index+=1
        for part in write_list:
            new_id_index=0
            new_id_index_dict={}
            for i in part:
                output_file_gene.write(total_syn_list[line_index][i]+'\t')
                new_id_index_dict[total_syn_list[line_index][i]]=new_id_index
                new_id_index+=1
            output_file_gene.write('\n')
            for modu_list in modularity_list[line_index]:
                # print (modu_list)
                if (modu_list[0] in part) and (modu_list[1] in part):
                    output_file_modu.write(str(new_id_index_dict[id_index_dict[modu_list[0]]])+' '+str(new_id_index_dict[id_index_dict[modu_list[1]]])+'\t')
            output_file_modu.write('\n')
        line_index +=1




if __name__ == "__main__":
    synteny(9)
