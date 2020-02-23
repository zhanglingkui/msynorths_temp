# -*- coding: utf-8 -*-
# 2020.01.6
# @zlk
# 找到每个物种的共线性基因，和tandem基因，然后到blast结果中将这些同源去掉，剩下的就是distance duplication
import os
import multiprocessing

def pick_syn(syn_file):
    file1=open(syn_file,'r')
    syn_dict={}
    for line in file1:
        line_list=line.strip().split('\t')
        for i in line_list:
            num_gene_list=i.strip().split('=')
            if num_gene_list[0] in syn_dict.keys():
                syn_dict[num_gene_list[0]].append(num_gene_list[1])
            else:
                syn_dict[num_gene_list[0]]=[]
                syn_dict[num_gene_list[0]].append(num_gene_list[1])
    return syn_dict
def pick_tandem(tandem_file):
    tandem_list=[]
    file2=open(tandem_file,'r')
    for line in file2:
        line_list=line.strip().split('\t')
        for i in line_list:
            tandem_list.append(i)
    return tandem_list

def distance_def(output_file,loop_num):
    syn_dict=pick_syn(output_file+'/total_syn_split')
    for file_index in range(loop_num):
        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        syn_list=syn_dict[str(file_index)]
        for f in files:
            if f.endswith('tandem_array'):
                tandem_file = work_path + f
                tandem_list=pick_tandem(tandem_file)
        distance_dict={}
        for f in files:
            if f=='self_blast':
                for line in open(work_path+f,'r'):
                    line_list=line.strip().split('\t')
                    if line_list[0][3:] in syn_list:
                        if (line_list[1][3:] not in syn_list) and (line_list[1][3:] not in tandem_list):
                            if line_list[0] in distance_dict.keys():
                                distance_dict[line_list[0]].append(line_list[1])
                            else:
                                distance_dict[line_list[0]] = []
                                distance_dict[line_list[0]].append(line_list[1])

        output_distance=open(work_path+'distance_array.txt','w')
        for key,value in distance_dict.items():
            output_distance.write(key+'\t')
            for i in value:
                output_distance.write(i+'\t')
            output_distance.write('\n')
def run_distance(output_file,loop_num,threads):
    pool = multiprocessing.Pool(threads)
    distance_def(output_file,loop_num)
    result = pool.apply_async(distance_def,(output_file,loop_num))
    result.get()
    pool.close()
    pool.join()
if __name__ == '__main__':
    import sys

    tandem_list = pick_tandem(sys.argv[3])
    syn_dict = pick_syn(sys.argv[2])
    syn_list = syn_dict[str(1)]
    distance_dict={}
    for line in open(sys.argv[1], 'r'):

        line_list = line.strip().split('\t')
        if line_list[0][3:] in syn_list:
            if (line_list[1][3:] not in syn_list) and (line_list[1][3:] not in tandem_list):
                if line_list[0][3:] in distance_dict.keys():
                    distance_dict[line_list[0]].append(line_list[1])
                else:
                    distance_dict[line_list[0]] = []
                    distance_dict[line_list[0]].append(line_list[1])

    output_distance = open('distance_array.txt', 'w')
    for key, value in distance_dict.items():
        output_distance.write(key + '\t')
        for i in value:
            output_distance.write(i + '\t')
        output_distance.write('\n')