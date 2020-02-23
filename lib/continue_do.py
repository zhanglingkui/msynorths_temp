# -*- coding: utf-8 -*-
import os
import re
import glob
import msynorths2
import flanking_gene_total9
import synteny6
import tandem5
import add_species
# 在进行m模式时中断
def continue_do1(fa,gff,output_file,threads,tools):
    flag_list=[]
    fa_file_list=fa.strip().split(',')
    latest_file=''
    latest_time=0
    latest_file2 = ''
    latest_time2 = 0
    latest_file2_index=0
    file_file_list=[]
    rm_list=[]
    for file_index in range(len(fa_file_list)):
        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        for f in files:
            if f.endswith('db2'):
                rm_list.append(work_path+'db/')
            if 'self_blast' in f:
                if os.path.getmtime(work_path+f)>latest_time2:
                    latest_file2=work_path+f
                    latest_file2_index=file_index
                    latest_time2=os.path.getmtime(work_path+f)
                flag_list.append(1)
            elif f.endswith('_syn'):
                flag_list.append(2)
            elif re.search(r'file\d*_file\d',f):
                file_file_list.append(f)
                if os.path.getmtime(work_path+f)>latest_time:
                    latest_file=work_path+f
                    latest_time=os.path.getmtime(work_path+f)
                flag_list.append(3)
    #代表已经跑到最后一步，共线性分析
    if 2 in flag_list:
        work_file_list = glob.glob(output_file + 'file?')
        file_num=len(work_file_list)
        flanking_gene_total9.run_flanking(threads,output_file,file_num)
        synteny6.synteny(file_num)
    # 代表跑到两两blast了，删除最近产生的一个blast结果，从它开始继续跑
    elif 3 in flag_list:
        file_name_list=latest_file.strip().split('/')
        os.system('rm '+latest_file)
        file_file_list.remove(file_name_list[-1])
        # 将建的库都删掉
        for i in rm_list:
            os.system('rm -r' + i)
        # 将未进行两两blast的进行blast
        flanking_gene_total9.resp_blast(len(fa_file_list), output_file, tools, threads, file_file_list)
        flanking_gene_total9.run_flanking(threads,output_file,len(fa_file_list))
        synteny6.synteny(len(fa_file_list))
    # 代表跑到了self_blast
    elif 1 in flag_list:
        os.system('rm ' + latest_file2)
        os.system('rm -rf ' + output_file + '/file' + str(latest_file2_index) + '/db')
        do_list=[]
        for file_index in range(len(fa_file_list)):
            work_path = output_file + '/file' + str(file_index) + '/'
            files = os.listdir(work_path)
            for f in files:
                if 'self_blast' in f:
                    do_list.append(file_index)
        # 自身blast，用作去除tandem
        tandem5.self_blast(len(fa_file_list), output_file, threads, tools, do_list)
        # 挑选tandem，并去除在基因序列文件中去除tandem
        tandem5.run_pick_tandem(output_file, threads, len(fa_file_list))
        # 两两物种进行blast
        flanking_gene_total9.resp_blast(len(fa_file_list), output_file, tools, threads, [])
        # 对两两的blast结果进行共线性分析
        flanking_gene_total9.run_flanking(threads, output_file, len(fa_file_list))
        # 进行多基因组共线性分析
        synteny6.synteny(len(fa_file_list))
    elif len(flag_list) == 0:
        for file_index in range(len(fa_file_list)):
            work_path = output_file + '/file' + str(file_index) + '/'
            os.system('rm -r ' + work_path)
        msynorths2.msynorths(fa, gff, threads, output_file,tools)
# 在a模式时中断
def continue_do2(fa,gff,output_file,threads,file_num,tools):
    flag_list = []
    fa_file_list = fa.strip().split(',')
    latest_file = ''
    latest_time = 0
    latest_file2 = ''
    latest_time2 = 0
    latest_file2_index = 0
    file_file_list = []
    rm_list=[]
    for file_index in range(file_num-1,file_num-len(fa_file_list)-1,-1):

        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        for f in files:
            if f.endswith('db2'):
                rm_list.append(work_path+'db/')
            if 'self_blast' in f:
                if os.path.getmtime(work_path + f) > latest_time2:
                    latest_file2 = work_path + f
                    latest_file2_index = file_index
                    latest_time2 = os.path.getmtime(work_path + f)
                flag_list.append(1)
            elif f.endswith('_syn'):
                flag_list.append(2)
            elif re.search(r'file\d*_file\d', f):
                file_file_list.append(f)
                if os.path.getmtime(work_path + f) > latest_time:
                    latest_file = work_path + f
                    latest_time = os.path.getmtime(work_path + f)
                flag_list.append(3)
        if 2 in flag_list:
            flanking_gene_total9.run_flanking(threads, output_file, file_num)
            synteny6.synteny(file_num)
        elif 3 in flag_list:
            file_file_list = []
            for file_index in range(file_num - 1, file_num - len(fa_file_list) - 1, -1):
                work_path = output_file + '/file' + str(file_index) + '/'
                files = os.listdir(work_path)
                for f in files:
                    if re.search(r'file\d*_file\d', f):
                        file_file_list.append(f)
            file_name_list = latest_file.strip().split('/')
            os.system('rm ' + latest_file)
            file_file_list.remove(file_name_list[-1])

            # 将建的库都删掉
            for i in rm_list:
                os.system('rm -r' + i)
            flanking_gene_total9.resp_blast(file_num, output_file, tools, threads, file_file_list)
            flanking_gene_total9.run_flanking(threads, output_file, file_num)
            synteny6.synteny(file_num)

        elif 1 in flag_list:
            os.system('rm ' + latest_file2)
            os.system('rm -rf ' + output_file + '/file' + str(latest_file2_index) + '/db')
            do_list = []
            for file_index in range(file_num):
                work_path = output_file + '/file' + str(file_index) + '/'
                files = os.listdir(work_path)
                for f in files:
                    if 'self_blast' in f:
                        do_list.append(file_index)
            # 自身blast，用作去除tandem
            tandem5.self_blast(file_num, output_file, threads, tools, do_list)
            # 挑选tandem，并去除在基因序列文件中去除tandem
            tandem5.run_pick_tandem(output_file, threads, file_num)
            # 两两物种进行blast
            flanking_gene_total9.resp_blast(file_num, output_file, tools, threads, [])
            # 对两两的blast结果进行共线性分析
            flanking_gene_total9.run_flanking(threads, output_file, file_num)
            # 进行多基因组共线性分析
            synteny6.synteny(file_num)
        elif len(flag_list) == 0:
            for file_index in range(file_num - 1, file_num - len(fa_file_list) - 1, -1):
                work_path = output_file + '/file' + str(file_index) + '/'
                os.system('rm -r ' + work_path)
            add_species.add_species(fa, gff, threads, output_file, file_num-len(fa),tools)