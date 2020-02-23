# -*- coding: utf-8 -*-
# 2020.01.06
# @zlk
# 将得到的共线性文件画出散点图
import pygal
import re
import glob
import sys
import os
from collections import OrderedDict
#
# syn_file=open(sys.argv[1],'r')
# pos1_file=open(sys.argv[2],'r')
# pos2_file=open(sys.argv[3],'r')

def sort_chr(file_name):
    file_name2=open(file_name,'r')
    total_list=[]
    for line in file_name2:
        line_list=line.strip().split('\t')
        total_list.append(line_list)
    total_length=0
    last_length=0
    chr_length_dict={}
    for i in total_list:
        total_length+=int(i[1])
        i[1]=total_length
        chr_length_dict[i[0]] = last_length
        last_length=i[1]
    file_name2.close()
    return total_list,total_length,chr_length_dict
def syn_pict(out_name,syn_file,pos1,pos2):
    x_list,x_length,x_length_dict=sort_chr(pos1)
    y_list,y_length,y_length_dict=sort_chr(pos2)
    point_list=[]
    syn_file_list=syn_file.strip().split('/')[-1].split('_')
    syn_file=open(syn_file,'r')
    for line in syn_file:
        line_list=line.strip().split('\t')
        line1_list=line_list[1].strip().split(':')
        line3_list=line_list[3].strip().split(':')
        point_list.append({'value':[x_length_dict[line1_list[0]]+int(line1_list[1]),y_length_dict[line3_list[0]]+int(line3_list[1])],'label':line_list[0][3:]+':'+line_list[2][3:]})
        # point_list.append([x_length_dict[line1_list[0]]+int(line1_list[1]),y_length_dict[line3_list[0]]+int(line3_list[1]),line_list[0]+':'+line_list[2]])
    xy_chart = pygal.XY(stroke=False,x_label_rotation=90,show_legend=False,dots_size=1,
                        x_title=syn_file_list[0],y_title=syn_file_list[1])
    xy_chart.config.style.title_font_size = 25
    xy_chart.config.style.label_font_size = 10
    xy_chart.config.show_x_guides =True
    xy_chart.config.width =1800
    xy_chart.config.height =1200
    xy_chart.title = 'SynOrths'
    xlable_list=[]
    for i in x_list:
        xlable_list.append({'label': i[0],'value': int(i[1])})
    xy_chart.x_labels = xlable_list
    ylable_list=[]
    for i in y_list:
        ylable_list.append({'label': i[0],'value': int(i[1])})
    xy_chart.y_labels = ylable_list
    xy_chart.add('', point_list)
    xy_chart.render_to_file(out_name+".svg")
def run_syn_pict(output_file):
    work_file_list = glob.glob(output_file + 'file?')
    file_num = len(work_file_list)
    for file_index in range(file_num):
        work_path = output_file + '/file' + str(file_index) + '/'
        files = os.listdir(work_path)
        for f in files:
            if f.endswith('_syn'):
                syn = work_path + f
                syn_list=f.strip().split('_')
                pos1=output_file+'/'+syn_list[0]+'/chr_length.txt'
                pos2=output_file+'/'+syn_list[1]+'/chr_length.txt'
                syn_pict(syn,syn,pos1,pos2)

if __name__ == '__main__':
    run_syn_pict('./')