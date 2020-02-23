# -*- coding: utf-8 -*-
import argparse
import sys
import os
import glob
import re
current_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(current_dir+"/lib")
import msynorths2
import flanking_gene_total9
import add_species
import continue_do
import syn_pict
import distance_dup

parser = argparse.ArgumentParser(description='msynorths')
parser.add_argument('f', type=str, help='input config file,The first line enters the gff file location,the second line enters the fasta file location, use , to split every file')
parser.add_argument('-o', type=str, default=os.getcwd(), help='output file')
parser.add_argument('-t', type=int, default=10,help='input num of threads')
parser.add_argument("-m", "--msynorths", action='store_true', help="msynoths")
parser.add_argument('-c', "--continuedo", action='store_true',help="If the program breaks,you can use this parameter to continue")
parser.add_argument('-a', "--add", action='store_true',help="add specie")
parser.add_argument('-s',  type=str, help='blast or diamond')
args = parser.parse_args()
fa_gff_file = args.f
compare_tools=args.s
threads = args.t
output_file = args.o+'/'
now_path = os.getcwd()
work_file_list = glob.glob(output_file + 'file?')
file_num = len(work_file_list)
fa=''
gff=''
for line in open(fa_gff_file,'r'):
    line_list=line.strip().split(',')
    if len(line_list)==2:
        fa_abs=os.path.abspath(line_list[0])
        fa+=(fa_abs+',')
        gff_abs=os.path.abspath(line_list[1])
        gff += (gff_abs + ',')
fa=fa[:-1]
gff=gff[:-1]
print (fa)
if args.continuedo == True:
    if args.msynorths == True:
        continue_do.continue_do1(fa,gff,output_file,threads,compare_tools)
    elif args.add == True:
        continue_do.continue_do2(fa, gff, output_file, threads,file_num,compare_tools)
elif args.msynorths == True:
    msynorths2.msynorths(fa,gff,threads,output_file,compare_tools)
elif args.add == True:
    add_species.add_species(fa, gff, threads, output_file, file_num,compare_tools)

syn_pict.run_syn_pict(output_file)
distance_dup.run_distance(output_file,file_num,threads)