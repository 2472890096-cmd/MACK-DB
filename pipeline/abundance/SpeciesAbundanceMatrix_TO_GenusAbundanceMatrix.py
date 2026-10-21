# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 17:32:08 2024

@author: Administrator
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from functools import reduce

def read_params():
    parser = argparse.ArgumentParser(usage='usage1:%(prog)s -l list.txt', description="")
    # add arguments
    parser.add_argument("-i", "--input", help="")
    parser.add_argument("-t", "--tool", help="")
    parser.add_argument("-d", "--dir", help="")
    parser.add_argument("-o", "--outdir", default='./', help="")
    # parser.add_argument("-d", "--dir", help="")
    parse = parser.parse_args()
    return parse

def get_proj_genus_sp_dic(args):
    with open(args.input, 'r') as f:
        proj_genus_sp_dic = {line.strip():{} for line in f}
    for proj in proj_genus_sp_dic:
        try:
            with open(os.path.join(os.path.abspath(args.dir), proj + '_species.txt'), 'r') as f:
                for line in f:
                    genus = line.strip().split('\t')[1].split(' ')[0].split('_')[1]
                    if genus in proj_genus_sp_dic[proj]:
                        proj_genus_sp_dic[proj][genus].append(line.strip().split('\t')[1])
                    else:
                        proj_genus_sp_dic[proj][genus] = [line.strip().split('\t')[1]]
        except FileNotFoundError:
            print('FileNotFoundError: '+proj)
    return proj_genus_sp_dic

def get_proj_sp_num_dic(args):
    with open(args.input, 'r') as f:
        proj_sp_num_dic = {line.strip():{} for line in f}
    for proj in proj_sp_num_dic:
        try:
            with open(os.path.join(os.path.abspath(args.dir), proj + '_species.txt'), 'r') as f:
                f_li = f.readlines()
                for i in range(len(f_li)):
                    proj_sp_num_dic[proj][f_li[i].strip().split('\t')[1]] = str(i)
        except FileNotFoundError:
            print('FileNotFoundError: '+proj)
    return proj_sp_num_dic

def get_sp_abundance(args):
    with open(args.input, 'r') as f:
        proj_num_abundance_dic = {line.strip():{} for line in f}
    for proj in proj_num_abundance_dic:
        with open(os.path.join(os.path.abspath(args.dir), proj + '_profile.txt'), 'r') as f:
            f.readline()
            for line in f:
                proj_num_abundance_dic[proj][line.split('\t')[0]] = list(map(float, line.strip().split('\t')[1:]))
                # print(proj_num_abundance_dic[proj][line.split('\t')[0]])
    return proj_num_abundance_dic
    
def merge_abundance(args, proj_genus_sp_dic, proj_sp_num_dic, proj_num_abundance_dic):
    for proj, genus_dic in proj_genus_sp_dic.items():
        with open(os.path.join(os.path.abspath(args.dir), proj + '_profile.txt'), 'r') as f:
            title = f.readline()
        f1 = open(os.path.join(os.path.abspath(args.outdir), 'genus', proj + '_genus.txt'), 'w')
        f2 = open(os.path.join(os.path.abspath(args.outdir), 'abundance', proj + '_abundance.txt'), 'w')
        f2.write(title)
        genus_num = 0
        for genus, sp_li in genus_dic.items():
            num_li = [proj_sp_num_dic[proj][sp] for sp in sp_li]
            sp_abundance_dic = {num:proj_num_abundance_dic[proj][num] for num in num_li}
            genus_abundance_li = [sum(x) for x in zip(*sp_abundance_dic.values())]
            f1.write(str(genus_num) + '\t' + genus + '\n')
            f2.write(str(genus_num) + '\t' + '\t'.join([format(num, '.4f') for num in genus_abundance_li]) + '\n')
            genus_num = genus_num + 1
        f1.close()
        f2.close()
        
if __name__=="__main__":
    args = read_params()
    proj_genus_sp_dic = get_proj_genus_sp_dic(args)
    proj_sp_num_dic = get_proj_sp_num_dic(args)
    proj_num_abundance_dic = get_sp_abundance(args)
    merge_abundance(args, proj_genus_sp_dic, proj_sp_num_dic, proj_num_abundance_dic)
