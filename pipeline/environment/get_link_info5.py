# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:25:09 2024

@author: Administrator
"""
import argparse
import os
import sys
import pandas as pd
from collections import defaultdict

def read_params():
    parser = argparse.ArgumentParser(description="Find all bioproject ID corresponding to a link")
    # add arguments
    parser.add_argument("-i", "--input", help="Use a file with '.belonging' format as the input file.")
    parser.add_argument("-t", "--tool", default="fastspar", help="Choose one of the following three options: fastspar; ccl_coun; ccl_freq.")
    # parser.add_argument("-e", "--env", default=None, help="macrogenomic environment. e.g. marine")
    parser.add_argument("-o", "--out", default="./", help="The folder where the output results are stored")
    parser.add_argument("-n", "--name", default="link_info.txt", help="output file name")
    parser.add_argument("-s", "--source", default="data/projects.source", help="Correspondence between the project and the environment.")
    parser.add_argument("-c", "--corr-dir", default="work/fastspar/corr_pv_csv", help="Directory containing PROJECT_corr_pvalue.csv files.")
    parser.add_argument("-a", "--abundance-dir", default="data/species/abundance", help="Directory containing PROJECT_abundance.matrix files.")
    args = parser.parse_args()
    return args

def get_dic(args):
    if args.input.endswith('belonging') == False:
        print("----------------\n\033[33mERROR\033[0m: '-i' parameter setting must be 'link_belonging_filtered.txt' or 'link_belonging.txt'. Program terminated.\n----------------")
        sys.exit(0)
    link_file = os.path.abspath(args.input)
    data_dict = defaultdict(dict)
    with open(link_file, "r") as file:
        file.readline()
        for line in file:
            species_pair, project_count, project_ids = line.strip().rsplit("\t",2)
            data_dict[species_pair] = project_ids.split(";")
            data_dict.pop('NA\tNA', None)
    return data_dict
    
def get_info(dic, args):
    info_dict = defaultdict(dict)
    
    # 1. 找出这个大环境下所有涉及到的独特 Project，防止重复读文件
    unique_projs = set()
    for proj_list in dic.values():
        unique_projs.update(proj_list)
        
    print(f"  [优化加速] 正在将 {len(unique_projs)} 个 Project 的数据载入内存缓存...")
    
    sam_cache = {}
    proj_data_cache = defaultdict(dict)
    
    # 2. 每个 Project 的文件只读取 1 次！
    for i in unique_projs:
        # 极速获取样本数（只读第一行，彻底抛弃巨慢的 pd.read_csv）
        profile_f = os.path.join(args.abundance_dir, i + '_abundance.matrix' )
        try:
            with open(profile_f, 'r') as f:
                header = f.readline()
                sam_cache[i] = str(len(header.strip().split('\t')) - 1)
        except:
            sam_cache[i] = "NA"
            
        # 极速读取 CSV 并存入字典
        if args.tool == 'fastspar':
            info_file = os.path.join(args.corr_dir, i + '_corr_pvalue.csv' )
        
        try:
            with open(info_file, "r") as file:
                file.readline() # 跳过表头
                for line in file:
                    species_pair, cor, pv = line.strip().rsplit(",", 2)
                    # 将 "物种A,物种B" 变成 frozenset，方便后续 O(1) 极速匹配
                    sp_set = frozenset(species_pair.split(","))
                    proj_data_cache[i][sp_set] = (round(float(cor), 5), round(float(pv), 5))
        except Exception as e:
            pass
            
    print("  [优化加速] 缓存完毕！开始极速匹配边数据...")
    
    # 3. O(1) 极速匹配，不再读硬盘
    for pair, proj_list in dic.items():
        cor_li = []
        pv_li = []
        sam_li = []
        
        target_set = frozenset(pair.split("\t"))
        
        for i in proj_list:
            sam_li.append(sam_cache.get(i, "NA"))
            
            # 直接从内存字典里取数据
            res = proj_data_cache[i].get(target_set)
            if res:
                cor_li.append(res[0])
                pv_li.append(res[1])
                
        if len(cor_li) > 0:
            aver_cor = str(round(sum(cor_li) / len(cor_li), 5))
            aver_pv = str(round(sum(pv_li) / len(pv_li), 5))
        else:
            aver_cor, aver_pv = "NA", "NA"
            
        cor_l = [str(num) for num in cor_li]
        pv_l = [str(num) for num in pv_li]
        info_dict[pair] = [aver_cor, aver_pv, cor_l, pv_l, sam_li]
        
    return info_dict
# def get_environment(args):
    
#     env = os.path.basename(os.path.dirname(os.path.dirname((os.path.abspath(args.input)))))
#     env_li = env.split("_")
#     if len(env_li) == 2:
#         sub_env = env
#     elif len(env_li) == 1:
#         sub_env = '-'
#     else:
#         print("----------------\n\033[33mERROR\033[0m: Environmental information is not compliant. Program terminated.\n----------------")
#         sys.exit(0)
#     super_env = env_li[0]
    
#     return super_env, sub_env

def get_environment(dic, args):
    env_dic = {k:[] for k in dic}
    with open(args.source, "r") as file:
        file.readline()
        for line in file:
            env, projects = line.strip().split("\t")
            for pair, proj_list in dic.items():
                for proj in proj_list:
                    # print(proj)
                    # print(projects)
                    # print(str(proj in projects))
                    if proj in projects:
                        env_dic[pair].append(env)

    env_dic_nodup = {k: list(set(v)) for k, v in env_dic.items()}
    del env_dic
    super_env_dic = defaultdict(dict)
    sub_env_dic = defaultdict(dict)
    
    for i,j in env_dic_nodup.items():
        if len(j) > 1:
            super_env_dic[i] = []
            sub_env_dic[i] = []
            for a in j:
                if len(a.split("_")) == 2:
                    sub_env_dic[i].append(a)
                elif len(a.split("_")) == 1:
                    sub_env_dic[i].append('-')
                else:
                    print("----------------\n\033[33mERROR\033[0m: Environmental information is not compliant. Program terminated.\n----------------")
                    sys.exit(0)
                super_env_dic[i].append(a.split("_")[0])
            super_env_dic[i] = ';'.join(super_env_dic[i])
            sub_env_dic[i] = ';'.join(sub_env_dic[i])
            
        elif len(j) == 1:
            if len(j[0].split("_")) == 2:
                sub_env_dic[i] = j[0]
            elif len(j[0].split("_")) == 1:
                sub_env_dic[i] = '-'
            else :
                print("----------------\n\033[33mERROR\033[0m: Environmental information is not compliant. Program terminated.\n----------------")
                sys.exit(0)
            super_env_dic[i] = j[0].split("_")[0]
            
    return super_env_dic, sub_env_dic
    
        
        
if __name__=="__main__":
    args = read_params()
    dic = get_dic(args)
    info_dict = get_info(dic, args)
    
    # super_env, sub_env = get_environment(args)
    # with open( os.path.join(os.path.abspath(args.out), args.name), 'w', encoding='utf-8') as f:
    #     f.write("species_a\tspecies_b\taver_corr\taver_pvalue\tproject_count\tproject_ID\tsuper_env\tsub_env\n")
    #     for i in dic:
    #         f.write(i+ "\t"+ '\t'.join(info_dict[i])+ "\t"+ str(len(dic[i]))+ "\t"+ ';'.join(dic[i])+ "\t"+ super_env+ "\t"+ sub_env+ "\n")
        
    # super_env_dic, sub_env_dic = get_environment(dic, args)
    # 修复版环境名称提取逻辑
    base_name = os.path.basename(args.input).replace('_core_links_filtered.belonging', '').replace('_core_links.belonging', '')
    env_parts = base_name.split('_')
    
    if len(env_parts) == 2:
        super_env = env_parts[0]
        sub_env = env_parts[1]
    elif len(env_parts) == 1:
        super_env = env_parts[0]
        sub_env = '-'
    else:
        super_env = base_name
        sub_env = 'NA'
    with open( os.path.join(os.path.abspath(args.out), args.name), 'w', encoding='utf-8') as f:
        f.write("species_a\tspecies_b\taver_corr\taver_pv\tproj_count\tproj_ID\tproj_corr\tproj_pv\tproj_sam\tsuper_env\tsub_env\n")
        for i in dic:
            # f.write(i+ "\t"+ '\t'.join(info_dict[i][0:2])+ "\t"+ str(len(dic[i]))+ "\t"+ ';'.join(dic[i])+ "\t"+ ';'.join(info_dict[i][2])+ "\t"+ ';'.join(info_dict[i][3])+ "\t"+ ';'.join(info_dict[i][4])+ "\t"+ super_env_dic[i]+ "\t"+ sub_env_dic[i]+ "\n")
            f.write(i+ "\t"+ '\t'.join(info_dict[i][0:2])+ "\t"+ str(len(dic[i]))+ "\t"+ ';'.join(dic[i])+ "\t"+ ';'.join(info_dict[i][2])+ "\t"+ ';'.join(info_dict[i][3])+ "\t"+ ';'.join(info_dict[i][4])+ "\t"+ super_env+ "\t"+ sub_env+ "\n")

                        
