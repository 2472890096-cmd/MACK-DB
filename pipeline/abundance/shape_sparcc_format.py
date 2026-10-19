#!/usr/bin/env python3

import os
import sys
import argparse

def read_params():
    parser = argparse.ArgumentParser(description="merge all sample's species into a table")
    # add arguments
    parser.add_argument("-l", "--list", action="store", help="Input a file list")
    parser.add_argument("-i", "--input_dir", action="store", help="Input a file dir")
    parser.add_argument("-f", "--filter", type=float, default=0.0, action="store", help="filter low abundance species, default: 0.0")
    parser.add_argument("-p", "--profile", action="store", help="output profile file name")
    parser.add_argument("-s", "--species", action="store", help="output species file name")
    args = parser.parse_args()
    return args

args = read_params()
run_li = open(args.list, "rt")
file_dir = args.input_dir

all_ncbi_id_dir = dict()
all_run_profile = dict()
n = 0
for run in run_li:
    n += 1
    run = run.rstrip()
    file_name = file_dir + "/" + run + "_profile.txt"
    profile_dic = dict()
    profile_f = open(file_name, "rt")
    for i in profile_f:
        li = i.rstrip().split("\t")
        ncbi_id = li[0]
        num = li[1]
        profile_dic[ncbi_id] = num
        if ncbi_id not in all_ncbi_id_dir:
            all_ncbi_id_dir[ncbi_id] = 1
        else:
            all_ncbi_id_dir[ncbi_id] += 1
    all_run_profile[run] = profile_dic
    profile_f.close()
run_li.close()


screen_ncbi_id_li = list()
filter_ratio = args.filter
screen_ratio = round(n * filter_ratio)
for i, j in all_ncbi_id_dir.items():
    if j >= screen_ratio:
        screen_ncbi_id_li.append(i)

#print("ncbi_id_id" + "\t" + "\t".join(list(all_run_profile.keys())))
profile_n = args.profile
profile_f = open(profile_n, "w")

profile_f.write("OTU_id" + "\t" + "\t".join([str(i) for i in range(n)]) + "\n")
for key, i in enumerate(screen_ncbi_id_li):
    profile_f.write(str(key))
    for run, profile_dic in all_run_profile.items():
        if i in profile_dic:
            profile_f.write("\t" + profile_dic[i])
        else:
            profile_f.write("\t"+"0")
    profile_f.write("\n")
profile_f.close()

ncbi_id_n = args.species
ncbi_id_f = open(ncbi_id_n, 'w')
for key, i in enumerate(screen_ncbi_id_li):
        ncbi_id_f.write(str(key) + "\t" + i + "\n")
ncbi_id_f.close()
