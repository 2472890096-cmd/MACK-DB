# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 09:48:09 2024

@author: Administrator
"""
import argparse
import os

def read_params():
    parser = argparse.ArgumentParser(usage='usage1:%(prog)s -l list.txt -d ./ -o ./links/\n\
       usage2:%(prog)s -l list.txt -d ./ -o ./links/ -s T\n\
       usage3:%(prog)s -l list.txt -d ./ -o ./links/ -s T -f 0.85',
                                    description="Find all bioproject ID corresponding to a link")
    # add arguments
    parser.add_argument("-l", "--list", default=None, help="Bioproject ID that needs to be analyzed. Requires '.txt' format with one ID per line")
    parser.add_argument("-d", "--dir", default="./", help="The path where the file is located")
    parser.add_argument("-o", "--out", default="./", help="The folder where the output results are stored")
    parser.add_argument("-n", "--name", default="link_belonging", help="output file name")
    parser.add_argument("-f", "--filter", type=float, default=2, help="Add a filtering target: \nIf 't' less than 1, Retain links that have appeared in at least 't'*100%% of all projects. \nIf parameter greater than 1, Retain links that have appeared in at least 't' projects. ")
    parser.add_argument("-s", "--switch", default="F", help="The switch to control whether to filter or not: 'T' to filter, 'F' to pass.")
    parse = parser.parse_args()
    return parse

def get_dic(parse):
    input_abs_path = os.path.abspath(parse.dir)
    # projectname_file = os.path.join(input_abs_path, parse.list)
    projectname_file = os.path.abspath(parse.list)
    dic_link_ID = {}
       
    with open( projectname_file, "rt") as ID_file:
        projectname_list = [ line.rstrip() for line in ID_file ]
    
    for i in projectname_list:
        input_file = os.path.join(input_abs_path, i + "_links.txt")
        if os.path.exists(input_file):
            with open(input_file, "rt") as file:
                for j in file.readlines():
                    if j.rstrip() not in dic_link_ID.keys():
                        dic_link_ID[ j.rstrip() ] = [i]
                    else:
                        dic_link_ID[ j.rstrip() ].append(i)
        else:
            print("----------------\n\033[31mWARNING\033[0m: '"+ input_file +"' is NOT exist.\n----------------")
            continue
    return dic_link_ID, len(projectname_list)

def dic_filter(dic, num, fil):
    if fil <= 1:
        filter_dic = {k:v for k, v in dic.items() if len(v) >= num*fil}
    else:
        filter_dic = {k:v for k, v in dic.items() if len(v) >= fil}
    return filter_dic

if __name__=="__main__":
    parse = read_params()
    dic, project_count = get_dic(parse)
    
    output_abs_path = os.path.abspath(parse.out)
    os.makedirs(output_abs_path, exist_ok=True)
    output_file = os.path.join(output_abs_path, parse.name)
    
    with open( output_file + ".belonging", 'w', encoding='utf-8') as f1:
        f1.write("species_a\tspecies_b\tproject_count\tproject_ID\n")
        for k, v_list in dic.items():
                f1.write(k+ "\t"+ str(len(v_list))+ "\t"+ ';'.join(v_list)+ "\n")
    
    if parse.switch == "T":
        print("----------------\n\033[33mNOTICE\033[0m: Filter Criteria is to retain only links that have appeared in at least "+ str(parse.filter*100) + "% of all projects\n----------------")
        filter_dic = dic_filter(dic, project_count, parse.filter)
        with open( output_file + "_filtered.belonging", 'w', encoding='utf-8') as f2:
            f2.write("species_a\tspecies_b\tproject_count\tproject_ID\n")
            for k, v_list in filter_dic.items():
                    f2.write(k+ "\t"+ str(len(v_list))+ "\t"+ ';'.join(v_list)+ "\n")
    else:
        print("----------------\n\033[33mNOTICE\033[0m: The results were not filtered this time, and the filtered results will not be output.\nIf filtering is required, please specify the '-s' parameter as 'T' and set the '-f' parameter according to the requirements.\n----------------")
    

    print("\033[32mUnfiltered results has been output to:\033[0m" + output_file + ".txt\n")
    print("\033[32mfiltered results has been output to:\033[0m" + output_file + "_filtered.txt\n")
