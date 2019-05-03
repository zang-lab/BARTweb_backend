# Time-stamp: <2017-08-10>
'''Module for calculating ROC-AUC values for all TF datasets

Copyright (c) 2017, 2018 Chongzhi Zang, Zhenjia Wang <zhenjia@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: $Id$
@author: Chongzhi Zang, Zhenjia Wang
@contact: zhenjia@virginia.edu

'''

import os,sys,time
import json
import multiprocessing # multiprocessing on dealing with TF datasets

from revised_bart import StatTest

def get_tf_file_data(tf_json):
    starttime = time.time()
    with open(tf_json, 'r') as fm:
        tf_file_map = json.load(fm)
    endtime = time.time()
    print("===loading tf file mapping list: {} seconds".format(endtime-starttime))
    return tf_file_map
    
# 84M memory taken
def get_matrix_data(overlap_json):
    starttime = time.time()
    with open(overlap_json, 'r') as fm:
        matrix_data = json.load(fm)
    endtime = time.time()
    print("===loading matrix file: {} seconds".format(endtime-starttime))
    return matrix_data

def cal_auc_for_all_tfs(positions, matrix_data, tf_file_len):
    udhs_len = len(positions)
    groupsize = 10000   # 2.7M / 10000 = 272 steps
    groups = int(udhs_len/groupsize)

    print('total groups: ' + str(groups))

    dict_tf_auc = {} # each tf file auc
    each_tf_t = {} # each tf file, how many 1s

    for i in range(1, tf_file_len+1): # initiate the axis for each tf file
        each_tf_t[i] = []
        dict_tf_auc[i] = 0.0
    
    print("===parsing data in group size: {}".format(groupsize))
    for group_id in range(groups+1):
        count_tf_t = {}
        for i in range(1, tf_file_len+1):
            count_tf_t[i] = 0

        for i in range(group_id*groupsize, min((group_id+1)*groupsize, udhs_len)):
            for tf in matrix_data[str(positions[i])].strip().split(): # count how many 1s in this group
                count_tf_t[int(tf)] += 1

        for tf, count_t in count_tf_t.items():
            if each_tf_t[tf]:
                each_tf_t[tf].append(each_tf_t[tf][-1]+count_t)
            else:
                each_tf_t[tf].append(count_t)

    for key, tf_t_value in each_tf_t.items():
        for i in range(len(tf_t_value)):
            cur_x = (min(groupsize*(i+1), udhs_len)-tf_t_value[i])/(udhs_len-tf_t_value[-1])
            cur_y = tf_t_value[i]/tf_t_value[-1]
            if i == 0:
                width = cur_x
                height = cur_y/2
                dict_tf_auc[key] += height*width
            else:
                width = cur_x - (min(groupsize*(i), udhs_len)-tf_t_value[i-1])/(udhs_len-tf_t_value[-1])
                height = (cur_y + tf_t_value[i-1]/tf_t_value[-1])/2
                dict_tf_auc[key] += height*width

    # TODO: ========below is for AUPR =====
    # for key, tf_t_value in each_tf_t.items():
    #     for i in range(len(tf_t_value)):
    #         cur_x = tf_t_value[i]/tf_t_value[-1]
    #         cur_y = tf_t_value[i]/min(groupsize*(i+1), udhs_len)
    #         if i == 0:
    #             width = cur_x
    #             height = cur_y
    #             dict_tf_auc[key] += height*width
    #         else:
    #             width = cur_x - tf_t_value[i-1]/tf_t_value[-1]
    #             height = cur_y
    #             dict_tf_auc[key] += height*width
    # ==========  end modification =========

    return dict_tf_auc

def get_position_list(enhancerfile):
    '''
    Get the ID list of DHS, according to the decreasingly sorted scores in MARGE enhancer profile
    ''' 
    fin = open(enhancerfile,'rb')
    line = fin.readline()  
    score = {}
    while line:
        line = line.strip().split()
        try:
            score[line[-2]]=float(line[-1])
        except:
            pass
        line = fin.readline()
    fin.close()
    return sorted(score.keys(),key=score.get,reverse=True)

def cal_auc(enhancerfile, args):
    tf_json = args.tffile
    overlap_json = args.tfoverlap
    output_name = args.ofilename
    normfile = args.normfile

    auc_file = output_name + '_auc.txt'
    stat_file = output_name + '_bart_results.txt'
    
    positions = get_position_list(enhancerfile)
    positions = [int(i) for i in positions]
    if len(positions) == 0:
        sys.stderr.write('Input file might not with right format!\n')
        sys.exit(1)
    sys.stdout.write("Calculating ROC-AUC values for all transcription factors:\n\n")

    tf_dict = get_tf_file_data(tf_json)
    overlap_dict = get_matrix_data(overlap_json)
    tf_auc = cal_auc_for_all_tfs(positions, overlap_dict, len(tf_dict))

    # output file of AUC-ROC values for all TFs
    with open(auc_file, 'w') as aucf:
        for key in sorted(tf_auc.keys(),key=tf_auc.get,reverse=True):
            tf_file_name = tf_dict[str(key)]
            aucf.write('{}\tAUC = {:.3f}\n'.format(tf_file_name, tf_auc[key]))
    print('\n--ROC-AUC calculation finished!\n--Results saved in file: {}\n'.format(auc_file))
    # StatTest.stat_test(AUCs,args)
    StatTest.stat_test(tf_auc, tf_dict, stat_file, normfile)

if __name__ == '__main__':
    margefile = '/nv/vol190/zanglab/wm9tr/software/BART-v1.0.1-py3-full/BART/hg38_library/hg38_test_data/hg38_AR_up_genes_enhancer_prediction.txt'
    tf_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_file_str_format.json'
    overlap_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_matrix_str_format.json'

    # tf_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_UAC_file_str_format.json'
    # overlap_json = '/nv/vol190/zanglab/wm9tr/test/test_bart_storage/human_UAC_matrix_str_format.json'

    output_name = 'test_bart/AR_aupr'
    norm_file = '/nv/vol190/zanglab/wm9tr/software/BART-v1.0.1-py3-full/BART/hg38_library/hg38_MSigDB.dat'

    cal_auc(margefile, tf_json, overlap_json, output_name, norm_file)
