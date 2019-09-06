# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 14:43:13 2019

@author: wangpeng884112
"""

import argparse
import os
import pandas as pd
import numpy as np
import Levenshtein

from tqdm._tqdm import trange
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(description='Give the necessary parameters')
    parser.add_argument('--A_file', type=str, default='seq_exp_94.txt', help = 'A compared with B, this is the path to A sequence file')
    parser.add_argument('--B_file', type=str, default='StorzG_TSS_Table_LB_2.0.txt', help = 'A compared with B, this is the path to B sequence file')
    parser.add_argument('--mode', type=int, default = 1, help = '0: Give the similarity score of each sequence  1: Besides, give the similarity distribution')
    args = parser.parse_args()
    
    return args

def data_loader_default(args):
    A_file = os.path.join('data',args.A_file)
    B_file = os.path.join('data',args.B_file)
    A_str = pd.read_table(A_file, header = None)
    A_str = A_str[11::]
    B_str = pd.read_table(B_file)
    B_str = B_str[1::]
    
    return A_str, B_str

def cut_string(str_A,str_B):
    for i in range(0,len(str_A)):
        str_A[i] = str_A[i][5:-1]
    for i in range(0,len(str_B)):
        str_B[i] = str_B[i][0:-1]
        str_B[i] = str_B[i].upper()
    return str_A,str_B

def get_dissim_distribution(str_A,str_B):
    i = 0
    dissim_distribution = np.zeros((50,))
    for i in trange(len(str_A)):
        for j in range(0,len(str_B)):    
            result = Levenshtein.editops(str_A[i], str_B[j]) # catch the editting steps of two sequences
            
            record = []  # Get the different region
            for item in result:
                if item [0] == 'delete' or item[0] == 'replace' or (item[0] == 'insert' and item[2] == 50):
                    record.append(item[1])
                    record = list(set(record)) #remove repeating location
                
            tmp_dissim_region = []  # Get the dissimilar region in the str_A
            for k in range(0,50):
                if k in record:
                    tmp_dissim_region.append(k)
            for item in tmp_dissim_region:
                dissim_distribution[item] = dissim_distribution[item] + 1
        print(i)
    return dissim_distribution


## Achieve by others code
if __name__ == "__main__":
    args = parse_args()
    A_str, B_str = data_loader_default(args)
    str_A = A_str[0].tolist()
    str_B = B_str['promoter_seq']
    str_B = str_B.drop_duplicates(keep='last').tolist()
    str_A,str_B = cut_string(str_A,str_B)
    
    # Plot the distribution
    dissim_distribution = get_dissim_distribution(str_A,str_B)
    x = list(range(50))
    plt.plot(x,dissim_distribution)
    plt.ylim([0,1000000])
    plt.xlabel('Location')
    plt.ylabel('dissimilarity counts')


## Achieve by my own code
def edit_distance_matrix(str_A, str_B):
    matrix = np.zeros((len(str_A) + 1,len(str_B) + 1))
    for i in range(0,len(str_A) + 1):
        matrix[i,0] = i
    for i in range(0,len(str_B) + 1):
        matrix[0,i] = i
#    for i in range(1,len(matrix)):
#        for j in range(1,len(matrix[0])):
#            if str_A[i-1] == str_B[j-1]:
#                
    
    
    
    return matrix
