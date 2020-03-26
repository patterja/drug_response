#!/usr/local/bin/python3

#program: calculate_GR_at_serum_conc_using_GRmetrics.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Thu Aug 23 12:43:08 2017

#description
#calculate the Growth ratio at serum concentration, and add it the GRmetrics

#import libraries
import argparse
from gr50 import logistic
import pandas as pd
import numpy as np

def GR_of_serum(serum_con, params):
    params[1] =np.log10(params[1])
    if params[2] == 0:
        #set GRserum to 1.0
        GRserum=1.0
    else:
        print (serum_con, params)
        GRserum = logistic(serum_con, params) 
    return GRserum

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--m', help='input the GRmetrics')
    parser.add_argument('--s', help='input the serum concentration file')
    parser.add_argument('--o', help='output the new matrix of GRmetrics')
    
    args = parser.parse_args()

    #read in the table of serum concentration of drugs
    df_serum_conc=pd.read_csv(args.s, sep="\t", index_col=0, header=0)

    #read in the GRmetrics
    df_GRmetrics=pd.read_csv(args.m, sep="\t", header=0, index_col=None)
    for i in range(len(df_GRmetrics.index)):
        lst_params=df_GRmetrics.ix[i, ["GRinf", 'GEC50', 'GR Hill Coefficient']].tolist()
        s_drug=df_GRmetrics.ix[i,'Small Molecule Name']
        print (df_serum_conc.ix[s_drug])
        serum_con=df_serum_conc.ix[s_drug].tolist()[0]
        GRserum = GR_of_serum(serum_con, lst_params)
        GRserum = ('%.4f' % GRserum)
        df_GRmetrics.ix[i,"GRserum"]=GRserum
        #df_GRmetrics.ix[i,'GRserum']=GRserum
    #write down the new GRmetrics
    df_GRmetrics.to_csv(args.o, sep="\t", index=False)

