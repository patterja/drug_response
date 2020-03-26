#program: add_max_drug_does_to_GRmetrics.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Tue Apr 18 15:50:40 2017

#description

#import libraries
import argparse
import pandas as pd

def add_max_conc_to_GRm(f_GRm, f_layout):
    df_layout=pd.read_csv(f_layout, sep="\t")
    df_conc=df_layout[['Compound', 'Concentration']]
    df_conc_max=df_conc.groupby('Compound').max()
    df_conc_max.index.name=None,
    df_conc_max.columns=['Maximal Does']
    #set the index to agent for df_GRm
    df_GRm=pd.read_csv(f_GRm, sep="\t")
    df_GRm=df_GRm.set_index('agent')
    #merge these two matrixs
    df_merge=pd.merge(df_GRm, df_conc_max, left_index=True, right_index=True)
    df_merge.index.name='agent'
    return df_merge

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--g', help='input the GRmetrics.txt')
    parser.add_argument('--l', help='input the layout.txt')
    parser.add_argument('--o', help='output the edited GRmetrics.txt')
    
    args = parser.parse_args()
    
    #read in the GRmetrics.txt
    df_GRm_edit=add_max_conc_to_GRm(f_GRm=args.g, f_layout=args.l)
    df_GRm_edit.to_csv(args.o, sep="\t")
