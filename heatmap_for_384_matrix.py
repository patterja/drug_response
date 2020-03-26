#program: generate_384_well_layout_table.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Wed Mar 15 10:10:35 2017

#description

#import libraries
import argparse
import pandas as pd
import csv
import seaborn as sns
import matplotlib.pyplot as plt

def generate_table(s_data, s_layout):
    df_data=pd.read_csv(s_data, sep="\t", index_col=1)
    df_layout=pd.read_csv(s_layout, sep="\t", index_col=1)
    df_merge=pd.merge(df_layout, df_data, left_index=True, right_index=True, how='outer')
    df_merge=df_merge[['PlateIndex', 'Objects']]
    df_merge=df_merge.set_index('PlateIndex')
    df_merge.columns=["CellCount"]
    d_cellcount=df_merge.to_dict(orient='index')
    d_data_layout={}
    i_row=0
    for i in range(384):
        i_num=i+1
        i_rem=i_num%24
        if (i_rem==1):
            i_row+=1
            d_data_layout.update({i_row: [d_cellcount[i_num]["CellCount"]]})
        else:
            d_data_layout[i_row].append(d_cellcount[i_num]["CellCount"])
    #convert the dictionary into dataframe
    df_data_layout=pd.DataFrame.from_dict(d_data_layout, orient='index')
    df_data_layout.columns=range(1,25)
    df_data_layout.index=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    #return the table
    return df_data_layout

def sns_heatmap(df_data):
    sns.set()
    fig=plt.figure(figsize=(18,8))
    ax=fig.add_subplot(1,1,1)
    sns.heatmap(df_data, ax=ax, annot=True, fmt='.0f', annot_kws={"size":7}, linewidths=.5)
    return fig

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the well level cell count')
    parser.add_argument('--l', help='input the layout.txt file')
    parser.add_argument('--h', help='output the heatmap figure')
    
    args = parser.parse_args()

    #convert the data into 384 well layout table
    df_table=generate_table(s_data=args.i, s_layout=args.l)
    #plot the heatmap
    hm_fig=sns_heatmap(df_table)
    hm_fig.savefig(args.h, bbox_inches='tight')
