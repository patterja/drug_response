#program: boxplot_for_GRmax_and_serum_conc.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Tue Apr 11 10:55:37 2017

#edit by: Janice Patterson
#edit on: 4/30/2020

#description

#import libraries
import argparse
import pandas as pd
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import math


parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', help='input the GRmetrics_GRserum_added.txt')
parser.add_argument('--c', help='input the cell line class and color')
parser.add_argument('--w', help='input the work directory ')

args = parser.parse_args()

#read in the GRmetrics.txt
df_data=pd.read_csv(args.i, sep="\t")
df_data=df_data.set_index('cell_line')
df_data.dropna(subset = ["GRserum"], inplace=True)

#read in cell line class and color
if args.c != None:
    df_color=pd.read_csv(args.c, sep="\t", index_col=0)
    df_merge = pd.merge(df_data, df_color, left_index=True, right_index=True)

else:
    col=["blue", "orange", "green", "red", "purple", "yellow", "brown", "gray", "pink"]

    #df_color = pd.DataFrame({'Class': df_data['Class'], 'color': [None]* len(df_data['Class'])})

    class_cond = []
    for subtype in df_data['Class'].unique():
        class_cond.append(df_data['Class'] == subtype)
    df_merge = df_data
    df_merge['Color'] = np.select(class_cond, col[0:len(df_data['Class'].unique())])




#combine these to matrix
#extract the sub-matrix
df_GRmax=df_merge[['agent', 'GRmax', 'GRserum', 'Class','Color']]
df_GRmax.index.name='cell_line'
df_GRmax=df_GRmax.reset_index()
df_GRmax=df_GRmax.set_index('agent')

lst_agents=df_GRmax.index.unique()
for s_agent in lst_agents:
    df_agent=df_GRmax.ix[s_agent]
    df_agent=df_agent.set_index('cell_line')
    df_agent=df_agent.sort_values('Class')
    ##barplot with the GRmax, and add the drug concentration in the serum
    fig=plt.figure(figsize=(16,8))
    ax=fig.add_subplot(1,1,1)
    fig1=df_agent[['GRmax', "GRserum"]].plot(kind='bar', width=0.7, color=[df_agent['Color'], df_agent['Color']], title=s_agent, ax=ax)
    fig1.set_xlabel("Cell Lines")
    fig1.set_ylabel("Growth Rate")
    fig1.set_ylim(-1,1.5) 
    #add hatch for each bar
    bars = ax.patches
    hatches = ''.join(h*len(df_agent) for h in '/.')
    for bar, hatch in zip(bars, hatches):
        bar.set_hatch(hatch)
    ######add line at 1.0 on y axis######
    cc_len=len(df_agent)
    fig2=plt.plot([-1, cc_len], [1.0, 1.0], color="black", linestyle="--", lw=1)  
    
    #add legend of cell line class
    df_class_color=df_agent[['Class', 'Color']].drop_duplicates()
    df_class_color=df_class_color.set_index('Class')
    lst_handles=[]
    for s_class in df_class_color.index:
        #if all the cell lines in this class didn't reach the GRmax, don't plot the class legend
        ds_GRmax=df_agent[df_agent['Class']==s_class]['GRmax']
        if (len(ds_GRmax.iloc[ds_GRmax.to_numpy().nonzero()[0]]) >0):
            patch=mpatches.Patch(color=df_class_color.ix[s_class, 'Color'], label=s_class)
            lst_handles.append(patch)
    #if no legend, don't plot
    if(len(lst_handles)>0):
        #add legend for GRmax and GRserum
        handle_GRmax=mpatches.Patch(facecolor="white", edgecolor="black", label='GRmax', hatch ='////')
        handle_GRserum=mpatches.Patch(facecolor="white", edgecolor="black", label='GRserum', hatch ='oo')
        lst_handles=[handle_GRmax, handle_GRserum] + lst_handles
        plt.legend(handles=lst_handles, prop={'size':8})
    else:
        print ('No cell line reach the GRmax for this compound')
    #save figure
    fig_name=args.w + '/' + s_agent + "_GRmax_and_GRserum.pdf"
    fig.savefig(fig_name, bbox_inches='tight')
