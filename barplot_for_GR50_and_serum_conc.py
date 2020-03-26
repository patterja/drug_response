#program: boxplot_for_GR50_and_serum_conc.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Tue Apr 11 10:55:37 2017

#description

#import libraries
import argparse
import pandas as pd
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import math

parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', help='input the GRmetrics.txt')
parser.add_argument('--c', help='input the cell line class and color')
parser.add_argument('--s', help='input the drug_serum_conc.txt')
parser.add_argument('--w', help='input the work directory ')

args = parser.parse_args()

#read in the serum concentration
#Drug list 	Serum concentration
#BIBW2992	12.74
#Vorinostat	1400-2200
d_serum={}
with open (args.s, 'r') as tsv_serum:
    h_reader=csv.reader(tsv_serum, delimiter="\t")
    s_header=next(h_reader)
    for lst_row in h_reader:
        (s_drug, s_conc) = lst_row
        lst_conc=s_conc.split("-")
        mean_conc=np.mean([float (x) for x in lst_conc])
        try:
            d_serum[s_drug]=mean_conc
        except KeyError:
            d_serum.update({s_drug: mean_conc})

#read in the GRmetrics.txt
df_data=pd.read_csv(args.i, sep="\t")
df_data=df_data.set_index('cell_line')
#read in cell line class and color
df_color=pd.read_csv(args.c, sep="\t", index_col=0)
#combine these to matrix
df_merge=pd.merge(df_data, df_color, left_index=True, right_index=True)
#extract the sub-matrix
df_gr50=df_merge[['agent', 'GR50', 'Class', 'Maximal Does', 'Color']]
df_gr50.index.name='cell_line'
df_gr50=df_gr50.reset_index()
df_gr50=df_gr50.set_index('agent')
#replace the infinate value into the 10*max does
for i in range(len(df_gr50)):
    if(abs(df_gr50.ix[i, 'GR50']) == np.inf):
        df_gr50.ix[i, 'GR50']=df_gr50.ix[i, 'Maximal Does']
#calculate the log10(GR50) 
df_gr50['GR50']=np.log10(df_gr50['GR50'])

lst_agents=df_gr50.index.unique()
for s_agent in lst_agents:
    df_agent=df_gr50.ix[s_agent]
    df_agent=df_agent.set_index('cell_line')
    df_agent=df_agent.sort_values('Class')
    ##barplot with the GR50, and add the drug concentration in the serum
    fig=plt.figure(figsize=(16,8))
    ax=fig.add_subplot(1,1,1)
    fig1=df_agent['GR50'].plot(kind='bar', color=df_agent['Color'], title=s_agent, ax=ax)
    fig1.set_xlabel("Cell Lines")
    fig1.set_ylabel("Log10(GR50) nM")
    
    ##line plot for the serum concentration
    #split the combine compounds
    lst_agent=s_agent.split('+')
    #calculate the mean serum concentration of compound
    lst_conc=[]
    for agent in lst_agent:
        try:
            lst_conc.append(d_serum[agent])
        except KeyError:
            print ('%s is does not have the serum concentration' % (agent))
            next
    if (len(lst_conc)>0):
        mean_conc=np.mean([float(x) for x in lst_conc])
        f_serum=np.log10(float(mean_conc))
        cc_len=len(df_agent)
        #line plot for the serum conc
        fig2=plt.plot([-1, cc_len], [f_serum, f_serum], color="black", linestyle="--")
    #plt.ylim(math.floor(min(df_agent['GR50'])), math.ceil(max(df_agent['GR50'])) + 1)
    #add legend of cell line class
    df_class_color=df_agent[['Class', 'Color']].drop_duplicates()
    df_class_color=df_class_color.set_index('Class')
    lst_handles=[]
    for s_class in df_class_color.index:
        #if all the cell lines in this class didn't reach the GR50, don't plot the class legend
        ds_gr50=df_agent[df_agent['Class']==s_class]['GR50']
        if (len(ds_gr50.iloc[ds_gr50.nonzero()[0]]) >0):
            patch=mpatches.Patch(color=df_class_color.ix[s_class, 'Color'], label=s_class)
            lst_handles.append(patch)
    #if no legend, don't plot
    if(len(lst_handles)>0):
        serum_line = mlines.Line2D([], [], color='black', linestyle="--",
                          label='serum concentration')
        lst_handles.append(serum_line)
        plt.legend(handles=lst_handles, prop={'size':8})
    else:
        print ('No cell line reach the GR50 for this compound')
    #save figure
    fig_name=args.w + '/' + s_agent + "_GR50_and_serum_concentration.pdf"
    fig.savefig(fig_name, bbox_inches='tight')
