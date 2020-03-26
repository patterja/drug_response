#program: boxplot_for_GRmax_and_serum_conc.py
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
import math
import seaborn as sns

parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', help='input the GRmetrics_GRserum_added.txt')
parser.add_argument('--c', help='input the cell line class and color')
parser.add_argument('--w', help='input the work directory ')

args = parser.parse_args()

#read in the GRmetrics.txt
df_data=pd.read_csv(args.i, sep="\t")
df_data=df_data.set_index('cell_line')
#read in cell line class and color
df_color=pd.read_csv(args.c, sep="\t", index_col=0)
#combine these to matrix
df_merge=pd.merge(df_data, df_color, left_index=True, right_index=True)
#extract the sub-matrix
df_sub=df_merge[['agent', 'GRmax', 'GRserum', 'Class','Color']]
df_sub.index.name='cell_line'
#df_sub=df_sub.reset_index()
#df_sub=df_sub.set_index('agent')
df_GRserum=df_sub[["GRserum", "Class", "agent"]]
df_GRserum=df_GRserum.sort_values(by=['agent', 'Class'])
df_GRmax=df_sub[["GRmax", "Class", "agent"]]
df_GRmax=df_GRmax.sort_values(by=['agent','Class'])
#lst_agents=df_GRmax.index.unique()

###plot the figure for GRmax
sns.set(style="ticks")
sns.set(font_scale=1.3)
grid = sns.FacetGrid(df_GRmax, col="agent", col_wrap=4, size=4)
grid.map(sns.boxplot, "Class", "GRmax", color='purple')
grid.map(sns.swarmplot, "Class", "GRmax", color='black')
plt.subplots_adjust(bottom=0.1, top=0.93)
grid.fig.suptitle("GRmax")
grid.set_xticklabels(rotation=90)
grid.set(ylim=(-1, 1.5))
figfile=args.w + "/GRmax_boxplot.pdf"
plt.savefig(figfile, format='pdf')

###plot the figure for GRserum
sns.set(style="ticks")
sns.set(font_scale=1.3)
grid = sns.FacetGrid(df_GRserum, col="agent", col_wrap=4, size=4)
grid.map(sns.boxplot, "Class", "GRserum", color="darkgreen")
grid.map(sns.swarmplot, "Class", "GRserum", color='black')
plt.subplots_adjust(bottom=0.1, top=0.93)
grid.fig.suptitle("GRserum")
grid.set_xticklabels(rotation=90)
grid.set(ylim=(-1, 1.5))
figfile=args.w + "/GRserum_boxplot.pdf"
plt.savefig(figfile, format='pdf')
