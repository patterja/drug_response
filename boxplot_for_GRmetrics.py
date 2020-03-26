#program: boxplot_for_GRmetrics.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Apr  7 17:33:33 2017

#description

#import libraries
import argparse
import pandas as pd
import numpy as np
import os.path as path
import matplotlib.pyplot as plt
import seaborn as sns

def sns_boxplot(df_data, y_label, log10):
    #drop inf and -inf
    df_data=df_data.replace([np.inf, -np.inf], np.nan).dropna()
    if(log10):
        df_data['value']=df_data['value'].map(lambda x: np.log10(x))
    #sort the index by median
    lst_sort_index=df_data.groupby('agent').median().sort_values('value').index
    df_data=df_data.set_index('agent')
    df_data=df_data.ix[lst_sort_index]
    #set the parameters for figure
    params = {'figure.figsize': (15, 10),
         'axes.labelsize': '12',
         'axes.titlesize':'15',
         'xtick.labelsize':'10',
         'ytick.labelsize':'10'}
    plt.rcParams.update(params)
    #set the figure size
    fig=plt.figure()
    ax=fig.add_subplot(1,1,1)
    fig1=sns.boxplot(x=df_data.index, y="value", data=df_data, ax=ax)
    #rotate the xtick lables
    for item in fig1.get_xticklabels():
        item.set_rotation(90)
    #Add some margin to the bottom so the labels aren't cut-off
    plt.subplots_adjust(bottom=0.25)
    #Use swarmplot() to show the datapoints on top of the boxes:
    fig2=sns.swarmplot(x=df_data.index, y="value", data=df_data, color=".25")
    #set x and y labels
    fig1.set_xlabel('Compound')
    fig1.set_ylabel(y_label)
    #return figure
    return fig

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the GRmetrics.txt')
    parser.add_argument('--w', help='input the work directory')
    
    args = parser.parse_args()

    #read in the data
    df_data=pd.read_csv(args.i, sep="\t")
    #get the submatrix of GR50
    df_gr50=df_data[["agent", "GR50"]]
    df_gr50.columns=["agent", "value"]
    fig_gr50=sns_boxplot(df_data=df_gr50, y_label='Log10(GR50) nM', log10=True)
    fig_name=path.join(args.w, 'GR50_boxplot.pdf')
    fig_gr50.savefig(fig_name, bbox_inches='tight')
    #boxplot for the GRmax
    df_grmax=df_data[["agent", "GRmax"]]
    df_grmax.columns=["agent", "value"]
    fig_grmax=sns_boxplot(df_data=df_grmax, y_label='GRmax', log10=False)
    fig_name=path.join(args.w, "GRmax_boxplot.pdf")
    fig_grmax.savefig(fig_name, bbox_inches='tight')
    #boxplot for the GR_AOC
    df_graoc=df_data[["agent", "GR_AOC"]]
    df_graoc.columns=["agent", "value"]
    df_graoc=df_graoc.set_index('agent')
    df_graoc[df_graoc<0]=0
    df_graoc=df_graoc.reset_index()
    fig_graoc=sns_boxplot(df_data=df_graoc, y_label='GR AOC', log10=False)
    fig_name=path.join(args.w, "GRAOC_boxplot.pdf")
    fig_graoc.savefig(fig_name, bbox_inches='tight')
    #boxplot for the hGR
    df_grhGR=df_data[["agent", "h_GR"]]
    df_grhGR.columns=["agent", "value"]
    fig_grhGR=sns_boxplot(df_data=df_grhGR, y_label='Hill coeff.(hGR)', log10=False)
    fig_name=path.join(args.w, "hGR_boxplot.pdf")
    fig_grhGR.savefig(fig_name, bbox_inches='tight')
