#program: gr50_and_curve_count_barplot.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Apr  7 11:41:34 2017

#description

#import libraries
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def curve_and_gr50_stat(df_data, pvalue):
    #plot by cell line
    df_fit_cell_line=df_data[['cell_line', 'pval_GR']]
    df_fit_cell_line =df_fit_cell_line[df_fit_cell_line['pval_GR']<=pvalue]
    df_fit_count_cell_line=df_fit_cell_line.groupby("cell_line").count()
    df_fit_count_cell_line.columns=["Sigmoidal Fit"]
    df_gr50_cell_line=df_data[['cell_line', 'GR50']]
    df_gr50_cell_line =df_gr50_cell_line.replace([np.inf, -np.inf], np.nan).dropna()
    df_gr50_count_cell_line=df_gr50_cell_line.groupby("cell_line").count()
    df_gr50_count_cell_line.columns=["GR50"]
    df_cell_line=pd.merge(df_fit_count_cell_line, df_gr50_count_cell_line, left_index=True, right_index=True, how='outer')
    df_cell_line.index.name=None
    df_cell_line=df_cell_line.sort_values(by='Sigmoidal Fit', ascending=False)
    cell_line_fig=plt.figure(figsize=(18,8))
    cell_line_ax=cell_line_fig.add_subplot(1,1,1)
    df_cell_line.plot(kind='bar', color=["red","green"], ax=cell_line_ax)
    cell_line_ax.set_xlabel("Cell Lines")
    cell_line_ax.set_ylabel("Number of Compounds")
    #plot by drug
    df_fit_agent=df_data[['agent', 'pval_GR']]
    df_fit_agent =df_fit_agent[df_fit_agent['pval_GR']<=pvalue]
    df_fit_agent_count=df_fit_agent.groupby("agent").count()
    df_fit_agent_count.columns=["Sigmoidal Fit"]
    df_gr50_agent=df_data[['agent', 'GR50']]
    df_gr50_agent =df_gr50_agent.replace([np.inf, -np.inf], np.nan).dropna()
    df_gr50_agent_count=df_gr50_agent.groupby("agent").count()
    df_gr50_agent_count.columns=["GR50"]
    df_agent=pd.merge(df_fit_agent_count, df_gr50_agent_count, left_index=True, right_index=True, how='outer')
    df_agent.index.name=None
    df_agent=df_agent.sort_values(by='Sigmoidal Fit', ascending=False)
    agent_fig=plt.figure(figsize=(16,8))
    agent_ax=agent_fig.add_subplot(1,1,1)
    df_agent.plot(kind='bar', color=["red", "green"], ax=agent_ax)
    agent_ax.set_xlabel("Compounds")
    agent_ax.set_ylabel("Number of Cell Lines")
    return (agent_fig, cell_line_fig)

     
if __name__=='__main__':

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the GRmetrics.txt')
    parser.add_argument('--p', help='input the pvalue cutoff')
    parser.add_argument('--w', help='input the work directory')
    
    args = parser.parse_args()

    #read in the data
    df_gr=pd.read_csv(args.i, sep="\t")

    #barplot for the curve and gr50 number
    (fig_agent, fig_cell_line) = curve_and_gr50_stat(df_gr, float(args.p))
    fig_name_agent=args.w + "/" + "stat_by_agent.pdf"
    fig_name_cell_line=args.w + "/" + "stat_by_cell_line.pdf"
    fig_agent.savefig(fig_name_agent, bbox_inches='tight')
    fig_cell_line.savefig(fig_name_cell_line, bbox_inches='tight')
