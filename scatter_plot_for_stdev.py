#program: scatter_plot_for_stdev.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Thu Mar  9 10:05:04 2017

#description

#import libraries
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import math

def scatter_plot(in_file):
    #read in the data
    df_data=pd.read_csv(in_file, sep="\t")
    df_stdev=df_data[['cv_raw', 'cv_filter']]
    max_value=math.ceil((df_stdev.max()).max())
    min_value=math.floor((df_stdev.min()).min())
    #scatter plot
    fig=plt.figure(figsize=(5,5))
    ax=fig.add_subplot(1,1,1)
    fig1=df_stdev.plot(kind="scatter", x='cv_raw', y='cv_filter', color='blue', ax=ax, xlim=(min_value, max_value), ylim=(min_value,max_value))
    fig2=plt.plot([min_value, max_value], [min_value, max_value], color="black", linestyle="--")
    fig1.set_xlabel("cv of raw value")
    fig1.set_ylabel("cv of filtered value")
    return fig

    
if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the stdev data file')
    parser.add_argument('--o', help='output the plotting')
    
    args = parser.parse_args()
    #plot the figure
    fig=scatter_plot(args.i)
    fig.savefig(args.o, bbox='tight')
