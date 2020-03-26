#program: boxplot_for_cv.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Mar 10 10:56:58 2017

#description

#import libraries
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def boxplot_for_raw_cv(cv_file, sam_file):
    #read in the sample information file
    df_sam=pd.read_csv(sam_file, sep="\t", index_col=None, header=0)
    df_sam=df_sam[['Image Barcode', 'Cell Line']]
    df_sam=df_sam.set_index('Image Barcode')
    #read in the data of filtered value
    df_data=pd.read_csv(cv_file, sep="\t")
    df_data=df_data.drop('plateID', 1)
    #transpose the data
    df_data=df_data.transpose()
    #merge the cell line name with the data
    df_merge=pd.merge(df_sam, df_data, left_index=True, right_index=True)
    df_merge=df_merge.set_index('Cell Line')
    df_merge=df_merge.transpose()
    #get the median value of df_data
    df_merge_median=df_merge.median()
    df_merge_median=df_merge_median.sort_values()
    #sort df_merge with cv median value of raw value
    df_merge=df_merge[df_merge_median.index] 
	#plot the figure
    fig=plt.figure(figsize=(18,9))
    ax=fig.add_subplot(1,1,1)
    color = dict(boxes='YellowGreen', whiskers='DarkOrange', medians='DarkBlue', caps='Gray')
    fig1=df_merge.plot(kind='box', rot=90, color=color, sym='r+',  ax=ax, fontsize=7)
    length=len(df_data.columns)
    fig2=plt.plot([-1, length], [1.0, 1.0], linestyle='--', color="black")
    fig1.set_ylabel('coefficient of variation')
    #return figure
    return fig

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the coefficient data')
    parser.add_argument('--s', help='input the sample information table:SM384_platesinfo.txt')
    parser.add_argument('--o', help='output the boxplot figure')
    
    args = parser.parse_args()
    #generate the boxplot
    fig=boxplot_for_raw_cv(cv_file=args.i, sam_file=args.s) 
    fig.savefig(args.o, bbox='tight')
