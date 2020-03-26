#program: boxplot_for_cv.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Mar 10 10:56:58 2017

#description

#import libraries
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def boxplot_for_filtered_cv(filter_cv_file, raw_cv_file, sam_file):
    #read in the sample information file
    df_sam=pd.read_csv(sam_file, sep="\t")
    df_sam=df_sam[['Image Barcode', 'Cell Line']]
    df_sam=df_sam.set_index('Image Barcode')
    
    #read in the data of filtered value
    df_data=pd.read_csv(filter_cv_file, sep="\t")
    df_data=df_data.drop('plateID', 1)
    #transpose the data
    df_data=df_data.transpose()
    #merge the cell line name with the data
    df_merge=pd.merge(df_sam, df_data, left_index=True, right_index=True)
    df_merge=df_merge.set_index('Cell Line')
    df_merge=df_merge.transpose()
    
    ###############sort by the cv median of raw value################
    #read in the cv data of raw value
    df_data_raw=pd.read_csv(raw_cv_file, sep="\t")
    df_data_raw=df_data_raw.drop('plateID', 1)
    #transpose the data
    df_data_raw=df_data_raw.transpose()
    #merge the cell line name with the data
    df_merge_raw=pd.merge(df_sam, df_data_raw, left_index=True, right_index=True)
    df_merge_raw=df_merge_raw.set_index('Cell Line')
    df_merge_raw=df_merge_raw.transpose()
    #get the median value of df_data
    df_merge_raw_median=df_merge_raw.median()
    df_merge_raw_median=df_merge_raw_median.sort_values()
    
    #sort df_data with cv median value of raw value
    df_merge=df_merge[df_merge_raw_median.index] 
    
    #plot the figure
    fig=plt.figure(figsize=(15,7))
    ax=fig.add_subplot(1,1,1)
    color = dict(boxes='DarkGreen', whiskers='DarkOrange', medians='DarkBlue', caps='Gray')
    fig1=df_merge.plot(kind='box', rot=90, color=color, sym='r+',  ax=ax, fontsize=7)
    length=len(df_data.columns)
    fig2=plt.plot([-1, length], [1.0, 1.0], linestyle='--', color="black")
    fig1.set_ylabel('coefficient of variation')
    #return fig
    return fig

if __name__=='__main__':    
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the filtered coefficient data')
    parser.add_argument('--i2', help='input the raw coefficient data')
    parser.add_argument('--s', help='input the sample info file: SM384_platesinfo.txt')
    parser.add_argument('--o', help='output the boxplot figure')
    
    args = parser.parse_args()
    #generate fig
    fig=boxplot_for_filtered_cv(filter_cv_file=args.i, raw_cv_file=args.i2, sam_file=args.s) 
    fig.savefig(args.o, bbox='tight')
