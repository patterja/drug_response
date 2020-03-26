#program: heatmap_for_GR50.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Apr 14 11:22:39 2017

#description

#import libraries
import argparse
import pandas as pd
import numpy as np
import os
import csv

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the GRmetrics.txt')
    parser.add_argument('--c', help='input the cell line class file: JWGray_BCCL_classifications_v5-1.txt')
    parser.add_argument('--w', help='input the work directory')
    
    args = parser.parse_args()

    #read in the data
    df_data=pd.read_csv(args.i, sep="\t")
    ###generate GR50 matrix for heatmap
    df_gr50=df_data[["agent",'cell_line', "GR50", 'Maximal Does']]
    #replace inf to maximal does
    for i in range(len(df_gr50)):
        if(abs(df_gr50.ix[i,'GR50']) == np.inf):
            df_gr50.ix[i,'GR50']=df_gr50.ix[i,'Maximal Does'] 
    #convert the value into log10
    df_gr50['GR50']=np.log10(df_gr50['GR50'])
    df_tmp=df_gr50[['agent', 'GR50']]
    df_tmp_median=df_tmp.groupby('agent').median()
    df_tmp_median.columns=['median']
    df_gr50=df_gr50.set_index('agent')
    df_merge=pd.merge(df_gr50, df_tmp_median, left_index=True, right_index=True)
    df_merge['median_center']=df_merge['GR50'] - df_merge['median']
    df_merge.index.name='agent'
    df_merge=df_merge.reset_index()
    #extract the serum_center column
    df_serum_center=df_merge[["agent",'cell_line', "median_center"]]
    df_serum_center=df_serum_center.set_index(['cell_line', 'agent'])
    df_serum_center=df_serum_center.unstack()
    df_serum_center.index.name=None
    df_serum_center.columns=df_serum_center.columns.droplevel(None)
    df_serum_center.columns.name=None
    #get the cell line index of df_serum_center
    ds_cell_line=df_serum_center.index
    df_serum_center=df_serum_center.transpose()
    gr50_file=os.path.join(args.w, "gr50_median_center_log10.txt")
    df_serum_center.to_csv(gr50_file, sep="\t")
    
    ###read in the data of cell line class
    df_class=pd.read_csv(args.c, sep='\t', index_col=0)
    df_cell_line_class=pd.DataFrame(df_class.ix[ds_cell_line, 'Class'])
    df_cell_line_class=df_cell_line_class.transpose()
    cc_class_file=os.path.join(args.w, 'cell_line_class.txt')
    df_cell_line_class.to_csv(cc_class_file, sep='\t', index=None, header=False)

    ###drug class
    s_drug_class=' Compound' * len(df_gr50)
    lst_drug_class=s_drug_class.split(' ')
    df_drug_class=pd.DataFrame(lst_drug_class[1:])
    df_drug_class=df_drug_class.transpose()
    drug_class_file=os.path.join(args.w, 'drug_class.txt')
    df_drug_class.to_csv(drug_class_file, sep='\t', header=False)
    
    #plot the heatmap for GR50 matrix
    heatmap_file=os.path.join(args.w, "GR50_median_center_heatmap.pdf")
    os.system('Rscript /Users/zhanho/Project/SMMART/scripts/heatmap3_new.R --data %s --t %s --r_side_label %s --c_side_label %s --w %d --h %d --min %f --max %f --out %s' % (gr50_file, "'Log10(GR50) nM'", drug_class_file, cc_class_file, 20, 10, -2 , 2, heatmap_file))
