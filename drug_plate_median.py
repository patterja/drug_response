#program: drug_plate_median.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Feb 24 13:46:39 2017

#description
#calculate the median of Cell count

#import libraries
import argparse
import pandas as pd
import re
import os

def CellCountMedian(datafile, layout):
    #read in the data
    df_data=pd.read_csv(datafile, sep="\t")
    df_data=df_data[["Name / Description", "Objects"]]
    df_data.columns=["WellNumber", "CellCount"]
    basename=os.path.basename(datafile) 
    plateID=re.search('(\w+).txt', basename)
    df_data["plateID"]=plateID.group(1)
    #read in the layout
    df_layout=pd.read_csv(layout, sep="\t")
    #merge the data and layou with WellNumber
    df_merge=pd.merge(df_data, df_layout, on="WellNumber")
    #reset the columns of df_merge
    df_merge=df_merge[['plateID', 'Compound', 'Concentration', 'Unit', 'CellCount']]
    #get the median
    df_merge_median=df_merge.groupby(['plateID', "Compound", "Concentration", 'Unit']).median()
    
    return df_merge_median

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the raw data of cell count')
    parser.add_argument('--l', help='input the layout file of plate')
    parser.add_argument('--o', help='output the median value of cell count')
    
    args = parser.parse_args()
    #calculate the median of cell count
    df_median=CellCountMedian(args.i, args.l)
    #write down
    df_median.to_csv(args.o, sep="\t")


