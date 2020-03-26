#program: scanR2tsv.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Thu Feb 16 09:53:31 2017

#description
#calculate the cell count from the ScanR data and output as a tsv file

#import libraries
import argparse
import pandas as pd
import re
import os

def scanR2tsv(infile):
    #get the plate ID
    filename=os.path.basename(infile)
    plateID=re.search("(\w+)_Main", filename)
    #read in the file
    df_data=pd.read_csv(infile, sep="\t", index_col=0)
    #get the well columns and calculate the count for each well
    df_well=df_data["Well "]
    df_well_count=df_well.value_counts()
    df_well_count=pd.DataFrame(df_well_count)
    #add PlateID for each row
    df_well_count["PlateID"]=plateID.group(1)
    #re-organize the columns
    df_well_count=df_well_count.reset_index()
    df_well_count=df_well_count[["PlateID", "index", "Well "]]
    df_well_count.columns=["PlateID", "Well", "CellCount"]
    df_well_count=df_well_count.sort_values("Well")
    
    return (df_well_count, plateID.group(1))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the scanR data')
    parser.add_argument('--o', help='input the work directory')

    args = parser.parse_args()

    #run the function
    (df_cell_count, plateID)=scanR2tsv(args.i)
    outfile=args.o + "/" + plateID + "_CellCount.txt"
    df_cell_count.to_csv(outfile, sep="\t", index=False)


