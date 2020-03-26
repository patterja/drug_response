#program: drug_plate_mean.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Feb 24 13:46:39 2017

#description
#calculate the mean of Cell count

#import libraries
import argparse
import pandas as pd
import re
import os
import math
import statistics
import csv

def CellCountMeanFilter(datafile, layout):
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
    #get the mean
    df_merge_mean=df_merge.groupby(['plateID', "Compound", "Concentration", 'Unit']).mean()
    #set index for df_merge
    df_merge=df_merge.set_index(['plateID', 'Compound', 'Concentration', 'Unit'])
    #filter the outlie data
    d_mean_filter={}
    d_cv={}
    i_remove_num=0
    for index in df_merge.index.unique():
        if index in df_merge_mean.index:
            mean_raw=df_merge_mean.ix[index]["CellCount"]
            lst_cc=list(df_merge.ix[index]["CellCount"])
            index=list(index)
            index[2]=str(index[2])
            s_index="\t".join(index)
            #calculate the stand deviation of raw list
            if (len(lst_cc)>1):
                f_stdev_raw=statistics.stdev(lst_cc)
                f_cv_raw=f_stdev_raw/mean_raw
            else:
                f_stdev_raw=0
                f_cv_raw=0
            #remove the value if the fold change between value and mean_value greater or smaller then 10 
            for value in lst_cc:
                if (value >0 and (abs(math.log10(value) - math.log10(mean_raw))>1)):
                    i_remove_num+=1
                    lst_cc.remove(value)
            mean_filter=statistics.mean(lst_cc)
            try:
                d_mean_filter[s_index]=mean_filter
            except KeyError:
                d_mean_filter.update({s_index: mean_filter})
            #calculate the stand deviation of filter list
            if (len(lst_cc)>1):
                f_stdev_filter=statistics.stdev(lst_cc)
                f_cv_filter=f_stdev_filter/mean_filter
            else:
                f_stdev_filter=0
                f_cv_filter=0
            #print the stand deviation out
            mean_raw="{0:.2f}".format(mean_raw)
            mean_filter="{0:.2f}".format(mean_filter)
            f_stdev_raw="{0:.2f}".format(f_stdev_raw)
            f_stdev_filter="{0:.2f}".format(f_stdev_filter)
            f_cv_raw="{0:.2f}".format(f_cv_raw)
            f_cv_filter="{0:.2f}".format(f_cv_filter)
            s_out=s_index + "\t" + str(mean_raw) + "\t" + str(f_stdev_raw) + "\t" + str(f_cv_raw) + "\t" + str(mean_filter) + "\t" + str(f_stdev_filter) + "\t" + str(f_cv_filter)
            try:
                d_cv[s_index]=s_out
            except KeyError:
                d_cv.update({s_index: s_out})
    
    return (d_mean_filter, d_cv, i_remove_num, plateID.group(1))

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the raw data of cell count')
    parser.add_argument('--l', help='input the layout file of plate')
    parser.add_argument('--m', help='output the mean value of cell count')
    parser.add_argument('--cv', help='output the standand deviation of cell count')
    
    args = parser.parse_args()
    #calculate the mean of cell count
    (d_mean, d_cv, i_remove_num, plateID)=CellCountMeanFilter(args.i, args.l)
    #write down the mean value of cell count
    with open (args.m, "w") as tsvout:
        h_out=csv.writer(tsvout, delimiter="\t")
        h_out.writerow(['plateID', "Compound", "Concentration", 'Unit', "Mean of CellCount"])
        for key in d_mean.keys():
            lst_mean=key.split("\t")
            lst_mean.append(d_mean[key])
            h_out.writerow(lst_mean)

    #write down the standand deviation of cell count
    with open (args.cv, "w") as tsvcv:
        h_cv=csv.writer(tsvcv, delimiter="\t")
        h_cv.writerow(['plateID', "Compound", "Concentration", 'Unit', 'mean_raw', 'stdev_raw', "cv_raw", 'mean_filter', 'stdev_filter', "cv_filter"])
        for key in d_cv.keys():
            lst_cv=d_cv[key].split("\t")
            h_cv.writerow(lst_cv)
