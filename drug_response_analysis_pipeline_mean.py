#program: drug_response_analysis_pipeline.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Feb 24 13:47:05 2017

#description
#this program will combine all the scripts together

#import libraries
import sys
import os
import argparse
import pandas as pd
import re
import csv
import matplotlib.pyplot as plt
from drug_response.drug_plate_mean import CellCountMeanFilter
from drug_response.generate_input_for_gr50_tools import time0_mean_filter, generate_input
from drug_response.scatter_plot_for_stdev import scatter_plot
from drug_response.generate_the_cv_table import generate_boxplot_table_for_cv
from drug_response.boxplot_for_raw_cv import boxplot_for_raw_cv
from drug_response.boxplot_for_filtered_cv import boxplot_for_filtered_cv

#arguments parse
parser = argparse.ArgumentParser(description='')
parser.add_argument('--d', help='input the data directory')
parser.add_argument('--w', help='input the work directory')
parser.add_argument('--f', help='input the file names of raw data')
parser.add_argument('--s', help='input the sample information file')
parser.add_argument('--l', help='input the layout of palte')

args = parser.parse_args()

##load path of scripts
#sys.path.append("/Users/zhanho/Project/SMMART/scripts/")
#from drug_plate_median import CellCountMedian
#from generate_input_for_gr50_tools import time0_median, generate_input

#mkdir directory
dir_var=os.path.join(args.w, "data_variance")
dir_CellCount=os.path.join(args.w, "CellCount")
dir_GR=os.path.join(args.w, "GrowthRate")
dir_cell_line_plot=os.path.join(args.w, "GrowthRate_plot/cell_line_plot")
dir_drug_plot=os.path.join(args.w, "GrowthRate_plot/drug_plot")
os.system("mkdir -p %s %s %s %s %s" % (dir_var, dir_CellCount, dir_GR, dir_drug_plot, dir_cell_line_plot))

#read in the file names, filter the data and calculate the mean of cell count
d_remove_num={}
with open (args.f, "r") as tsvfile:
    h_file=csv.reader(tsvfile, delimiter="\t")
    for lst_row in h_file:
        s_file_name=lst_row[0]
        s_file_path=os.path.join(args.d, s_file_name)
        #filter the data and calculate the mean value of cell count
        (d_mean, d_cv,i_remove_num, plateID)=CellCountMeanFilter(s_file_path, args.l)
        d_remove_num.update({plateID: i_remove_num})
        #write out the mean value
        s_plateID=re.search("(\w+).txt", s_file_name)
        s_out_mean=os.path.join(dir_CellCount, s_plateID.group(1) + "_CellCount_mean.txt")
        with open (s_out_mean, "w") as tsvout:
            h_out=csv.writer(tsvout, delimiter="\t")
            h_out.writerow(['plateID', "Compound", "Concentration", 'Unit', "Mean of CellCount"])
            for key in d_mean.keys():
                lst_mean=key.split("\t")
                lst_mean.append(d_mean[key])
                h_out.writerow(lst_mean)
        #write out the variance
        s_out_var=os.path.join(dir_var, s_plateID.group(1) + "_CellCount_variance.txt")
        with open (s_out_var, "w") as tsvcv:
            h_cv=csv.writer(tsvcv, delimiter="\t")
            h_cv.writerow(['plateID', "Compound", "Concentration", 'Unit', 'mean_raw', 'stdev_raw', "cv_raw", 'mean_filter', 'stdev_filter', "cv_filter"])
            for key in d_cv.keys():
                lst_cv=d_cv[key].split("\t")
                h_cv.writerow(lst_cv)
        ##scatter plot for the coefficient variance
        #s_fig_file=os.path.join(dir_var, s_plateID.group(1) + "_CellCount_variance.pdf")
        #fig=scatter_plot(s_out_var)
        #fig.savefig(s_fig_file, bbox='tight')

#barplot for the remove well number
df_remove_num=pd.DataFrame(list(d_remove_num.items()), columns=['PlateID', 'RemoveNum'])
df_sam=pd.read_csv(args.s, sep="\t")
df_sam=df_sam[['Image Barcode', 'Cell Line']]
df_sam.columns=['PlateID', 'CellLine']
df_sam_remove_num=pd.merge(df_sam, df_remove_num, on='PlateID')
df_sam_remove_num=df_sam_remove_num[['CellLine', 'RemoveNum']]
df_sam_remove_num=df_sam_remove_num.set_index('CellLine')
f_data_file=os.path.join(dir_var, 'remove_number.txt')
df_sam_remove_num=df_sam_remove_num.sort_values(by='RemoveNum', ascending=False)
fig_remove_num=plt.figure(figsize=(15,7))
ax=fig_remove_num.add_subplot(1,1,1)
#fig_remove_num.subplots_adjust(bottom=6)
fig1=df_sam_remove_num.plot(kind='bar', color='orange', ax=ax, legend=False, fontsize=7)
fig1.set_ylabel('Deleted Well Number')
f_fig_name=os.path.join(dir_var, 'remove_number.pdf')
fig_remove_num.savefig(f_fig_name, bbox='tight')

#########boxplot for the coefficient of variation########
os.chdir(dir_var)
#cat the cv of raw cell count into one file
os.system('cut -f1,7 *.txt >tmp')
os.system('sort tmp>tmp_raw')
#prepare the cv table for box plotting
d_raw_cv=generate_boxplot_table_for_cv('./tmp_raw')
#write out the table
with open ('./tmp_raw_cv.tsv', 'w') as tsv_out:
    h_out=csv.writer(tsv_out, delimiter="\t")
    for key in sorted(d_raw_cv.keys()):
        #lst_cv=d_cv[key].split("\t")
        #print (type(d_cv[key]))
        h_out.writerow(d_raw_cv[key])
#boxplot with the cv
fig_raw=boxplot_for_raw_cv(cv_file='./tmp_raw_cv.tsv', sam_file=args.s)
raw_cv_fig='./cv_of_raw_CellCount_boxplot.pdf'
fig_raw.savefig(raw_cv_fig, bbox='tight')

#cat the cv of filtered cell count into one file
os.system('cut -f1,10 *.txt >tmp')
os.system('sort tmp>tmp_filter')
#prepare the cv table for box plotting
d_filter_cv=generate_boxplot_table_for_cv('./tmp_filter')
#write out the table
with open ('./tmp_filter_cv.tsv', 'w') as tsv_out:
    h_out=csv.writer(tsv_out, delimiter="\t")
    for key in sorted(d_filter_cv.keys()):
        h_out.writerow(d_filter_cv[key])
#boxplot with the cv
fig_filter=boxplot_for_filtered_cv(filter_cv_file='./tmp_filter_cv.tsv', raw_cv_file='./tmp_raw_cv.tsv', sam_file=args.s)
filter_cv_fig='./cv_of_filtered_CellCount_boxplot.pdf'
fig_filter.savefig(filter_cv_fig, bbox='tight')


#generate the input of gr50_tools
#######read in the sample#########
#Plate Barcode   Image Barcode   Exp Time        Time Unit       Image Barcode of Time0  Well Number of Time0    Cell Line       Cell Number     Unit of Cell Number
#SM384000011     SM300011        72      h       SM300010        B2;B3;B4;B5;B6;B7;B8;B9;B10;B11;B12;B13;B14     HCC1569 4000    cell/well
with open (args.s, "r") as tsvsamp:
    h_samp=csv.reader(tsvsamp, delimiter="\t")
    header=next(h_samp)
    i_num=0
    d_output={}
    for lst_row in h_samp:
        (drug_plateID, time, time0_plateID, wells, cellLine)=(lst_row[1], lst_row[2], lst_row[4],lst_row[5], lst_row[6])
        t0_mean=time0_mean_filter(plateID=time0_plateID, wells=wells, path=args.d)
        drug_file=os.path.join(dir_CellCount, drug_plateID + "_CellCount_mean.txt")
        if (os.path.exists(drug_file)):
            (i_num, d_output)=generate_input(fname=drug_file, testTime=time, t0_median=t0_mean, CL=cellLine, num=i_num, d_out=d_output)
        else:
            print ("{} not exists".format(drug_file))
########write out the result#########  
#cell_line  agent   perturbation    replicate   time    concentration   cell_count  cell_count__ctrl    cell_count__time0
#MCF10A drugA   0   1   48  0.001   1131    1212.5  299.5
s_gr50_input=os.path.join(dir_CellCount, "input_file_of_gr50_tools.txt")
with open (s_gr50_input, "w") as tsvout:
    h_out=csv.writer(tsvout, delimiter="\t")
    lst_header=["cell_line", "agent", "time", "concentration", "cell_count", "cell_count__ctrl", "cell_count__time0"]
    h_out.writerow(lst_header)
    for key in d_output.keys():
        h_out.writerow(d_output[key])

###########calculate the GR values and metrics###########
s_GRvalue=os.path.join(dir_GR, "GRvalues.txt")
s_GRmetrics=os.path.join(dir_GR, "GRmetrics.txt")
os.system("python3 /Users/zhanho/src/github/gr50_tools/SRC/python/scripts/add_gr_column.py %s > %s" % (s_gr50_input, s_GRvalue))
os.system("python3 /Users/zhanho/src/github/gr50_tools/SRC/python/scripts/compute_gr_metrics.py %s > %s" % (s_GRvalue, s_GRmetrics))
os.system("python3 /Users/zhanho/src/github/gr50_tools/SRC/python/scripts/plot_figure.py --v %s --m %s --d %s --d2 %s" % (s_GRvalue, s_GRmetrics, dir_drug_plot, dir_cell_line_plot))
