#!/usr/bin/env python

#program: drug_response_analysis_pipeline.py
#Author: Hongmei Zhang
#data:Fri Feb 24 13:47:05 2017 v1
#Edited: Janice Patterson
#v2 edits 2020-04-15

#description
#Glue code for pipeline

#import libraries
import sys
import os
import shutil
import re
import csv
import argparse
from drug_plate_median import CellCountMedian
from generate_input_for_gr50_tools import time0_median, generate_input
from drug_plate_mean import CellCountMeanFilter
from boxplot_for_raw_cv import boxplot_for_raw_cv
from generate_the_cv_table import generate_boxplot_table_for_cv
from add_max_drug_does_to_GRmetrics import add_max_conc_to_GRm

#arguments parse
parser = argparse.ArgumentParser(description='')
parser.add_argument('--d', help='input the directory of raw data')
parser.add_argument('--w', help='input the work directory')
parser.add_argument('--f', help='input the file names of raw data')
parser.add_argument('--s', help='input the sample information file')
parser.add_argument('--l', help='input the layout of palte')
parser.add_argument('--c', help='input the cell line class: JWGray_BCCL_classifications_v5-1.txt')
parser.add_argument('--sc', help='input the drug_serum_conc.txt')

args = parser.parse_args()

##load path of scripts
#from drug_plate_median import CellCountMedian
#from generate_input_for_gr50_tools import time0_median, generate_input

#mkdir directory
dir_var=os.path.join(args.w, "data_variance")
dir_layout_heatmap=os.path.join(args.w, "raw_data_heatmap")
dir_CellCount=os.path.join(args.w, "CellCount")
dir_GR=os.path.join(args.w, "GrowthRate")
dir_cell_line_plot=os.path.join(args.w, "GrowthRate_plot/GR_plot_by_cell_line")
dir_drug_plot=os.path.join(args.w, "GrowthRate_plot/GR_plot_by_drug")
dir_stat_plot=os.path.join(args.w, "GR_stat_plot")
dir_serum_plot=os.path.join(args.w, "GR50_and_serum_concentration_plot")
dir_GRmax_plot=os.path.join(args.w, "GRmax_plot")
os.system("mkdir -p %s %s %s %s %s %s %s %s %s" % (dir_var, dir_layout_heatmap, dir_CellCount, dir_GR, dir_cell_line_plot, dir_drug_plot, dir_stat_plot, dir_serum_plot, dir_GRmax_plot))

#read in the files one by one
with open (args.f, "r") as tsvfile:
    h_file=csv.reader(tsvfile, delimiter="\t")
    for lst_row in h_file:
        s_file_name=lst_row[0]
        s_file_path=os.path.join(args.d, s_file_name)
        #calculate the median value of cell count for each dose of drug
        df_median=CellCountMedian(s_file_path, args.l)
        s_plateID=re.search("(\w+).txt", s_file_name)
        s_out_median=os.path.join(dir_CellCount, s_plateID.group(1) + "_CellCount_median.txt")
        df_median.to_csv(s_out_median, sep="\t")
        #plot the heatmap for the raw cell count
        heatmap_file=os.path.join(dir_layout_heatmap, s_plateID.group(1) + "_heatmap.pdf")
        os.system('python3 /Users/patterja/Workspace/hongmei/drug_response/heatmap_for_384_matrix.py --i %s --l %s --h %s' % (s_file_path, args.l, heatmap_file))
        #remove the outlier of triplates and calculate the coefficient of variance
        (d_mean, d_cv,i_remove_num, plateID)=CellCountMeanFilter(s_file_path, args.l)
        #write down the coefficient of variance
        s_out_var=os.path.join(dir_var, s_plateID.group(1) + "_CellCount_variance.txt")
        with open (s_out_var, "w") as tsvcv:
            h_cv=csv.writer(tsvcv, delimiter="\t")
            h_cv.writerow(['plateID', "Compound", "Concentration", 'Unit', 'mean_raw', 'stdev_raw', "cv_raw", 'mean_filter', 'stdev_filter', "cv_filter"])
            for key in d_cv.keys():
                lst_cv=d_cv[key].split("\t")
                h_cv.writerow(lst_cv)

#########boxplot for the coefficient of variation########
os.chdir(dir_var)
#cat the cv of raw cell count into one file
os.system('cut -f1,7 *.txt >tmp')
os.system('sort tmp>tmp_raw')
#prepare the cv table for box plotting
d_raw_cv=generate_boxplot_table_for_cv('./tmp_raw')
os.system('rm tmp*')
#write out the table
with open ('./cv_for_raw_cell_count.tsv', 'w') as tsv_out:
    h_out=csv.writer(tsv_out, delimiter="\t")
    for key in sorted(d_raw_cv.keys()):
        h_out.writerow(d_raw_cv[key])
#boxplot with the cv
fig_raw=boxplot_for_raw_cv(cv_file='./cv_for_raw_cell_count.tsv', sam_file=args.s)
raw_cv_fig='./boxplot_for_cv_of_raw_CellCount.pdf'
fig_raw.savefig(raw_cv_fig, bbox='tight')

#generate the input of gr50_tools
#######read in the sample information#########
#Plate Barcode   Image Barcode   Exp Time        Time Unit       Image Barcode of Time0  Well Number of Time0    Cell Line       Cell Number     Unit of Cell Number
#SM384000011     SM300011        72      h       SM300010        B2;B3;B4;B5;B6;B7;B8;B9;B10;B11;B12;B13;B14     HCC1569 4000    cell/well
with open (args.s, "r") as tsvsamp:
    h_samp=csv.reader(tsvsamp, delimiter="\t")
    header=next(h_samp)
    i_num=0
    d_output={}
    for lst_row in h_samp:
        (drug_plateID, time, time0_plateID, wells, cellLine)=(lst_row[1], lst_row[2], lst_row[4],lst_row[5], lst_row[6])
        t0_median=time0_median(plateID=time0_plateID, wells=wells, path=args.d)
        drug_file=os.path.join(dir_CellCount, drug_plateID + "_CellCount_median.txt")
        if (os.path.exists(drug_file)):
            (i_num, d_output)=generate_input(fname=drug_file, testTime=time, t0_median=t0_median, CL=cellLine, num=i_num, d_out=d_output)
        else:
            print ("{} not exists".format(drug_file))
########write down the result#########  
#cell_line  agent   perturbation    replicate   time    concentration   cell_count  cell_count__ctrl    cell_count__time0
#MCF10A drugA   0   1   48  0.001   1131    1212.5  299.5
s_gr50_input=os.path.join(dir_CellCount, "input_file_of_gr50_tools.txt")
with open (s_gr50_input, "w") as tsvout:
    h_out=csv.writer(tsvout, delimiter="\t")
    lst_header=["cell_line", "agent", "timepoint", "concentration", "cell_count", "cell_count__ctrl", "cell_count__time0"]
    h_out.writerow(lst_header)
    for key in d_output.keys():
        h_out.writerow(d_output[key])

###########calculate the GR values and metrics###########
print("running gr50 tools")
s_GRvalue=os.path.join(dir_GR, "GRvalues.txt")
s_GRmetrics=os.path.join(dir_GR, "GRmetrics.txt")

# Run below code separately after installing gr50_tools from hafner and using his tool which is built with python2
os.system("/usr/bin/python /Users/patterja/Workspace/gr50_tools_repos/gr_metrics/SRC/python/scripts/add_gr_column.py %s > %s" % (s_gr50_input, s_GRvalue))
os.system("/usr/bin/python /Users/patterja/Workspace/gr50_tools_repos/gr_metrics/SRC/python/scripts/compute_gr_metrics.py %s > %s" % (s_GRvalue, s_GRmetrics))
os.system("/usr/bin/python /Users/patterja/Workspace/hongmei/drug_response/plot_curves.py --metrics %s --values %s --outdir_bydrug %s --outdir_bycellline %s" % (s_GRmetrics, s_GRvalue, dir_drug_plot, dir_cell_line_plot))
#os.system("python3 /Users/patterja/Workspace/hongmei/scripts/plot_figure.py --v %s --m %s --d %s" % (s_GRvalue, s_GRmetrics, dir_drug_plot))

############add max does of each drug into the GRmetrics.txt############
# make copy of original GRmetrics file to prevent confusion
shutil.copy(s_GRmetrics, os.path.join(dir_GR, "GRmetrics_original.txt"))
df_GRm_edit=add_max_conc_to_GRm(f_GRm=s_GRmetrics, f_layout=args.l)
df_GRm_edit.to_csv(s_GRmetrics, sep="\t")

###########GR matrics statistics plot####################

os.system('python /Users/patterja/Workspace/hongmei/drug_response/boxplot_for_GRmetrics.py --i %s --w %s' % (s_GRmetrics, dir_stat_plot))
os.system('python /Users/patterja/Workspace/hongmei/drug_response/gr50_and_curve_count_barplot.py --i %s --p %s --w %s' % (s_GRmetrics, "0.05", dir_stat_plot))
print("\nFinished GR matrics statistics plots\n")
###########GR50 and drug serum concentration plot#############

print("\nStarting GR50 and drug serum concentration plot\n")
os.system('python /Users/patterja/Workspace/hongmei/drug_response/barplot_for_GR50_and_serum_conc.py --i %s --c %s --s %s --w %s' % (s_GRmetrics, args.c, args.sc, dir_serum_plot))
print('python /Users/patterja/Workspace/hongmei/drug_response/barplot_for_GR50_and_serum_conc.py --i %s --c %s --s %s --w %s' % (s_GRmetrics, args.c, args.sc, dir_serum_plot))


###########GRmax plot#############

print("\nStarting GRmax plot\n")

#calculate the GR values of drugs at serum concentration, and add the values into GRmetrics
s_GRmetrics_new=os.path.join(dir_GR, "GRmetrics_GRserum_added.txt")
os.system('python /Users/patterja/Workspace/hongmei/drug_response/calculate_GR_at_serum_conc.py --m %s --gr %s --s %s --o %s' % (s_GRmetrics, s_GRvalue, args.sc, s_GRmetrics_new))
os.system('python /Users/patterja/Workspace/hongmei/drug_response/barplot_for_GRmax_and_GRserum.py --i %s --c %s  --w %s' % (s_GRmetrics_new, args.c, dir_GRmax_plot))
print(('python /Users/patterja/Workspace/hongmei/drug_response/barplot_for_GRmax_and_GRserum.py --i %s --c %s  --w %s' % (s_GRmetrics_new, args.c, dir_GRmax_plot)))

###########GR50 heatmap plot############

print("\nGR50 heatmap plot\n")

os.system('python /Users/patterja/Workspace/hongmei/drug_response/heatmap_for_GR50_median_center.py --i %s --c %s --w %s' % (s_GRmetrics, args.c, dir_stat_plot))
os.system(
    'Rscript /Users/zhanho/Project/SMMART/scripts/heatmap3_new.R --data %s --t %s --r_side_label %s --c_side_label %s --w %d --h %d --min %f --max %f --out %s' % (
    gr50_file, "'Log10(GR50) nM'", drug_class_file, cc_class_file, 20, 10, -2, 2, heatmap_file))

print('python /Users/patterja/Workspace/hongmei/drug_response/heatmap_for_GR50_median_center.py --i %s --c %s --w %s' % (s_GRmetrics, args.c, dir_stat_plot))
print("The End")

