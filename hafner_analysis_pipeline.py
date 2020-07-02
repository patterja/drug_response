#!/usr/bin/env python

# by: Janice Patterson

# Usage:
# python drug_response_analysis_pipeline.py --datdir /Users/patterja/Workspace/hongmei/data


#import libraries
import os
import argparse

#arguments parse
parser = argparse.ArgumentParser(description='')
parser.add_argument('--dataset', help='input metrics file from hafner')
parser.add_argument('--serum_conc', help='input the drug_serum_conc.txt')
parser.add_argument('--classif', help='JWGray_BCCL_classifications_v5-1_edit2020.txt')
parser.add_argument('--outdir', help='output file name')

args = parser.parse_args()

## load path of scripts
gr50_tools_scripts = "/Users/patterja/Workspace/gr50_tools_repos/gr_metrics/SRC/python/scripts"
drug_response_scripts ="/Users/patterja/Workspace/hongmei/drug_response"


## mkdir
dir_stat_plot=os.path.join(args.outdir, "GR_stat_plot")
dir_serum_plot=os.path.join(args.outdir, "GR50_and_serum_concentration_plot")
dir_GRmax_plot=os.path.join(args.outdir, "GRmax_plot")
os.system("mkdir -p %s %s %s %s" % (args.outdir, dir_stat_plot, dir_serum_plot, dir_GRmax_plot))


#from drug_plate_median import CellCountMedian
#from generate_input_for_gr50_tools import time0_median, generate_input

datasetout =os.path.join(outdir, "dataset_serum_added.txt")
os.system("/usr/bin/python %s/GRserum_calculation.py --m %s --s %s --o %s" % (drug_response_scripts, args.dataset, args.serum_conc, datasetout))

########### GR matrics statistics plot ####################

os.system('/usr/bin/python %s/boxplot_for_GRmetrics.py --i %s --w %s' % (drug_response_scripts, datasetout, dir_stat_plot))
os.system('/usr/bin/python %s/gr50_and_curve_count_barplot.py --i %s --p %s --w %s' % (drug_response_scripts, datasetout, "0.05", dir_stat_plot))
print("\nFinished GR matrics statistics plots\n")
########### GR50 and drug serum concentration plot #############

# This seems to require a serum conc file. While GRmax plots have the serum conc added to the GRmetric.txt.
# I think I should just change this file to act like barplot_for_GRmax_and_GRserum.py
print("\nStarting GR50 and drug serum concentration plot\n")
#os.system('/usr/bin/python %s/barplot_for_GR50_and_serum_conc.py --i %s --c %s --s %s --w %s' % (drug_response_scripts,datasetout, args.classification, args.serum_conc, dir_serum_plot))


###########GRmax plot#############

print("\nStarting GRmax plotting\n")

#calculate the GR values of drugs at serum concentration, and add the values into GRmetrics
os.system('/usr/bin/python %s/barplot_for_GRmax_and_GRserum.py --i %s --w %s' % (drug_response_scripts,datasetout, dir_GRmax_plot))

###########GR50 heatmap plot############

print("\nStarting GR50 heatmap plotting\n")

rscript = os.path.join(drug_response_scripts, "heatmap3_new.R")
os.system('/usr/bin/python %s/heatmap_for_GR50_median_center.py --metrics %s --c %s --w %s --heatmap_script %s' % (drug_response_scripts,s_GRmetrics, args.classification, dir_stat_plot, rscript))

print('/usr/bin/python %s/heatmap_for_GR50_median_center.py --metrics %s --c %s --w %s --heatmap_script %s' % (drug_response_scripts,s_GRmetrics, args.classification, dir_stat_plot, rscript))
print("The End")


