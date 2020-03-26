#program: generate_input_for_gr50_tools.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Tue Feb 21 16:32:02 2017

#description
#this program will generate the input file for gr50_tools

#import libraries
import argparse
import csv
import os
import statistics
import math
import sys

def time0_median(plateID, wells, path):
    filename=os.path.join(path, plateID + ".txt")
    if (os.path.exists(filename)):
        with open (filename, "r") as tsvfile:
            h_t0=csv.reader(tsvfile, delimiter="\t")
            header=next(h_t0)
            lst_count=[]
            for lst_row in h_t0:
                (well, count)=lst_row[1:3]
                count=float(count)
                if well in wells:
                    lst_count.append(count)
            print (plateID)
            print (lst_count)
            median=statistics.median(lst_count)
            return (median)
    else:
        print ("{} not exists".format(filename))
        return None

def time0_mean_filter(plateID, wells, path):
    filename=os.path.join(path, plateID + ".txt")
    if (os.path.exists(filename)):
        with open (filename, "r") as tsvfile:
            h_t0=csv.reader(tsvfile, delimiter="\t")
            header=next(h_t0)
            lst_count=[]
            for lst_row in h_t0:
                (well, count)=lst_row[1:3]
                count=float(count)
                if well in wells:
                    lst_count.append(count)
            mean_raw=statistics.mean(lst_count)
            #filter the data
            for value in lst_count:
                if (abs(math.log10(value) - math.log10(mean_raw))>1):
                    lst_count.remove(value)
            mean_filter=statistics.mean(lst_count)
            return (mean_filter)
    else:
        print ("{} not exists".format(filename))
        return None

def generate_input(fname, testTime,  t0_median, CL, num, d_out):
        #open the input handle and find the cell count of negtive control
        #PlateID Compound        Concentration   Unit    CellCount
        #SM300049         BYL719 160     nM      1759.0
        f_cc_ctrl=0.0
        ctrl_num=0.0
        with open (fname, "r")  as tsvin:
            h_in=csv.reader(tsvin, delimiter="\t")
            for lst_row in h_in:
                if (lst_row[1] == "PBS" or lst_row[1] == "DMSO"):
                    f_cc_ctrl+=float(lst_row[4])
                    ctrl_num+=1
            if(ctrl_num ==0):
                sys.exit("No value for PBS or DMSO")
            else:
                f_cc_ctrl=f_cc_ctrl/ctrl_num
        #read in the drug data and write down   
        with open (fname, "r")  as tsvin:
            h_in=csv.reader(tsvin, delimiter="\t")
            header=next(h_in)
            for lst_row in h_in:
                if (lst_row[1] not in ["DMSO", 'PBS']):
                    lst_out=[CL, lst_row[1], testTime, lst_row[2], lst_row[4], f_cc_ctrl, t0_median]
                    num+=1
                    try:
                        d_out[num]=lst_out
                    except KeyError:
                        d_out.update(num, lst_out)
        return(num, d_out)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--d0', help='input the path of time 0 paltes')
    parser.add_argument('--d', help='input the path of drug paltes')
    parser.add_argument('--s', help='input the sample information file')
    parser.add_argument('--o', help='output the merge file')
    
    args = parser.parse_args()
    
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
            t0_median=time0_median(plateID=time0_plateID, wells=wells, path=args.d0)
            drug_file=os.path.join(args.d, drug_plateID + "_CellCount_median.txt")
            if (os.path.exists(drug_file)):
                (i_num, d_output)=generate_input(fname=drug_file, testTime=time, t0_median=t0_median, CL=cellLine, num=i_num, d_out=d_output)
            else:
                print ("{} not exists".format(drug_file))
    ########write out the result#########  
    #cell_line	agent	perturbation	replicate	time	concentration	cell_count	cell_count__ctrl	cell_count__time0
    #MCF10A	drugA	0	1	48	0.001	1131	1212.5	299.5
    with open (args.o, "w") as tsvout:
        h_out=csv.writer(tsvout, delimiter="\t")
        lst_header=["cell_line", "agent", "time", "concentration", "cell_count", "cell_count__ctrl", "cell_count__time0"]
        h_out.writerow(lst_header)
        for key in d_output.keys():
            h_out.writerow(d_output[key])


