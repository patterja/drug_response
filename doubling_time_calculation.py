#program: doubling_time_calculation.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Tue Mar  7 12:50:01 2017

#description

#import libraries
import argparse
import csv
import math

parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', help='input the input_file_of_gr50_tools.txt file')
parser.add_argument('--o', help='output the doubling time of cell line')

args = parser.parse_args()

#open the output file handle
with open (args.o, "w") as tsvout:
    h_out=csv.writer(tsvout, delimiter="\t")
    h_out.writerow(["CellLine", "DoublingTime"])
    #open the input file
    with open (args.i, "r") as tsvin:
        h_in=csv.reader(tsvin, delimiter="\t")
        header=next(h_in)
        d_dt={}
        for lst_row in h_in:
            (cellLine, cn0, cnt, duration)=(lst_row[0], float(lst_row[6]), float(lst_row[5]), float(lst_row[2]))
            dt=(duration*math.log2(2))/(math.log2(cnt)-math.log2(cn0))
            dt='{:.2f}'.format(dt)
            try:
                d_dt[cellLine]=dt
            except  KeyError:
                d_dt.update({cellLine: dt})
            
    #write out
    for key in d_dt.keys():
        h_out.writerow([key, d_dt[key]])
