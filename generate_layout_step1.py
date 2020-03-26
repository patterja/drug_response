#!/usr/local/bin/python3

#program: generate_layout.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Tue May 30 10:41:31 2017

#description

#import libraries
import argparse
import csv

parser = argparse.ArgumentParser(description='')
parser.add_argument('--i', help='input the layout_primary.txt')
parser.add_argument('--o', help='output the first step layout file')

args = parser.parse_args()

#open the output file
with open (args.o, 'w') as outfile:
    h_out=csv.writer(outfile, delimiter="\t")
    lst_header_out=["PlateIndex", 'WellNumber', 'Compound', 'Row', 'Column']
    h_out.writerow(lst_header_out)
    #open the input file
    with open (args.i, "r") as infile:
        h_in=csv.reader(infile, delimiter="\t")
        lst_column_num=next(h_in)
        row_num=0
        PlateIndex=0
        for lst_row in h_in:
            row_num+=1
            row=lst_row[0]
            lst_drugs=lst_row[1:]
            column_num=0
            for s_drug in lst_drugs:
                column_num+=1
                PlateIndex+=1
                WellNumber=str(row) + str(column_num)
                Compound=s_drug
                lst_out=[PlateIndex, WellNumber, Compound, row_num, column_num]
                h_out.writerow(lst_out)
