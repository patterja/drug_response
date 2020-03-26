#program: generate_boxplot_table_for_cv.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Fri Mar 10 10:06:18 2017

#description

#import libraries
import argparse
import csv

def generate_boxplot_table_for_cv(infile):
    #read in the data
    d_cv={}
    with open (infile, "r") as tsv_in:
        h_in=csv.reader(tsv_in, delimiter="\t")
        #header=next(h_in)
        num=0
        old_plate=''
        i_plateNum=0
        for lst_row in h_in:
            num+=1
            (plate, cv)=lst_row
            if (plate != old_plate):
                num=1
                try:
                   d_cv[i_plateNum].append(plate)
                   d_cv[num].append(cv)
                except KeyError:
                   d_cv.update({i_plateNum: [plate]})
                   d_cv.update({num: [cv]})
            else:
                try:
                   d_cv[num].append(cv)
                except KeyError:
                   d_cv.update({num: [cv]})
            old_plate=plate
    #return the dictornary
    return d_cv

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the variance data of all cell lines')
    parser.add_argument('--o', help='output the unstack variance data')
    
    args = parser.parse_args()
    
    #read in the raw table and get the dictornary of cv
    d_cv=generate_boxplot_table_for_cv(args.i) 
    #write out
    with open (args.o, "w")  as tsv_out:
        h_out=csv.writer(tsv_out, delimiter="\t")
        #print (sorted(d_cv.keys()))
        for key in sorted(d_cv.keys()):
            h_out.writerow(d_cv[key])
