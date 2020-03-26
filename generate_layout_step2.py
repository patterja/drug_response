#!/usr/local/bin/python3

#program: generate_layout.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Tue May 30 10:41:31 2017

#description

#import libraries
import argparse
import csv

def generate_new_layout(outfile, d_drugs):
    #output the layout file of plate3
    with open (outfile, 'w') as OUT:
        h_out=csv.writer(OUT, delimiter="\t")
        lst_header_out=["PlateIndex", 'WellNumber', 'DrugNumber', 'Compound', 'Concentration', 'Unit', 'Row', 'Column']
        h_out.writerow(lst_header_out)
        #read in the layout1.txt
        with open (args.i, "r") as i_file:
            h_i=csv.reader(i_file, delimiter="\t")
            lst_i_header=next(h_i)
            for lst_row in h_i:
                if(lst_row[2] == "Media" or lst_row[2] == 'DMSO' or lst_row[2] == 'PBS'):
                    drug=lst_row[2]
                    conc=""
                    drug_num=""
                    unit="nM"
                    lst_out=[lst_row[0], lst_row[1], drug_num, drug, conc, unit, lst_row[3], lst_row[4]]
                else:
                    lst_drug_info=lst_row[2].split("_")
                    drug=d_drugs[lst_drug_info[1]]['Drug']
                    conc=d_drugs[lst_drug_info[1]][lst_drug_info[2]]
                    drug_num=lst_drug_info[2]
                    unit="nM"
                    lst_out=[lst_row[0], lst_row[1], drug_num, drug, conc, unit, lst_row[3], lst_row[4]]
                h_out.writerow(lst_out)

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--i', help='input the layout_first_step.txt')
    parser.add_argument('--i2', help='input the drug_info.txt')
    parser.add_argument('--w', help='input the path of work directory')
    
    args = parser.parse_args()
    
    #read in the drug information file
    d_p1_drugs={}
    d_p2_drugs={}
    d_p3_drugs={}
    with open (args.i2, 'r') as i2_file:
        h_i2=csv.reader(i2_file, delimiter="\t")
        lst_i2_header=next(h_i2)
        for lst_row in h_i2:
            if(lst_row[1] == 'P1'):
                d_p1_drugs[lst_row[2]]={}
                d_p1_drugs[lst_row[2]].update({'Drug':lst_row[0]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[3]:lst_row[3]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[4]:lst_row[4]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[5]:lst_row[5]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[6]:lst_row[6]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[7]:lst_row[7]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[8]:lst_row[8]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[9]:lst_row[9]})
                d_p1_drugs[lst_row[2]].update({lst_i2_header[10]:lst_row[10]})
            if(lst_row[1] == 'P2'):
                d_p2_drugs[lst_row[2]]={}
                d_p2_drugs[lst_row[2]].update({'Drug':lst_row[0]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[3]:lst_row[3]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[4]:lst_row[4]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[5]:lst_row[5]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[6]:lst_row[6]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[7]:lst_row[7]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[8]:lst_row[8]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[9]:lst_row[9]})
                d_p2_drugs[lst_row[2]].update({lst_i2_header[10]:lst_row[10]})
            if(lst_row[1] == 'P3'):
                d_p3_drugs[lst_row[2]]={}
                d_p3_drugs[lst_row[2]].update({'Drug':lst_row[0]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[3]:lst_row[3]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[4]:lst_row[4]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[5]:lst_row[5]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[6]:lst_row[6]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[7]:lst_row[7]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[8]:lst_row[8]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[9]:lst_row[9]})
                d_p3_drugs[lst_row[2]].update({lst_i2_header[10]:lst_row[10]})

    #write down the final layout file
    if (len(d_p1_drugs) >0):
        f_p1_layout=args.w + "/layout_p1.txt"
        generate_new_layout(outfile=f_p1_layout, d_drugs=d_p1_drugs)  
    if (len(d_p2_drugs) >0):
        f_p2_layout=args.w + "/layout_p2.txt"
        generate_new_layout(outfile=f_p2_layout, d_drugs=d_p2_drugs)  
    if (len(d_p3_drugs) >0):
        f_p3_layout=args.w + "/layout_p3.txt"
        generate_new_layout(outfile=f_p3_layout, d_drugs=d_p3_drugs)
