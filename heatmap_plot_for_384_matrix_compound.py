#program: generate_384_well_layout_table.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Wed Mar 15 10:10:35 2017

#description

#import libraries
import argparse
import pandas as pd
import csv
import seaborn as sns
import matplotlib.pyplot as plt
import math

def generate_conc_table(s_layout):
    df_layout=pd.read_csv(s_layout, sep="\t")
    df_conc=df_layout[['PlateIndex', 'Concentration']]
    df_conc=df_conc.set_index('PlateIndex')
    d_conc=df_conc.to_dict(orient='index')
    d_conc_layout={}
    i_row=0
    for i in range(384):
        i_num=i+1
        i_rem=i_num%24
        if (d_conc[i_num]["Concentration"]>0):
            f_conc=math.log10(d_conc[i_num]["Concentration"])
        else:
            f_conc=-1
        if (i_rem==1):
            i_row+=1
            d_conc_layout.update({i_row: [f_conc]})
        else:
            d_conc_layout[i_row].append(f_conc)
    #convert the dictionary into dataframe
    df_conc_layout=pd.DataFrame.from_dict(d_conc_layout, orient='index')
    df_conc_layout.columns=range(1,25)
    df_conc_layout.index=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    #return the table
    return df_conc_layout

def generate_comp_table(s_layout):
    df_layout=pd.read_csv(s_layout, sep="\t")
    df_comp=df_layout[['PlateIndex', 'Compound']]
    df_comp=df_comp.set_index('PlateIndex')
    d_comp=df_comp.to_dict(orient='index')
    d_comp_layout={}
    i_row=0
    for i in range(384):
        i_num=i+1
        i_rem=i_num%24
        s_comp=d_comp[i_num]["Compound"]
        if (i_rem==1):
            i_row+=1
            d_comp_layout.update({i_row: [s_comp]})
        else:
            d_comp_layout[i_row].append(s_comp)
    #convert the dictionary into dataframe
    df_comp_layout=pd.DataFrame.from_dict(d_comp_layout, orient='index')
    df_comp_layout.columns=range(1,25)
    df_comp_layout.index=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    #return the table
    return df_comp_layout

def sns_heatmap(df_conc, df_comp):
    sns.set()
    fig=plt.figure(figsize=(18,8))
    ax=fig.add_subplot(1,1,1)
    sns.heatmap(df_conc, ax=ax, annot=df_comp, fmt='', cbar_kws={'label': 'Log10 (Drug Concentration)'}, annot_kws={"size":7}, linewidths=.5)
    legend='''
Abir: Abiraterone
Afin: Afinitor
BIBW: BIBW2992
BYL719: BYL719
Cabo: Cabozantinib
Carb: Carboplatin
DMSO: DMSO
Enza: Enzaludemine
Gemc: Gemcitabine
Ibup: Ibuprofen
ILL: Ibuprofen+Levoxyl+Lapatinib
Lenv : Lenvatinib 
Levo: Levoxyl
LOSI: Linsitinib OSI-906
Media: Media
Olap: Olaparib
PD: PD-0332991
PH: Pertuzumab+Herceptin
Pona: Ponatinib
Tram: Tramestinib
Tret: Tretinoin
Vori: Vorinostat
'''
    plt.text(28,4, legend, fontsize=12)
    return fig

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--l', help='input the layout.txt file')
    parser.add_argument('--o', help='output the figure')
    
    args = parser.parse_args()

    #convert the data into 384 well layout table
    df_conc_table=generate_conc_table(s_layout=args.l)
    df_comp_table=generate_comp_table(s_layout=args.l)
    #plot the heatmap
    fig=sns_heatmap(df_conc=df_conc_table, df_comp=df_comp_table)
    #save figure
    fig.savefig(args.o, bbox='tight')
