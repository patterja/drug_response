#!/usr/local/bin/python3

#program: calculate_GR_at_serum_conc.py
#Author: Hongmei Zhang
#email:zhanho@ohsu.edu
#data:Thu Jul 20 12:43:08 2017

#description
#calculate the Growth ratio at serum concentration, and add it the GRmetrics

#import libraries
import argparse
from gr50 import logistic, _fit, _flat, _calculate_pval
import pandas as pd
import numpy as np

def GR_of_serum(df, alpha, serum_con):
    conc_min = df.concentration.min() / 100
    conc_max = df.concentration.max() * 100
    bounds = np.array([[-1, 1], np.log10([conc_min, conc_max]), [0.1, 5]])
    prior = np.array([0.1, np.log10(np.median(df.concentration)), 2])
    logistic_result = _fit(logistic, df.concentration, df.GRvalue,
                           prior, bounds)
    flat_result = _fit(_flat, df.concentration, df.GRvalue,
                       prior[[0]], bounds[[0]])
    pval = _calculate_pval(logistic_result, flat_result, len(df.concentration))
    if pval > alpha or not logistic_result.success:
        #set GRserum to max GR value
        GRserum=max(df.GRvalue)
        GRserum=1.0
    else:
        GRserum = logistic(serum_con, logistic_result.x)
    
    return [GRserum]

def mklist(values):
    """Convert tuple to list, and anything else to a list with just that thing.

    This is a helper to fix an inconsistency with the group keys in a
    pandas.groupby object. When grouping by multiple columns, the keys are
    tuples of values. When grouping by a single column, even if specified as a
    single-element list, the keys are single values. We always want a list,
    which is what this function accomplishes."""
    if isinstance(values, tuple):
        return list(values)
    else:
        return [values]

def serum_metrics(data, alpha, df_serum_conc):
    non_keys = set(('concentration', 'cell_count', 'cell_count__ctrl',
                    'cell_count__time', 'GRvalue'))
    metric_columns = ['GRserum']
    keys = list(set(data.columns) - non_keys)
    drug_index=keys.index("agent")
    gb = data.groupby(keys)
    data = [mklist(k) + GR_of_serum(v, alpha, df_serum_conc.ix[k[drug_index], "Serum concentration"]) for k, v in gb]
    df_serum = pd.DataFrame(data, columns=keys + metric_columns)
    return df_serum


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--m', help='input the GRmetrics')
    parser.add_argument('--gr', help='input the GRvalues')
    parser.add_argument('--s', help='input the serum concentration file')
    parser.add_argument('--o', help='output the new matrix of GRmetrics')
    
    args = parser.parse_args()

    #read in the table of serum concentration of drugs
    df_serum_conc=pd.read_csv(args.s, sep="\t", index_col=0, header=0)
    #read in the matrix of GRvalues
    df_GR=pd.read_csv(args.gr, sep="\t", header=0)
    #calculate the GR values for drug at serum concentration
    df_serum_GR=serum_metrics(data=df_GR, alpha=0.05, df_serum_conc=df_serum_conc)
    df_serum_GR=df_serum_GR.set_index(["agent", "timepoint", "cell_line"])

    #read in the GRmetrics
    df_GRmetrics=pd.read_csv(args.m, sep="\t", header=0, index_col=None)
    df_GRmetrics=df_GRmetrics.set_index(["agent", "timepoint", "cell_line"])
    df_merge=pd.merge(df_GRmetrics, df_serum_GR, left_index=True, right_index=True)

    df_merge.to_csv(args.o, sep="\t")
