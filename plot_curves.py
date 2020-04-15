import argparse
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from gr50 import logistic


#This code is based on the plot_curves.py function found here https://github.com/datarail/gr_metrics
# It has been significantly altered from a function to runnable script that outputs data aggregated
# into specified plots.

VERSION = "0.1.0"


def supply_args():
    """
    Input arguments
    """
    parser = argparse.ArgumentParser(description='plotting curves')
    parser.add_argument('--metrics', type=str, help='GR_metrics.txt')
    parser.add_argument('--values', type=str,help='GR_values.txt')
    parser.add_argument('--outdir', help='input the work directory for curves figure')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args




def plot_curves(gr_metrics, gr_values, colorvar=None, colvar=None, rowvar=None):
    '''Produce a trellis plot showing the fitted curves and some of the metrics
     across the different cell lines and drugs.'''

    # padding with dummy columns to have general variables
    grm = gr_metrics.copy()
    grv = gr_values.copy()
    if rowvar is None:
        grm['_r'] = 0
        grv['_r'] = 0
        rowvar = '_r'
    if colvar is None:
        grm['_c'] = 0
        grv['_c'] = 0
        colvar = '_c'
    if colorvar is None:
        grm['_cl'] = 0
        grv['_cl'] = 0
        colorvar = '_cl'

    # display parameters
    colorkeys = grv[colorvar].drop_duplicates().values
    cmap = matplotlib.cm.get_cmap('nipy_spectral')
    colors = cmap(np.array(range(len(colorkeys)))/(1.*len(colorkeys)))
    colors = np.maximum(colors**1.3-.1,0)
    x_min = grv.concentration.min() / 10
    x_max = grv.concentration.max() * 10
    fit_x = np.logspace(np.log10(x_min), np.log10(x_max))

    #### construction of the plot grid
    sns.set(style="ticks")
    # two cases based on rowvar
    if rowvar is '_r':
        grid = sns.FacetGrid(grv, col_wrap=np.ceil(np.sqrt(len(grv[colvar].drop_duplicates()))),
                             col=colvar, margin_titles=True)
    else:
        grid = sns.FacetGrid(grv, row=rowvar, col=colvar, margin_titles=True)
    # format the axes
    grid.set(xscale="log")
    grid.set_axis_labels(x_var='Concentration (uM)', y_var='GR value')
    grid.set(ylim=(-1, 1.15))
    grid.set(xlim=(x_min, x_max*3))
    grid.fig.tight_layout(w_pad=1)
    # set the positions of the axis labels
    if rowvar is not '_r':
        if colvar is '_c':
            grid.set_titles(row_template='{row_var}={row_name}', col_template='')
        else:
            grid.set_titles(template='{row_var}={row_name}, {col_var}={col_name}',
                            row_template='{row_var}={row_name}',
                            col_template='{col_var}={col_name}')
    else:
        if colvar is not '_c':
            grid.set_titles(template='{col_var}={col_name}',row_template='')
    # small hack to get uniform structure if rowvar was not specified
    if rowvar is '_r':
        grid.row_names = [0]
        grid.axes[0] = grid.axes.copy()

    #### loop to plot the data
    for rowval, row_axes in zip(grid.row_names, grid.axes):
        for colval, ax in zip(grid.col_names, row_axes):
            ax.hlines(0, x_min, x_max, '#707070', lw=0.5)
            h = []
            for i,colorval in enumerate(colorkeys):
                htemp, = ax.plot(np.nan, np.nan, 'o-', ms=5, color=colors[i,], lw=1, label=colorval)
                h.append(htemp)
                v = grv[(grv[colvar] == colval) & (grv[colorvar] == colorval) &
                        (grv[rowvar] == rowval)]
                ax.plot(x_max*np.array([1.5, 3]), [v.GRvalue.min()]*2, '-', color=colors[i,], lw=2)
                #ax.plot(v.concentration, v.GRvalue, marker='o', ms=5, color=colors[i,], ls='None', label='_nolegend_')
                for m in grm[(grm[colvar] == colval) & (grm[colorvar] == colorval) &
                                (grm[rowvar] == rowval)].itertuples():
                    fit_y = logistic(fit_x, [m.GRinf, np.log10(m.GEC50 +1), m.h_GR])
                    ax.plot([m.GR50]*2, [-1, -.9], '-', color=colors[i,], lw=2)
                    ax.plot(fit_x, fit_y, lw=1, color=colors[i,])
                    ax.legend(bbox_to_anchor=(1, 0.5), loc=6, prop={'size': 2}, ncol=1, borderaxespad=3, markerscale=0.1)

    # legend (needs to be better positioned)

    #plt.legend(tuple(h), tuple(colorkeys), loc='upper left', prop={'size': 3})

    return plt.gcf()

def main():
    args = supply_args()
    gr_metrics = pd.read_csv(args.metrics, sep='\t')
    gr_values = pd.read_csv(args.values, sep='\t')
    #write by drug
    for drug in gr_metrics.agent.drop_duplicates().values:
        grm_subset = gr_metrics[gr_metrics['agent']==drug]
        grv_subset = gr_values[gr_values['agent']==drug]
        pltobj = plot_curves(grm_subset, grv_subset, colorvar="cell_line", colvar="agent",rowvar=None)
        pltobj.savefig(os.path.join(args.outdir, drug +"_GRplot.pdf"), bbox_inches='tight')
        #pltobj = plot_curves(grm_subset, grv_subset, colorvar=None, colvar="agent",rowvar="cell_line")
        plt.close()
        #write by cell line
    for cellline in gr_metrics.cell_line.drop_duplicates().values:
        grm_subset = gr_metrics[gr_metrics['cell_line']==cellline]
        grv_subset = gr_values[gr_values['cell_line']==cellline]
        #pltobj = plot_curves(grm_subset, grv_subset, colorvar=None, colvar="agent",rowvar="cell_line")
        pltobj = plot_curves(grm_subset, grv_subset, colorvar="agent", colvar="cell_line",rowvar=None)
        pltobj.savefig(os.path.join(args.outdir, cellline +"_GRplot.pdf"), bbox_inches='tight')
        plt.close()





        #giant plot of figures in a grid
    #pltobj = plot_curves(gr_metrics, gr_values, colorvar=None, colvar="agent", rowvar="cell_line")
    #pltobj.savefig('plot_curve.pdf')


if __name__ == "__main__":
    main()
