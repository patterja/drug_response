
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gr50
import argparse

parser = argparse.ArgumentParser(description='')
parser.add_argument('--v', help='input the GRvalue file')
parser.add_argument('--m', help='input the GRmetrics file')
parser.add_argument('--d', help='input the work directory for the GR50 barplot')

args = parser.parse_args()

# Compute the GR metrics from the data.
df_data=pd.read_csv(args.m, sep='\t')

# Produce a trellis plot showing the fitted curves and some of the metrics
# across the different cell lines and drugs.
sns.set(style="ticks")
sns.set(font_scale=1.3)
grid = sns.FacetGrid(df_data, col="agent", col_wrap=5)
grid.set(yscale="log")
grid.map(plt.hist, "GR50", color='red')
        #ax.hlines(0, x_min, x_max, '#707070', lw=1.0)
        ##ax.hlines(m.GRinf, x_min, x_max, '#ff00ff', linestyles='dashed',lw=1.0)
        #ax.vlines(m.GR50, -1, 1, 'b', linestyles='dashed', lw=1.0)
        ##ax.vlines(m.GEC50, -1, 1, '#CC6600', linestyles='dashed', lw=1.0)
        #ax.plot(fit_x, fit_y, 'r', lw=1.5)
grid.set(ylim=(-1, 1.5))
grid.set(xlim=(x_min, x_max))
figfile=args.d + "/" + s_agent + "_GR_plot.pdf"
