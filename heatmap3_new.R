#!/usr/bin/R

#load packages
library(argparse)
library("gplots")
library("devtools")

#Load latest version of heatmap.3 function
source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#######median center#######
## by row
#medianCenter <- function(y){
#  y <- apply(y, 1, function(x){x-median(x)})
#  scale(t(y), center=FALSE)
#}
## by col
#medianCenter <- function(y){
#  y <- apply(y, 2, function(x){x-median(x)})
#  scale(y, center=FALSE)
#}

#zoom limit for the color 
zlim=function(x, lim) { 
  x[which(x < min(lim))] = min(lim)
  x[which(x > max(lim))] = max(lim)
  return(x)
}

#argpase
parser <- ArgumentParser(description='Process some integers')
parser$add_argument('--data', help="input the data matrix")
parser$add_argument('--t', help="input the title of color key")
parser$add_argument('--r_side_label', help='input the lable of row side colorbar')
parser$add_argument('--c_side_label', help='input the lable of column side colorbar')
parser$add_argument('--w', type="double", help='input the width of the heatmap')
parser$add_argument('--h', type="double", help='input the height of the heatmap')
parser$add_argument('--min', type="double", help='input the min value of zoom limit')
parser$add_argument('--max', type="double", help='input the max value of zoom limit')
parser$add_argument('--out', help='out the heatmap')
args <- parser$parse_args()

#set output
pdf(file=args$out, width=args$w, height=args$h, family="ArialMT")

#read data
DATA=read.table(args$data, sep = "\t", header=T, row.names=1)
##median center
#median_norm=medianCenter(DATA)
#zoom limit
limit_DATA=zlim(data.matrix(DATA),c(args$min,args$max))

#define custom dist and hclust functions
mydist=function(x) {dist(x,method="euclidian")}
myclust=function(x) {hclust(x,method="average")}

#create color for cell lines
CellLine=read.table(args$c_side_label)
v_CellLine=as.vector(t(CellLine))
v_CellLine_uniq=unique(v_CellLine)

#create color for cell lines
Gene=read.table(args$r_side_label)
v_Gene=as.vector(t(Gene))
v_Gene_uniq=unique(v_Gene)

#combine the Cell line and gene uniq items
v_uniq=unique(v_CellLine_uniq)
v_uniq_sort=sort(v_uniq)

#define the color for cell line 
v_CellLine[v_CellLine=="Basal"]="red"
v_CellLine[v_CellLine=="Claudin-low"]="blue"
v_CellLine[v_CellLine=="ERBB2AMP"]="orange"
v_CellLine[v_CellLine=="Luminal"]="purple"
v_CellLine[v_CellLine=="Non-malignant"]="green"
CellLineColor=matrix(v_CellLine)

#define the color for gene class 
v_Gene[v_Gene=="Compound"]="yellow"
GeneColor=t(v_Gene)

#rownames(GeneColor)=c("Gene Class")
#colnames(CellLineColor)=c("Cell Lines")
print (nrow(CellLineColor))
#draw the heatmap
heatmap.3(limit_DATA,
        dendrogram="both",
        hclustfun=myclust,
        distfun=mydist,
        scale="none",
        col = colorpanel(100, "blue","gray","red"),
        #RowSideColors=GeneColor,
        #RowSideColorsSize=0.8,
        ColSideColors=CellLineColor,
        ColSideColorsSize=1,
        margins = c(8,14),
        cexRow = 1.2, 
        cexCol = 1.2,
        key=TRUE,
        keysize=0.9,
        KeyValueName=args$t,
        trace=c("none"),
)

legend("topright", legend=v_uniq_sort, fill=c('red', 'blue', 'orange', 'purple', 'green'), border=FALSE, bty="n", y.intersp = 0.8, cex=0.8) #fill=v_color[1:length(v_uniq_sort)]
dev.off()

