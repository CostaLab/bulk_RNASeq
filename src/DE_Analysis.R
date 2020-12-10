library(limma, quietly=TRUE)
library(edgeR, quietly=TRUE)
library(RColorBrewer)
library(FactoMineR)
library(optparse)
source("helper_functions.R")  #import helper function limma_pipeline

# example usage: Rscript DE_Analysis.R --fo csv -c ../testdata/test_counts.csv -d ../testdata/data.csv -o ../test/ --vs F-E

option_list = list(
  make_option(c("-c", "--counts"), type="character", default="", help="file containing the gene counts with gene identifiers in first column", metavar="character"),
  make_option(c("-d", "--data"), type="character", default="", help="file containing information about the samples with sample identifiers in first column", metavar="character"),
  make_option("--fi", type="logical", default=TRUE, help="Determine whether lowly expressed genes should be filtered out before DE analysis"),
  make_option("--fo", type="character", default="tsv", help="format of the input files. Choose one from: \"tsv\", \"csv\"", metavar="character"),
  make_option(c("-g", "--groupcolumn"), type="integer", default=2, help="number of the experiment dataframe column which indicates the group"),
  make_option(c("-n", "--normalize"), type="character", default="TMM", help="the normalization that should be applied before DE analysis. Choose one from: \"none\", \"TMM\", \"TMMwsp\", \"RLE\" and \"upperquartile\""),
  make_option(c("-o", "--output"), type="character", default="./", help="directory for output files", metavar="character"),
  make_option(c("-s", "--significance"), type="double", default=0.05, help="Significance threshold for adjusted p-Value (significant DE genes are those having an adj. p-Value below this threshold)"),
  make_option(c("-t", "--tool"), type="character", default="limma", help="name of the tool you want to use. Choose one from: \"DESeq2\", \"limma\", \"limma_trend\" and \"limma_voom\"", metavar="character"),
  make_option(c("-u", "--unit"), type="logical", default=FALSE, help="Determine whether to scale to unit variance when plotting the counts"),
  make_option("--vs", type="character", default="", help="Determine the conditions which should be compared (used for building contrast matrix) - usually: treatment - wt", metavar="character")
)

# TODO maybe add possibility to directly include gene information in output file later

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check tool parameter
print(">> check passed parameters", quote=FALSE)
valid_tools = c("DESeq2", "limma", "limma_trend", "limma_voom")
valid_formats = c("tsv", "csv")

# exit script if invalif tool or input file format is passed
if(! opt$tool %in% valid_tools) {
  stop("Invalid tool passed. Please choose one from: \"DESeq2\", \"limma\", \"limma_trend\", \"limma_voom\"", quote=FALSE)
}
if(! opt$fo %in% valid_formats){
  stop("Invalid input file format passed. Please choose one from:  \"tsv\", \"csv\"", quote=FALSE)
}

# import data
print(">> Import data", quote=FALSE)
if (opt$fo == "tsv") {
  count_data = read.table(file=opt$counts, sep="\t", row.names=1, header=None)
  exp_data = read.table(file=opt$data, sep="\t", row.names=1, header=None)
} else {
  count_data = read.csv(file=opt$counts, sep=",", row.names=1)
  exp_data = read.csv(file=opt$data, sep=",", row.names=1)
}

# get group labels
labels = as.character(exp_data[,opt$groupcolumn])

# plot raw counts (PCA)
print(">> PCA Analysis", quote=FALSE)
PCA_counts <- as.data.frame(t(count_data))  # count matrix needs to be transposed before PCA
PCA_counts <- cbind(exp_data, PCA_counts)  # enrich counts with information about samples
cairo_pdf(file=paste(c(opt$output, "raw_counts_PCA.pdf"), sep="", collapse=""), width=11.43, height=8.27, onefile=T)
res.PCA = PCA(PCA_counts, scale.unit=opt$unit, quali.sup=seq(1,ncol(exp_data)), graph=F)
plot(res.PCA, cex=0.9, habillage=opt$groupcolumn, col.hab=c("black","blue"), invisible="quali", unselect=0, axes=c(1,2), title="PC1 vs. PC2 raw counts")
dump <- dev.off()

# call tool pipeline
# limma
if (opt$tool=="limma"){
  print(">> Start analysis with limma", quote=FALSE)
  result = limma_pipeline(expression_matrix = count_data, group_labels = labels, group_contrasts = opt$vs, significance_threshold=opt$significance, reference_group_labels=NULL,
    removeIntercept=TRUE, filter=opt$fi, normalize=opt$normalize, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
}
#limma trend
else if (opt$tool=="limma_trend"){
  print(">> Start analysis with limma trend", quote=FALSE)
  result = limma_pipeline(expression_matrix = count_data, group_labels = labels, group_contrasts = opt$vs, significance_threshold=opt$significance, reference_group_labels=NULL,
    removeIntercept=TRUE, filter=opt$fi, normalize=opt$normalize, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
}
# limma voom
else if (opt$tool=="limma_voom"){
  print(">> Start analysis with limma voom", quote=FALSE)
  sample_weights = askYesNo(msg="Do you want to apply sample specific weights?")
  if (sample_weights){
    result = limma_pipeline(expression_matrix = count_data, group_labels = labels , group_contrasts = opt$vs, significance_threshold=opt$significance, reference_group_labels=NULL,
      removeIntercept=TRUE, filter=opt$fi, normalize=opt$normalize, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
  }
  else {
    result = limma_pipeline(expression_matrix = count_data, group_labels = labels, group_contrasts = opt$vs, significance_threshold=opt$significance, reference_group_labels=NULL,
      removeIntercept=TRUE, filter=opt$fi, normalize=opt$normalize, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
  }
}
# DESeq2
else {
  print(">> Start analysis with DESeq2", quote=FALSE)
  # TODO put DESeq2 pipeline call here
}

# output results
# TODO
print(">> finished analysis", quote=FALSE)
