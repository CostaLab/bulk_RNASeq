library(limma, quietly=TRUE)
library(edgeR, quietly=TRUE)
library(RColorBrewer)
library(FactoMineR)
library(optparse)
source("helper_functions.R")  #import helper function limma_pipeline

option_list = list(
  make_option(c("-t", "--tool"), type="character", default="limma",
              help="name of the tool you want to use. Choose one from: \"DESeq2\", \"limma\", \"limma_trend\", \"limma_voom\"", metavar="character"),
  make_option(c("-f", "--format"), type="character", default="tsv", help="format of the input files. Choose one from: \"tsv\", \"csv\"", metavar="character"),
  make_option(c("-c", "--counts"), type="character", default="", help="file containing the gene counts", metavar="character"),
  make_option(c("-d", "--data"), type="character", default="", help="file containing information about the samples", metavar="character"),
  make_option(c("-g", "--groupcolumn"), type="integer", default=2, help="number of the experiment dataframe column which indicates the group"),
  make_option(c("-o", "--output"), type="character", default="./", help="directory for output files", metavar="character"),
  make_option(c("-v", "--contrast"), type="character", default="",
              help="Determine the conditions which should be compared (used for building contrast matrix) - usually: treatment - wt", metavar="character"),
  make_option(c("-s", "--significance"), type="double", default=0.05,
              help="Significance threshold for adjusted p-Value (significant DE genes are those having an adj. p-Value below this threshold)"),
  make_option(c("-u", "--unit"), type="logical", default=FALSE, help="Determine whether to scale to unit variance when plotting the counts")
)

# TODO maybe add possibility to directly include gene information in output file later

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check tool parameter
print(">> check passed parameters", quote=FALSE)
valid_tools = c("DESeq2", "limma", "limma_trend", "limma_voom")
valid_formats = c("tsv", "csv")
if(! opt$tool %in% valid_tools) {
  print("Invalid tool passed. Please choose one from: \"DESeq2\", \"limma\", \"limma_trend\", \"limma_voom\"", quote=FALSE)
} else if(! opt$format %in% valid_formats){
  print("Invalid input file format passed. Please choose one from:  \"tsv\", \"csv\"", quote=FALSE)
} else {

  # import data
  print(">> Import data", quote=FALSE)
  if (opt$format == "tsv") {
    count_data = read.table(file=opt$counts, sep="\t", row.names=1, header=None)
    exp_data = read.table(file=opt$data, sep="\t", row.names=1, header=None)
  }
  else {
    count_data = read.csv(file=opt$counts, sep=",", row.names=1)
    exp_data = read.csv(file=opt$data, sep=",", row.names=1)
  }

  # get group labels
  labels = as.character(exp_data[,opt$groupcolumn-1])

  # plot raw counts (PCA)
  print(">> PCA Analysis", quote=FALSE)
  PCA_counts <- as.data.frame(t(count_data))  # count matrix needs to be transposed before PCA
  PCA_counts <- cbind(exp_data, PCA_counts)  # enrich counts with information about samples
  cairo_pdf(file=paste(c(opt$output, "raw_counts_PCA.pdf"), sep="", collapse=""), width=11.43, height=8.27, onefile=T)
  res.PCA = PCA(PCA_counts, scale.unit=opt$unit, quali.sup=seq(1,ncol(exp_data)), graph=T)
  plot(res.PCA, cex=0.9, habillage=opt$groupcolumn, col.hab=c("black","blue"), invisible="quali", unselect=0, axes=c(1,2), title="PC1 vs. PC2 raw counts")
  dev.off()

  # call tool pipeline
  # limma
  if (opt$tool=="limma"){
    print(">> Start analysis with limma", quote=FALSE)
    limma_pipeline(expression_matrix = count_data, group_labels = labels, group_contrasts = opt$contrast, significance_threshold=opt$significance, reference_group_labels=NULL,
      removeIntercept=TRUE, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
  }
  #limma trend
  else if (opt$tool=="limma_trend"){
    print(">> Start analysis with limma trend", quote=FALSE)
    limma_pipeline(expression_matrix = count_data, group_labels = labels, group_contrasts = opt$contrast, significance_threshold=opt$significance, reference_group_labels=NULL,
      removeIntercept=TRUE, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
  }
  # limma voom
  else if (opt$tool=="limma_voom"){
    print(">> Start analysis with limma voom", quote=FALSE)
    sample_weights = askYesNo(msg="Do you want to apply sample specific weights?")
    if (sample_weights){
      limma_pipeline(expression_matrix = count_data, group_labels = labels , group_contrasts = opt$contrast, significance_threshold=opt$significance, reference_group_labels=NULL,
        removeIntercept=TRUE, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
    }
    else {
      limma_pipeline(expression_matrix = count_data, group_labels = labels, group_contrasts = opt$contrast, significance_threshold=opt$significance, reference_group_labels=NULL,
        removeIntercept=TRUE, include_unadj_pval=FALSE, limma_trend=FALSE, limma_voom=FALSE, limma_voom_weight_samples=FALSE)
    }
  }
  # DESeq2
  else {
    print(">> Start analysis with DESeq2", quote=FALSE)
    # TODO put DESeq2 pipeline call here
  }

  # output results
  # TODO
}
