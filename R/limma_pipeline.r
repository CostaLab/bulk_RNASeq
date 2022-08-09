#' limma pipeline function
#'
#' @param expression_matrix A matrix with features as rows and samples as columns
#' @param group_labels TODO
#' @param group_contrasts TODO
#' @param significance_threshold TODO
#' @param reference_group_labels TODO
#' @param removeIntercept TODO
#' @param include_unadj_pval TODO
#' @param filter TODO
#' @param normalize TODO
#' @param limma_trend TODO
#' @param limma_voom TODO
#' @param limma_voom_weight_samples TODO
#' @param output_dir TODO
#'
#' @export
limma_pipeline <- function(
  expression_matrix,
  group_labels,
  group_contrasts,
  significance_threshold=0.05,
  reference_group_labels=NULL,
  removeIntercept=FALSE,
  include_unadj_pval=FALSE,
  filter=FALSE,
  normalize="none",
  limma_trend=FALSE,
  limma_voom=FALSE,
  limma_voom_weight_samples=FALSE,
  output_dir="./"
){
  # features as rows, samples as columns

  # term to skip intercept col index later on
  skip_I <- 0
  group_labels <- factor(group_labels)
  if(!is.null(reference_group_labels)){
    group_labels <- relevel(group_labels, ref=reference_group_labels)
  }

  # TODO should probably always remove intercept...
  print(">>> Create design matrix", quote=FALSE)
  if(removeIntercept){
    # build design matrix
    design <- model.matrix(~group_labels-1) # NOTE removing intercept!
    colnames(design) <- levels(group_labels)
    skip_I <- 0
  } else {
    design <- model.matrix(~group_labels)
    colnames(design) <- c(colnames(design)[1],levels(group_labels)[-1])
    skip_I <- 1
  }

  # create DGEList object
  dge <- edgeR::DGEList(counts=expression_matrix)

  # filter lowly expressed genes if this is specified
  if(filter){
    print(">>> Filter out lowly expressed genes", quote=FALSE)
    keep <- edgeR::filterByExpr(dge, design)
    print(paste(c(">>>> Number of genes filtered out: ", length(keep[keep == FALSE])), sep="", collapse=""), quote=FALSE)
    print(paste(c(">>>> Number of genes still considerd: ", length(keep[keep == TRUE])), sep="", collapse=""), quote=FALSE)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    # TODO PCA plot after filtering
  }

  # TODO add possibility to normalize (also plot normalized counts if this option is chosen)
  print(">>> Normalize Counts", quote=FALSE)
  dge <- edgeR::calcNormFactors(dge, method=normalize)
  if(normalize!="none"){
    # TODO PCA plot after normalization
  }

  print(">>> Fit model", quote=FALSE)
  if(removeIntercept){
    # voom was shown to have the edge when the sequencing depths
    # were very inconsistent between replicates but, otherwise,
    # limma-trend was just as good
    if(limma_voom){
      vv <- limma::voom(expression_matrix,design,plot=TRUE)
      # TODO plot tranformed counts
      fit <- limma::lmFit(vv,design)
      limma_trend=FALSE # if using limma-voom don't use limma-trend
    } else if (limma_voom_weight_samples){
      vv <- limma::voomWithQualityWeights(expression_matrix,design,normalization="none",plot=TRUE)
      # TODO plot transformed counts
      fit <- limma::lmFit(vv,design)
      limma_trend=FALSE # if using limma-voom don't use limma-trend
    } else {
      fit <- limma::lmFit(expression_matrix,design)
    }
    contrast_matrix <- limma::makeContrasts(
      contrasts=group_contrasts,
      levels=levels(group_labels)
    )

    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2,trend=limma_trend)
    
    grDevices::cairo_pdf(file=paste(c(output_dir, "Model_plots.pdf"), sep="", collapse=""), width=11.43, height=8.27, onefile=T)
    limma::plotSA(fit2, main="Final model: Mean-variance trend")
    dump <- grDevices::dev.off()

  } else {
    fit <- limma::lmFit(expression_matrix,design)
    fit2 <- limma::eBayes(fit,trend=limma_trend)
  }

  print(">>> Evaluate results", quote=FALSE)
  limma_res <- list()

  # TODO: this but without using a for loop
  for(i in seq_len(ncol(fit2$coefficients)-skip_I)){

    i <- i + skip_I
    # NOTE not applying FC thresholding, not filtering out any genes/probes
    results <- limma::decideTests(
      fit2,
      coef=i,
      p.value=significance_threshold,
      adjust.method="BH",
      lfc=0
    )
    # see notes from documentation ie
    # "Users wanting to use fold change
    # thresholding are usually recommended
    # to use treat and topTreat instead"
    aux <- limma::topTable(
      fit2,
      coef=i,
      number="Inf",
      adjust.method="BH",
      sort.by="none"
    )

    contrast_name <- colnames(fit2$coefficients)[i]
    contrast_name_fc <- paste0(contrast_name,"_FC")
    contrast_name_pval <- paste0(contrast_name,"_Pval")
    contrast_name_de <- paste0(contrast_name,"_DE")
    
    limma_res[[contrast_name_fc]] <- aux[,"logFC"]
    limma_res[[contrast_name_pval]] <- aux[,"adj.P.Val"]

    # -1, 0 or 1 depending on whether each t-statistic is classified as
    # significantly negative, not significant or significantly positive, respectively
    limma_res[[contrast_name_de]] <- results@.Data[,i]

    # ID col only appears in the case of duplicated rownames in fit
    if("ID" %in% colnames(aux)){
      limma_res[["Gene_Symbol"]] <- aux[,"ID"]
    } else {
      limma_res[["Gene_Symbol"]] <- rownames(aux)
    }
    
    if(include_unadj_pval){
      contrast_name_unadjpval <- paste0(contrast_name,"_unadjPval")
      limma_res[[contrast_name_unadjpval]] <- aux[,"P.Value"]
    }
  }

  limma_res <- as.data.frame(limma_res)
  # limma_res$Gene_Symbol <- rownames(limma_res)

  return(list(
    topGenes=limma_res,
    fit=fit2
  ))
}


