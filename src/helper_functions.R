
# Helper functions for genomics analysis
doSimpleGSEA <- function(gsea_fc_data,genesetCollection){
  library(data.table)
  library(piano)

  gsea_dt<-data.table()
  for (i_comp in seq_len(ncol(gsea_fc_data))) {

    fc <- gsea_fc_data[,i_comp]
    names(fc) <- rownames(gsea_fc_data)

    rgsa.mean <- runGSA(
      gsc = genesetCollection,
      geneLevelStats = fc,
      geneSetStat = "mean",
      nPerm=10000,
      ncpus = 1,
      verbose=FALSE
    )

    gsasumm <- GSAsummaryTable(gsaRes = rgsa.mean)

    keep_cols<-c(
      "Name",
      #"Genes (tot)", "Genes (up)", "Genes (down)",
      "p adj (dist.dir.up)", "p adj (dist.dir.dn)"#,
      #"p (dist.dir.up)","p (dist.dir.dn)"#,
    )

    gsasumm<-gsasumm[,keep_cols]

    # Get pvalue columns and Log them
    gsasumm[,grep("dist.dir",keep_cols,value=TRUE)] <- log(gsasumm[,grep("dist.dir",keep_cols,value=TRUE)])

    # Make dir.up pvalues absolute
    gsasumm[,grep("dir.up",keep_cols,value=TRUE)] <- abs(gsasumm[,grep("dir.up",keep_cols,value=TRUE)])

    gsasumm[,"exp"] <- colnames(gsea_fc_data)[i_comp]

    gsea_dt<-rbind(gsea_dt,gsasumm)
  }

  gsasumm_melt <- melt(
    gsea_dt, id.vars=c("Name","exp"),
    measure.vars=grep("dist.dir",keep_cols,value=TRUE)
  )

  gsasumm_melt$exp_cond <- ifelse(
    gsasumm_melt$variable %in% grep("dir.up",gsasumm_melt$variable,value=TRUE),
    paste0(gsasumm_melt$exp,"_UP"),
    paste0(gsasumm_melt$exp,"_DOWN")
  )

  gsasumm_melt$exp_cond <- as.factor(gsasumm_melt$exp_cond)

  return(gsasumm_melt)
}

doSimpleGSEA_wList <- function(list_gsea_fc,genesetCollection){
  # tmp_var <- cmlhc_limma_res$CML.Normal_FC
  # names(tmp_var)<-cmlhc_limma_res$Gene_Symbol
  # list_gsea_fc<-list("cmlhc"=tmp_var)
  # genesetCollection<-gsCollection
  # genesetCollection$gsc$IL2_STAT5_signaling<-NULL
  # str(list_gsea_fc)
  # summary(tmp_var)

  gsea_dt<-data.table()

  for (i_comp in seq_len(length(list_gsea_fc))) {

    fc <- list_gsea_fc[[i_comp]]

    #names(fc) <- rownames(list_gsea_fc)

    rgsa.mean <- runGSA(
      gsc = genesetCollection,
      geneLevelStats = fc,
      geneSetStat = "mean",
      nPerm=10000,
      ncpus = 1,
      verbose=FALSE
    )

    gsasumm <- GSAsummaryTable(gsaRes = rgsa.mean)

    keep_cols<-c(
      "Name",
      #"Genes (tot)", "Genes (up)", "Genes (down)",
      "p adj (dist.dir.up)", "p adj (dist.dir.dn)"#,
      #"p (dist.dir.up)","p (dist.dir.dn)"#,
    )

    gsasumm<-gsasumm[,keep_cols]

    # Get pvalue columns and Log them
    gsasumm[,grep("dist.dir",keep_cols,value=TRUE)] <- log(gsasumm[,grep("dist.dir",keep_cols,value=TRUE)])

    # Make dir.up pvalues absolute
    gsasumm[,grep("dir.up",keep_cols,value=TRUE)] <- abs(gsasumm[,grep("dir.up",keep_cols,value=TRUE)])


    gsasumm[,"exp"] <- names(list_gsea_fc)[i_comp]

    gsea_dt<-rbind(gsea_dt,gsasumm)
  }

  gsasumm_melt <- melt(
    gsea_dt, id.vars=c("Name","exp"),
    measure.vars=grep("dist.dir",keep_cols,value=TRUE)
  )

  gsasumm_melt$exp_cond <- ifelse(
    gsasumm_melt$variable %in% grep("dir.up",gsasumm_melt$variable,value=TRUE),
    paste0(gsasumm_melt$exp,"_UP"),
    paste0(gsasumm_melt$exp,"_DOWN")
  )

  gsasumm_melt$exp_cond <- as.factor(gsasumm_melt$exp_cond)

  return(gsasumm_melt)
}


consensus_GSA = function(
    gene_statistics,
    gene_set_collection
){
  # TODO implement
}




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
  dge <- DGEList(counts=expression_matrix)

  # filter lowly expressed genes if this is specified
  if(filter){
    print(">>> Filter out lowly expressed genes", quote=FALSE)
    keep <- filterByExpr(dge, design)
    print(paste(c(">>>> Number of genes filtered out: ", length(keep[keep == FALSE])), sep="", collapse=""), quote=FALSE)
    print(paste(c(">>>> Number of genes still considerd: ", length(keep[keep == TRUE])), sep="", collapse=""), quote=FALSE)
    dge <- dge[keep,,keep.lib.sizes=FALSE]
    # TODO PCA plot after filtering
  }

  # TODO add possibility to normalize (also plot normalized counts if this option is chosen)
  print(">>> Normalize Counts", quote=FALSE)
  dge <- calcNormFactors(dge, method=normalize)
  if(normalize!="none"){
    # TODO PCA plot after normalization
  }

  print(">>> Fit model", quote=FALSE)
  if(removeIntercept){
    # voom was shown to have the edge when the sequencing depths
    # were very inconsistent between replicates but, otherwise,
    # limma-trend was just as good
    if(limma_voom){
      vv = limma::voom(expression_matrix,design,plot=TRUE)
      # TODO plot tranformed counts
      fit <- limma::lmFit(vv,design)
      limma_trend=FALSE # if using limma-voom don't use limma-trend
    } else if (limma_voom_weight_samples){
      vv = limma::voomWithQualityWeights(expression_matrix,design,normalization="none",plot=TRUE)
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
    
    cairo_pdf(file=paste(c(output_dir, "Model_plots.pdf"), sep="", collapse=""), width=11.43, height=8.27, onefile=T)
    limma::plotSA(fit2, main="Final model: Mean-variance trend")
    dump <- dev.off()

  } else {
    fit <- limma::lmFit(expression_matrix,design)
    fit2 <- limma::eBayes(fit,trend=limma_trend)
  }

  print(">>> Evaluate results", quote=FALSE)
  limma_res <- list()
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
    
    limma_res[[contrast_name_fc]]<-aux[,"logFC"]
    limma_res[[contrast_name_pval]]<-aux[,"adj.P.Val"]

    # -1, 0 or 1 depending on whether each t-statistic is classified as
    # significantly negative, not significant or significantly positive, respectively
    limma_res[[contrast_name_de]] <- results@.Data[,i]

    # ID col only appears in the case of duplicated rownames in fit
    if("ID" %in% colnames(aux)){
      limma_res[["Gene_Symbol"]] = aux[,"ID"]
    } else {
      limma_res[["Gene_Symbol"]] = rownames(aux)
    }
    
    if(include_unadj_pval){
      contrast_name_unadjpval <- paste0(contrast_name,"_unadjPval")
      limma_res[[contrast_name_unadjpval]]<-aux[,"P.Value"]
    }
  }

  limma_res <- as.data.frame(limma_res)
  # limma_res$Gene_Symbol <- rownames(limma_res)

  return(list(
    topGenes=limma_res,
    fit=fit2
  ))
}



makeVolcanoPlot <- function(
  limma_res,
  name_tag,
  highlight_genes=NULL,
  geneset_name=NULL,
  plot_dir=NULL,
  hg_text=TRUE,
  plt_subtitle="Text highlights DE genes in geneset with FC>1",
  create_plot_subdir=TRUE,
  formats=c("png","pdf","tiff")
){
  # Makes volcano plot from results obtained from my limma_pipeline function
  # note that volcanoplot from limma uses normal pvalue while im using
  # adjusted.pvalues hence the difference when visualizing it.

  # TODO add functionality: option to automatically label highest fold change
  # and significant probes

  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ggplot2))

  fc_col <- grep("_FC",colnames(limma_res),value=TRUE)
  pval_col <- grep("_Pval",colnames(limma_res),value=TRUE)
  gs_col <- grep("Gene_Symbol",colnames(limma_res),value=TRUE)
  de_col <- grep("_DE",colnames(limma_res),value=TRUE)
  facet_col = NULL

  if(length(fc_col)>1){

    fc_col = "FC"
    pval_col = "Pval"
    de_col = "DE"
    facet_col = "VS"

    molten_pval = reshape2::melt(
    limma_res,id.vars = c("Gene_Symbol"),
    measure.vars = grep("_Pval$",colnames(limma_res)),
    value.name=pval_col
    )
    molten_fc = reshape2::melt(
      limma_res,id.vars = c("Gene_Symbol"),
      measure.vars = grep("_FC$",colnames(limma_res)),
      value.name=fc_col
    )
    molten_de = reshape2::melt(
      limma_res,id.vars = c("Gene_Symbol"),
      measure.vars = grep("_DE$",colnames(limma_res)),
      value.name=de_col
    )
    limma_res = dplyr::bind_cols(molten_pval,molten_fc,molten_de)
    limma_res[,facet_col] = gsub("_Pval$","",limma_res$variable)

  }

  limma_res[,"minusLog10pval"] <- -log10(limma_res[,pval_col])

  xlim_max <- max(abs(limma_res[,fc_col]))+0.8
  ylim_max <- ifelse(
    max(abs(limma_res[,"minusLog10pval"]))<(-log10(0.05)),
    (-log10(0.01)),
    max(abs(limma_res[,"minusLog10pval"]))
  )

  plt_title = paste0(fc_col," ",geneset_name)
  volc_plt <- ggplot(
    limma_res,
    aes_(
      x=as.name(fc_col),
      y=as.name("minusLog10pval"),
      label=as.name(gs_col)
    )
  )+geom_point(color="grey50")+
  geom_point(
    data=subset(limma_res,limma_res[,"Gene_Symbol"] %in% highlight_genes),
    color="red")+
  scale_color_brewer(palette="Set1")+
  geom_vline(xintercept=1,linetype="dashed", color = "grey80")+
  geom_vline(xintercept=-1,linetype="dashed", color = "grey80")+
  geom_hline(yintercept=-log10(0.05),linetype="dashed", color = "grey80")
  if(hg_text){
    volc_plt <-volc_plt + geom_text_repel(
      # genes in HL gene list that are differentially expressed with a FC > 1
      data = subset(
        limma_res,
        limma_res[,"Gene_Symbol"] %in% highlight_genes & limma_res[,de_col]!=0 & abs(limma_res[,fc_col])>1
      ),
      segment.size = 0.5,
      segment.color = "darkred",
      direction = "both",
      size = 8,
      box.padding = 2,
      ylim=c(1.5,NA),
      max.iter=5000,force=1.5
    )
    if(is.null(highlight_genes)){
      # if highlight text but no highlighted genes, then highlight DE + abs(FC)>1 genes
      volc_plt <- volc_plt +
        geom_point(
          data = subset(limma_res, limma_res[,de_col]!=0 & abs(limma_res[,fc_col])>1),
          color="red"
        )+
        geom_text_repel(
          data = subset(limma_res, limma_res[,de_col]!=0 & abs(limma_res[,fc_col])>1),
          segment.size = 0.5,
          segment.color = "darkred",
          direction = "both",
          size = max(2,(8 * (1/length(which(limma_res[,de_col]!=0 & abs(limma_res[,fc_col])>1))))),
          box.padding = .5,
          max.iter=5000,
          force=1.
        )
    }
  }
  volc_plt = volc_plt +
  geom_label_repel(
    data=data.frame(lab="P-value = 0.05",pval=-log10(0.05),xpos=xlim_max),
    mapping=aes(label=lab,x=xpos,y=pval),
    color = "grey60",
    alpha=0.8
  )+
  labs(
    title=plt_title,
    subtitle=plt_subtitle,
    y="-log10(P-value)",
    x="log2 Fold Change"
  )+
  coord_cartesian(xlim=c(-xlim_max,xlim_max),ylim=c(0,ylim_max))+
  theme_bw(base_size=16)+
  scale_x_continuous(
    minor_breaks=NULL,
    breaks=c(-ceiling(xlim_max):ceiling(xlim_max)))+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(),
    panel.border = element_blank()
  )

  if(!is.null(facet_col)) volc_plt=volc_plt+facet_grid(VS~.)

  save_different_plot_format(
    plt=volc_plt,plot_dir=plot_dir,
    create_plot_subdir=create_plot_subdir,save_device="ggplot",
    type_name="volcano_plot",name_tag=name_tag,formats=formats,
    units="cm",width=20,height=20
  )

  return(volc_plt)
}

save_different_plot_format = function(
  plt=NULL,plot_dir=NULL,create_plot_subdir=FALSE,
  save_device=c("ggplot","grDevice","complexheatmap"),
  type_name="",name_tag="",formats=c("png","pdf","tiff"),
  units="cm",width=20,height=20,...
){
  if(!is.null(plot_dir) & !is.null(plt)){
    save_device = match.arg(save_device)

    f_name = paste0(type_name,"-",name_tag)

    for(fmt in formats){
      f_path_fmt = file.path(plot_dir,paste0(f_name,".",fmt))
      if(create_plot_subdir)  dir.create(file.path(plot_dir,fmt),recursive=TRUE,showWarnings=FALSE)
      if(dir.exists(file.path(plot_dir,fmt))) f_path_fmt = file.path(plot_dir,fmt,paste0(f_name,".",fmt))
      if(save_device=="ggplot"){
        ggplot2::ggsave(filename=f_path_fmt,plot = plt,device = fmt,units = units,width = width,height = height,...)
      }
      if(save_device=="complexheatmap"){
        save_f = get(fmt)
        if(fmt=="png"|fmt=="tiff")  save_f(f_path_fmt,width=width,height=height,units="in",res=300,...)
        else  save_f(f_path_fmt,width=width,height=height)
        ComplexHeatmap::draw(plt)
        dev.off()
      }
    }
  }
}




makePCA <- function(
  data_matrix, # matrix to run the pca on
  label_factor, # label for each sample
  color_palette="lancet", # color palette to be used
  addEllipses=FALSE,
  ellipse_level=0.95,
  title="PCA",
  subtitle=NULL,
  colour_label=waiver(),
  name_tag=NULL,
  plot_dir=NULL,...
){
  # Samples as rows
  # Features as columns
  suppressPackageStartupMessages(library(ggplot2))

  assertthat::assert_that(
    class(label_factor) %in% c("factor","numeric")
  )
  assertthat::assert_that(
    length(label_factor) == nrow(data_matrix)
  )

  pca_res<-FactoMineR::PCA(data_matrix,scale.unit=FALSE,graph=FALSE)

  if(is.factor(label_factor)){
    pca_plt<-factoextra::fviz_pca_ind(
      pca_res,
      palette=color_palette,
      habillage=label_factor,
      pointshape=19,
      mean.point=FALSE,
      addEllipses=addEllipses,
      ellipse.level=ellipse_level
    )
  }else{
    pca_plt<-factoextra::fviz_pca_ind(
      pca_res,
      palette=color_palette,
      col.ind=label_factor,
      pointshape=19,
      mean.point=FALSE,
      addEllipses=addEllipses,
      ellipse.level=ellipse_level
    )
  }

  pca_plt = pca_plt+
    ggplot2::labs(title=title,subtitle=subtitle,colour=colour_label)+
    theme_minimal(base_size=16)

  save_different_plot_format(
    plt=pca_plt,plot_dir=plot_dir,
    create_plot_subdir=TRUE,save_device="ggplot",
    type_name="pca",name_tag=name_tag,
    units="cm",...
  )
  return(pca_plt)
}


pltFCHeatmapWrapper <- function(
  topGenesTbl,
  #genesetTbl,
  name_tag,
  genes_of_int,
  #peak_height_filter=50,
  column_names=NULL,
  plot_dir=NULL,
  create_plot_subdir=TRUE,
  row_fontsize=5,
  plot_width=12,
  plot_height=20
)
{
  suppressPackageStartupMessages(library(ComplexHeatmap))

  hm_data<-topGenesTbl[,grep("_FC",colnames(topGenesTbl)),drop=F]
  #hm_data<-hm_data[intersect(rownames(hm_data),genesetTbl[genesetTbl$"Nearest gene"%in%genes_of_int,]$"Nearest gene"),,drop=F]


  # hmcol <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(256)
  maxfcscaleval<-ceiling(max(abs(hm_data)))
  red_shade<-RColorBrewer::brewer.pal(9,"RdBu")[1]
  blue_shade<-RColorBrewer::brewer.pal(9,"RdBu")[9]
  hmcol<-circlize::colorRamp2(breaks=c(-maxfcscaleval, 0, maxfcscaleval), colors=c(blue_shade, "#FFFFFF", red_shade))

  #tbp<-genesetTbl[grep(paste0(rownames(hm_data),collapse="|"),genesetTbl$"Nearest gene"),]
  #tbp<-tbp[!duplicated(tbp$"Nearest gene"),]
  #tbp<-as.data.frame(tbp)
  #row.names(tbp) <- tbp$"Nearest gene"
  #tbp<-tbp[rownames(hm_data),,drop=F]

  # row_ha <- HeatmapAnnotation("Peak Height"=anno_barplot(tbp$Peak,which="row"),which = "row")

  hmplt<-Heatmap(
    matrix=as.matrix(hm_data),
    name="FC",
    column_title=paste0(name_tag,"\nFold change gene expression \nGene of interest in red"),
    col=hmcol,
    #right_annotation=row_ha,
    column_dend_reorder=FALSE,
    row_dend_reorder=FALSE,
    row_dend_width=unit(2, "cm"),
    column_names_rot=90,
    column_labels=if(!is.null(column_names)){column_names} else {colnames(hm_data)},
    row_names_gp=gpar(col = ifelse(rownames(hm_data)%in%genes_of_int,"red","black"),fontsize = row_fontsize),
    heatmap_width=unit(plot_width, "cm")
  )

  save_different_plot_format(
    plt=hmplt,plot_dir=plot_dir,
    create_plot_subdir=create_plot_subdir,save_device="complexheatmap",
    type_name="fc_heatmap",name_tag=name_tag,
    units="cm",width=plot_width,height=plot_height
  )

  return(hmplt)

}


run_consensus_GSA = function(
  fc_gene_mat,
  gene_sets_fpath,
  method_stat=c(
    "mean","median","sum",
    "gsea","fgsea","maxmean",
    "wilcoxon","page"
  ),
  n_perm=100,
  n_cpus=1,
  output_dir
){
  # TODO add enrichment score plot for top significant sets
  #fgsea::plotEnrichment(
  #  lgsc$gsc$DAZARD_RESPONSE_TO_UV_SCC_DN,
  #  fc_vec
  #)

  # library("piano")
  library("snow")
  library("snowfall")

  tictoc::tic("Consensus GSA run")

  assertthat::assert_that(is.matrix(fc_gene_mat))

  lgsc = piano::loadGSC(gene_sets_fpath)

  res = lapply(colnames(fc_gene_mat),function(col_n){
    tictoc::tic(paste0(basename(gene_sets_fpath)," ",col_n))
    fc_vec=fc_gene_mat[,col_n]
    names(fc_vec) = rownames(fc_gene_mat)
    rgsa_mean=NULL
    rgsa_median=NULL
    rgsa_sum=NULL
    rgsa_gsea=NULL
    rgsa_fgsea=NULL
    rgsa_maxmean=NULL
    rgsa_wilcoxon=NULL
    rgsa_page=NULL
    if("mean"%in%method_stat){
      rgsa_mean <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "mean",
        adjMethod="BH",gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }
    if("median"%in%method_stat){
      rgsa_median <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "median",
        adjMethod="BH",gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }
    if("sum"%in%method_stat){
      rgsa_sum <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "sum",
        adjMethod="BH",gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }
    if("gsea"%in%method_stat){
      rgsa_gsea <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "gsea",
        adjMethod="BH",# FIXME only accepts fdr but fails in some cases when this is set, check source
        gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }
    if("fgsea"%in%method_stat){
      rgsa_fgsea <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "fgsea",signifMethod="geneSampling",
        adjMethod="BH",gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }
    if("maxmean"%in%method_stat){
      rgsa_maxmean <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "maxmean",
        adjMethod="BH",gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }
    if("wilcoxon"%in%method_stat){
      rgsa_wilcoxon <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "wilcoxon",
        adjMethod="BH",gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }
    if("page"%in%method_stat){
      rgsa_page <- piano::runGSA(geneLevelStats = fc_vec,
        geneSetStat = "page",
        adjMethod="BH",gsc = lgsc,nPerm=n_perm,ncpus = n_cpus,verbose=FALSE)
    }

    res_list = list(
      rgsa_mean,
      rgsa_median,
      rgsa_sum,
      rgsa_gsea,
      rgsa_fgsea,
      rgsa_maxmean,
      rgsa_wilcoxon,
      rgsa_page
    )

    res_list = res_list[!sapply(res_list,is.null)]

    res_mat = piano::consensusHeatmap(
      resList=res_list,
      cutoff=Inf, # ensures every geneset is present
      # only applies to the rank and not other
      # stats (pval/nGenes) which are always median
      method="mean",
      adjusted=TRUE,
      plot=FALSE
    )

    f_res = list(
      "consensus_res" = res_mat,
      "gsea_res" = res_list
    )

    saveRDS(
      object = f_res,
      file = file.path(
        output_dir,
        paste0("GSEA_",basename(gene_sets_fpath),"_results_",col_n,".RDS")
      ),
      compress = FALSE
    )

    tictoc::toc()
    return(f_res)
  })
  names(res) = colnames(fc_gene_mat)

  saveRDS(
    object = res,
    file = file.path(
      output_dir,
      paste0("GSEA_",basename(gene_sets_fpath),"_results_","all_comparisons_list",".RDS")
    ),
    compress = FALSE
  )

  tictoc::toc()
  return(res)
}

# Experimental function
plot_consensus_heatmap_res = function(
  consensus_res,
  name_tag,
  n_top = 10,
  plot_dir=NULL,
  create_plot_subdir=TRUE,
  plot_width=15,
  plot_height=10
){
  library("ComplexHeatmap")

  consensus_pval_mat = as.data.frame(consensus_res$pMat[,c(1,5)])
  colnames(consensus_pval_mat) = c("distinct_dirUP","distinct_dirDOWN")
  top_consensus_pval_df = rbind(
    head(consensus_pval_mat[order(consensus_pval_mat$distinct_dirUP),],n_top),
    tail(consensus_pval_mat[order(consensus_pval_mat$distinct_dirDOWN,decreasing=TRUE),],n_top)
  )
  #summary(consensus_pval_mat)

  consensus_rankscore_mat = as.data.frame(consensus_res$rankMat[,c(1,5)])
  colnames(consensus_rankscore_mat) = c("distinct_dirUP","distinct_dirDOWN")
  top_consensus_rankscore_df = rbind(
    head(consensus_rankscore_mat[order(consensus_rankscore_mat$distinct_dirUP),],n_top),
    tail(consensus_rankscore_mat[order(consensus_rankscore_mat$distinct_dirDOWN,decreasing=TRUE),],n_top)
  )
  #summary(consensus_rankscore_mat)

  table(rownames(top_consensus_rankscore_df)%in%rownames(top_consensus_pval_df))

  max_pval_leg = ifelse(
    max(-log10(as.matrix(top_consensus_pval_df)))<(-log10(0.05)),
    2,
    ceiling(max(-log10(as.matrix(top_consensus_pval_df))))
  )

  col_mid_rank = max(
    max(head(sort(consensus_rankscore_mat$distinct_dirUP),n_top)),
    max(head(sort(consensus_rankscore_mat$distinct_dirDOWN),n_top))
  )
  max_rank = ceiling(max(as.matrix(top_consensus_rankscore_df)))
  col_rank = circlize::colorRamp2(
    c(1,col_mid_rank,max_rank),
    c("red","white","green")
  )

  hm_cs_pval_plt = ComplexHeatmap::Heatmap(
    -log10(as.matrix(top_consensus_pval_df)),
    col=circlize::colorRamp2(
      c(0, -log10(0.05), max_pval_leg),
      c("blue", "white", "red")
    ),
    heatmap_legend_param=list(direction="horizontal"),
    name="-log10(adjusted pval) (higher is better)",
    column_title = paste0(name_tag,"-",col_n,"\nConsensus median adjusted p-values and rank of gene-sets\n(significantly",
    " affected) by gene regulation\nTop ",n_top," up and down regulated"),
    heatmap_width=grid::unit(plot_width, "cm"),
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_rot=45,
    cluster_columns=FALSE,
    cluster_rows=FALSE
  )
  hm_cs_rank_plt = ComplexHeatmap::Heatmap(
    as.matrix(top_consensus_rankscore_df),
    col=col_rank,
    heatmap_legend_param=list(
      direction="horizontal",
      at=c(1,max_rank)
    ),
    name="Rank (lower is better)",
    heatmap_width=grid::unit(plot_width, "cm"),
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_rot=45,
    cluster_columns=FALSE,
    cluster_rows=FALSE
  )
  hm_cs_list = hm_cs_pval_plt%v%hm_cs_rank_plt
  hmplt = ComplexHeatmap::draw(hm_cs_list, heatmap_legend_side = "bottom")

  save_different_plot_format(
    plt=hmplt,plot_dir=plot_dir,
    formats=c("png"),
    create_plot_subdir=create_plot_subdir,save_device="complexheatmap",
    type_name="cs_heatmap",name_tag=paste0(name_tag,"-",col_n),
    units="cm",width=plot_width,height=plot_height
  )
}

plot_limma_topGenes_tbl = function(
  limma_topTable,
  tbl_name,
  comparison,
  output_dir,
  width_px="100%",
  sort_by=c("Pval","FC"),
  n_top=10
){

  suppressPackageStartupMessages(library(formattable))
  #webshot::install_phantomjs()

  export_formattable <- function(
    f_table,
    file,
    width = "100%",
    height = NULL,
    background = "white",
    delay = 0.2){
      # from https://github.com/renkun-ken/formattable/issues/26
      w <- formattable::as.htmlwidget(f_table, width = width, height = height)
      path <- htmltools::html_print(w, background = background, viewer = NULL)
      url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
      webshot::webshot(
        url,
        file = file,
        selector = ".formattable_widget",
        delay = delay
      )
  }

  customGreen = "#71CA97"
  customRed = "#ff7f7f"

  # Select desired data
  limma_topTable = limma_topTable %>%
    select(contains(comparison),"Gene_Symbol") %>%
    arrange_at(vars(ends_with(sort_by)))%>%head(n_top)

  # Columns in proper order
  limma_topTable = limma_topTable %>%
    select("Gene_Symbol",contains("_DE"),contains("_FC"),contains("_Pval"))

  # Proper column names for table
  colnames(limma_topTable) = c(
    "Gene Symbol",
    "DE direction",
    "Log2 Fold-Change",
    "Adjusted P-value"
  )

  format_list = list(
    "Gene Symbol" = formattable::formatter("span", style = ~ formattable::style(color = "grey",font.weight = "bold")),
    "DE direction"=formattable::formatter("span", style = x ~ style(color = ifelse(x > 0, "red", ifelse(x < 0, "blue", "black")))),
    "Log2 Fold-Change"=formattable::color_bar("lightgrey",fun=function(x){formattable::proportion(abs(x))}),
    "Adjusted P-value"= formattable::color_tile(customGreen, customRed)
  )

  # theres a bug with this viz where the minus appears on the right for negative numbers
  # making it abs to avoid this, direction can be seen by DE dir column anyway
  limma_topTable[,"Log2 Fold-Change"] = format(
    abs(limma_topTable[,"Log2 Fold-Change"]),
    digits=2,
    nsmall=2
  )

  ft_tbl = formattable::formattable(
    limma_topTable,
    align =c("l","r","r","r"),
    format_list
  )
  rownames(ft_tbl) = NULL

  export_formattable(
    f_table=ft_tbl,
    file=file.path(
      output_dir,
      paste0("limma_topGenes_table-",comparison,"-",tbl_name,".png")
    ),
    width=width_px
  )
}

