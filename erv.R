each_group_samples <- Tumor_samples_annot %>% dplyr::filter_(.dots =paste0("grepl(\"","ASPS","\",","DIAGNOSIS.Substatus.Tumor.Normal.Tissue",")") )
each_group_samplesDF <- erv.data.max.cytolytic.Tumor %>% dplyr::select(one_of(each_group_samples[,"Sample.Data.ID"]))


## Flatten rcorr results in to a table
flattenCorrMatrix <- function(cormat, pmat, group_name="") {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    group = group_name
  )
}

## Do Correlation for each diagnosis groups
rcorr_groups = function(annotDF = NULL, dataDF = NULL, GroupID =NULL, samplesID= NULL) {
  unique_groups <- unique(as.character(annotDF[,GroupID]))
  mergedDF <- do.call(rbind, lapply(unique_groups, function(x) {
    each_group_samples <- annotDF %>% dplyr::filter_(.dots =paste0("grepl(\"",x,"\",",GroupID,")") )
    each_group_samplesDF <- dataDF %>% dplyr::select(one_of(each_group_samples[,samplesID]))
    print(paste("group name", x,"  ", dim(each_group_samplesDF)[1], " ",dim(each_group_samplesDF)[2] ))
    each_group_samplesDF.corr <- rcorr(t(each_group_samplesDF), type = "spearman");
    CorrDF <- flattenCorrMatrix(each_group_samplesDF.corr$r, each_group_samplesDF.corr$P, x)
    median_exp <- apply(each_group_samplesDF,1,median)
    exp_df <- data.frame(ERV_median_exp=rep(median_exp[1:483],241), cytolyticScore=rep(median_exp[483],241))
    finaDF <- cbind(CorrDF,exp_df)
    return(finaDF)
  } ))
  return(mergedDF)
}

## Read ERV expression file

## Need to improve the below code
## Log the expression and saved
erv.data.max <- read.csv("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/data.geneNames.merged.max.v2.txt", sep="\t", header = T)
colnames(erv.data.max) <- gsub(".rsem.ENS.FPKM","",colnames(erv.data.max))
erv.data.max <- erv.data.max %>% tibble::column_to_rownames(var="geneName")
erv.data.max <- log(erv.data.max+1,2)
write.table(erv.data.max, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/data.geneNames.merged.max.v2.log2.txt", sep="\t", quote = F)
## Get cytolytic score and save
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)
write.table(cytolyticScore, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/cytolyticScore.txt", sep="\t", quote = F)

## Read the final expression file
erv.data.max.cytolytic <- read.csv("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/data.geneNames.merged.max.cytolytic.v2.log2.txt",
                                   sep="\t", header = T)
# erv.data.max.cytolytic.complete <- erv.data.max.cytolytic[complete.cases(erv.data.max.cytolytic),] %>% data.frame()
erv.data.max.cytolytic <- tibble::column_to_rownames(erv.data.max.cytolytic, var = "GeneID")

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Data.ID"=colnames(erv.data.max.cytolytic)), 
                                                   rnaseqProject$validMetaDataDF[,c("Sample.Biowulf.ID.GeneExp", "LibraryPrep",
                                                                                    "Sample.ID.Alias", "Sample.Data.ID", "DIAGNOSIS.Alias" ,rnaseqProject$factorName)] )
## Remove Normal samples
Tumor_samples_annot <- AliasNames_df %>% dplyr::filter(!grepl('Teratoma',DIAGNOSIS.Substatus.Tumor.Normal.Tissue))
Tumor_samples_annot <- Tumor_samples_annot[complete.cases(Tumor_samples_annot),]
erv.data.max.cytolytic.Tumor <- erv.data.max.cytolytic %>% dplyr::select(one_of(Tumor_samples_annot$Sample.Data.ID)); dim(erv.data.max.cytolytic.Tumor)

## Perform rcorr for each tumor groups
Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue <- gsub('NS.*','NS',Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue) 
correlationDF <- rcorr_groups(Tumor_samples_annot, erv.data.max.cytolytic.Tumor, "DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Sample.Data.ID")
correlationDF_cytolytic <- correlationDF %>% filter(grepl('CytolyticScore', column))

## Check ERVs highly correlated in normal
correlationDF_cytolytic.NS <- correlationDF_cytolytic %>% dplyr::filter(grepl('NS', group))
correlationDF_cytolytic.NS.ervs <- correlationDF_cytolytic.NS %>% filter(p <= 0.05 & cor > 0) %>% dplyr::select(row)
correlationDF_cytolytic <- correlationDF_cytolytic %>% dplyr::mutate(coorelated_with_Normal= ifelse(row %in% correlationDF_cytolytic.NS.ervs$row, "Y", "N"))
write.table(correlationDF_cytolytic, 
            "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/erv.data.max.cytolytic.Tumor.each.corr.withNormal.annot.medianExp.v2.txt", quote = FALSE, sep="\t")



################## Working version 1 ##########################
## Flatten rcorr results in to a table
flattenCorrMatrix <- function(cormat, pmat, group_name="") {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    group = group_name
  )
}

## Do Correlation for each diagnosis groups
rcorr_groups = function(annotDF = NULL, dataDF = NULL, GroupID =NULL, samplesID= NULL) {
  unique_groups <- unique(as.character(annotDF[,GroupID]))
  mergedDF <- do.call(rbind, lapply(unique_groups, function(x) {
    each_group_samples <- annotDF %>% dplyr::filter_(.dots =paste0("grepl(\"",x,"\",",GroupID,")") )
    each_group_samplesDF <- dataDF %>% dplyr::select(one_of(each_group_samples[,samplesID]))
    print(paste("group name", x,"  ", dim(each_group_samplesDF)[1], " ",dim(each_group_samplesDF)[2] ))
    each_group_samplesDF.corr <- rcorr(t(each_group_samplesDF), type = "spearman");
    return(flattenCorrMatrix(each_group_samplesDF.corr$r, each_group_samplesDF.corr$P, x))
  } ))
  return(mergedDF)
}
###############################################################

################## Working version 2 ##########################
## Flatten rcorr results in to a table
flattenCorrMatrix <- function(cormat, pmat, group_name="") {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
    group = group_name
  )
}

## Do Correlation for each diagnosis groups
rcorr_groups = function(annotDF = NULL, dataDF = NULL, GroupID =NULL, samplesID= NULL) {
  unique_groups <- unique(as.character(annotDF[,GroupID]))
  mergedDF <- do.call(rbind, lapply(unique_groups, function(x) {
    each_group_samples <- annotDF %>% dplyr::filter_(.dots =paste0("grepl(\"",x,"\",",GroupID,")") )
    each_group_samplesDF <- dataDF %>% dplyr::select(one_of(each_group_samples[,samplesID]))
    print(paste("group name", x,"  ", dim(each_group_samplesDF)[1], " ",dim(each_group_samplesDF)[2] ))
    each_group_samplesDF.corr <- rcorr(t(each_group_samplesDF), type = "spearman");
    CorrDF <- flattenCorrMatrix(each_group_samplesDF.corr$r, each_group_samplesDF.corr$P, x)
    median_exp <- apply(each_group_samplesDF,1,median)
    exp_df <- data.frame(ERV_median_exp=rep(median_exp[1:483],241), cytolyticScore=rep(median_exp[483],241))
    finaDF <- cbind(CorrDF,exp_df)
    return(finaDF)
  } ))
  return(mergedDF)
}
###############################################################










