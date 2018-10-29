rm(list=ls())
## Source all classes and packages

source("./utilityPackages.R")
source("./statisticalPackages.R")
source("./class.R")

## Project Title: Expression Analysis for Landscape paper

## ## Instantiate a new Object of type ProjectSetUp
rnaseqProject <- ProjectSetUp$new(
  
  date                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[1]],
  time                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[2]],
  projectName             = "Nitya.CART.Human",
  annotationRDS           = "C:/Users/sindiris/R Scribble/Annotation RDS/annotation_06302016_ensembl.CAR_gene.RDS",
  pcRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/pc.other.HGNCTableFlat.rds",
  tfRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/TFs_no_epimachines.RDS",
  csRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/CellSurface.RDS",
  cgaRDS                  = "C:/Users/sindiris/R Scribble/Annotation RDS/cancerGermlineAntigens.rds",
  outputPrefix            = "landscape",
  filterGenes             = TRUE,
  filterGeneMethod        = "bySum",
  factorName              = "DIAGNOSIS.Substatus",
  metaDataFileName        = "MetadataMapper.txt",
  outputdirRDSDir         = "GeneRDSOutput",
  outputdirTXTDir         = "GeneTXTOutput",
  gseaDir                 = "GSEA",
  plotsDir                = "Figures",
  plotsDataDir            = "FigureData",
  DiffGeneExpAnaDir       = "DiffExpResults",
  DiffGeneExpRDS          = "DiffGeneExpRDSOutput",
  #factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", 
  #                                                                                                  "LibraryPrep" = "Ribozero"))
)

## Add utility functions to the project
corUtilsFuncs <- CoreUtilities$new(  ProjectSetUpObject = rnaseqProject )

## Generate expression matrix
rm(mergeObjectsNoDup)
mergeObjectsNoDup <- corUtilsFuncs$getMergedMatrix(dir               = "TPM_Genes", 
                                                   fileFormat        = "txt", 
                                                   colNameSelect     = "expected_count", 
                                                   isRowNames        = TRUE, 
                                                   rowNamesColInFile = 1,
                                                   fileSuffix        = ".genes.results",
                                                   primaryID         = "gene_id",
                                                   metadata          = rnaseqProject$metaDataDF,
                                                   metadataFileRefCol= "SAMPLE_ID")

## Evaluate presence of duplicate features and consolidate them
setDT(mergeObjectsNoDup, keep.rownames = TRUE)
mergeObjectsConso <- corUtilsFuncs$consolidateDF(mergeObjectsNoDup, funcName = "sum", featureName = "rn")
mergeObjectsConso <- mergeObjectsConso %>% data.frame() %>% tibble::column_to_rownames(var = "rn") %>% as.matrix()

## Subset metaDataDF by the number of samples in the folder
colnamesDF    <- data.frame( "SAMPLE_ID"= colnames(mergeObjectsNoDup))
corUtilsFuncs$subsetMetaData(colnamesDF=colnamesDF)

## Instantiate a new Object of type GeneExpNormalization
expressionObj <- GeneExpNormalization$new(
  
  countObj       = as.matrix(mergeObjectsConso), 
  featureType    = "Gene", 
  packageRNAseq  = "edgeR", 
  annotationDF   = rnaseqProject$annotationDF, 
  design         = rnaseqProject$metaDataDF[,rnaseqProject$factorName], 
  #design         = newMetaDataDF[,rnaseqProject$factorName],
  proteinCodingOnly = FALSE
)

## Get expression in desired units
expressionTMM.RPKM = expressionObj$edgeRMethod("TMM-RPKM")
write.table(expressionTMM.RPKM, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,"/",rnaseqProject$outputdirTXTDir,"/","expressionTMM.RPKM.txt"),
            sep="\t", quote = FALSE, row.names = FALSE)

expressionRawCounts = expressionObj$edgeRMethod("RawCounts")
write.table(expressionRawCounts, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,"/",rnaseqProject$outputdirTXTDir,"/","rawCounts.txt"),
            sep="\t", quote = FALSE, row.names = FALSE)

## Perfrom unsupervised clustering
expressionTMM.RPKM.ForHC  <- expressionTMM.RPKM %>% tibble::column_to_rownames(var="GeneID")
expressionTMM.RPKM.ForHC.zscore <- t(apply( expressionTMM.RPKM.ForHC[,-c(1:6)], 1, corUtilsFuncs$zscore_All))
corUtilsFuncs$performClustering(expressionTMM.RPKM.ForHC.zscore)

## generate ssGSEA input file
pcList <- corUtilsFuncs$selectGenes(geneMatrix=expressionTMM.RPKM, selectGeneDF=rnaseqProject$pcDF, 
                                                                annotationDF=rnaseqProject$annotationDF, featureType=c("GeneName","GeneID") ) 
expressionTMM.RPKM.PC     <- pcList[[1]] %>% tibble::column_to_rownames(var= "GeneName") %>% dplyr::select(-c(1:9))
expressionTMM.RPKM.ssGSEA <- corUtilsFuncs$createBroadGCTFile(expressionTMM.RPKM.PC)
write.table(expressionTMM.RPKM.ssGSEA, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,"/",rnaseqProject$gseaDir,"/rnk/","expressionTMM.RPKM.ssGSEA.Input.gct"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## Perform Differential gene expression analysis

## Control groups ##
IM_tumor_IM_UTD       <- c("IM_tumor_IM_UTD")
IM_tumor_IV_UTD       <- c("IM_tumor_IV_UTD")
IV_tumor_lung_met_UTD <- c("IV_tumor_lung_met_UTD")
IM_tumor              <- c("IM_tumor") 


##Condition group
IM_tumor_IM_CART        <- c("IM_tumor_IM_CART")
IM_tumor_IV_CART        <- c("IM_tumor_IV_CART")
IV_tumor_lung_met_CARs  <- c("IV_tumor_lung_met_CART_EF1a", "IV_tumor_lung_met_CART_MSCV")
IV_tumor_lung_met       <- c("IV_tumor_lung_met")


## Perform Differential gene expression analysis
dgeObj  <- DifferentialGeneExp$new (
  countObj          = expressionObj$edgeRMethod("NormFactorDF")$counts,
  group1            = list(list("IM_tumor_IM_UTD"=IM_tumor_IM_UTD,each=FALSE), list("IM_tumor_IV_UTD"=IM_tumor_IV_UTD,each=FALSE),
                           list("IV_tumor_lung_met_UTD"=IV_tumor_lung_met_UTD, each=FALSE), list("IM_tumor"=IM_tumor, each=FALSE) ),
  group2            = list(list("IM_tumor_IM_CART"=IM_tumor_IM_CART, each=FALSE), list("IM_tumor_IV_CART"=IM_tumor_IV_CART, each=FALSE),
                           list("IV_tumor_lung_met_CARs"=IV_tumor_lung_met_CARs, each=TRUE), list("IV_tumor_lung_met"=IV_tumor_lung_met, each=FALSE)),
  OneToOne          = TRUE,
  packageRNAseq     = "edgeR",
  groupColumnName   = rnaseqProject$factorName,
  metadataDF        = rnaseqProject$metaDataDF,
  samplesColumnName = "SAMPLE_ID",
  expressionUnit    = "TMM-RPKM",
  featureType       = "Gene",
  subsetGenes       = TRUE,
  writeFiles        = TRUE
)

DiffExpObj <- dgeObj$performDiffGeneExp()




