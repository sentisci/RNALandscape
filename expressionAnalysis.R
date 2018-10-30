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
  projectName             = "RNASeq.RSEM",
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

## parse ssGSEA output file and zscore
expressionTMM.RPKM.ssGSEA.output <- corUtilsFuncs$parseBroadGTCOutFile( fileName = paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                                            "/",rnaseqProject$gseaDir,"/results/","expressionTMM.RPKM.ssGSEA.Output.gct"))

## append Cytolytic score
cytolyticScore <- corUtilsFuncs$cytolyticScore(expDF = expressionTMM.RPKM.PC)
expressionTMM.RPKM.ssGSEA.output.Cytolytic <- rbind(expressionTMM.RPKM.ssGSEA.output, cytolyticScore)
expressionTMM.RPKM.ssGSEA.output.Cytolytic.zscore <- t( apply(expressionTMM.RPKM.ssGSEA.output.Cytolytic, 1, corUtilsFuncs$zscore_All) ) %>% data.frame()
expressionTMM.RPKM.ssGSEA.output.Cytolytic.zscore <- expressionTMM.RPKM.ssGSEA.output.Cytolytic.zscore %>% tibble::rownames_to_column(var="GeneSignatures")
write.table(expressionTMM.RPKM.ssGSEA.output.Cytolytic.zscore, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,"/",rnaseqProject$gseaDir,
                      "/results/","expressionTMM.RPKM.ssGSEA.Output.Cytolytic.zscore.gct"),sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Perform Differential gene expression analysis

## Control groups ##
Brain          <- c("NS.cerebellum","NS.cerebrum")
Heart          <- c("NS.heart", "NS.heart")
Kidney         <- c("NS.kidney")
Liver          <- c("NS.liver")
Lung           <- c("NS.lung")
germline       <- c("NS.testis","NS.ovary")
vitalNormals   <- c("NS.heart","NS.kidney","NS.liver","NS.lung")
vital.Brain.Normals   <- c("NS.cerebellum","NS.cerebrum", "NS.heart","NS.kidney","NS.liver","NS.lung")
othersNormals  <- c("NS.adrenalgland","NS.bladder","NS.colon","NS.ileum","NS.ovary","NS.pancreas","NS.prostate", 
                    "NS.skeletalmuscle","NS.spleen", "NS.stomach","NS.testis", "NS.ureter", "NS.uterus")
Normals        <- c("NS.adrenalgland","NS.bladder","NS.cerebellum","NS.cerebrum","NS.colon","NS.heart",
                    "NS.ileum","NS.kidney","NS.liver","NS.lung","NS.ovary","NS.pancreas","NS.prostate", 
                    "NS.skeletalmuscle","NS.spleen", "NS.stomach","NS.testis", "NS.ureter", "NS.uterus")
NormalsNoGermLine <- c("NS.adrenalgland","NS.bladder","NS.cerebellum","NS.cerebrum","NS.colon","NS.heart",
                       "NS.ileum","NS.kidney","NS.liver","NS.lung","NS.pancreas","NS.prostate", 
                       "NS.skeletalmuscle","NS.spleen", "NS.stomach", "NS.ureter", "NS.uterus")

tumorSubStatus.polyA <- c("ASPS", "DSRCT", "EWS" ,"HBL", "ML", "NB.MYCN.NA","NB.MYCN.A", "NB.Unknown", "OS", "RMS.FP" , "RMS.FN", 
                          "SS", "Teratoma" ,"UDS" ,"YST")
tumorSubStatus.ribozero <-  c("WT" ,"CCSK")

Tumors         <-  c("ASPS","DSRCT", "EWS" ,"HBL", "ML", "NB" ,"OS", "RMS", "SS", "Teratoma" ,"UDS" ,"YST","WT", "CCSK")


## Testing 
dgeObj  <- DifferentialGeneExp$new(
  countObj          = expressionObj$edgeRMethod("NormFactorDF")$counts,
  group1            = list(list("NormalsNoGermLine"=NormalsNoGermLine,each=FALSE)),
  group2            = list(list("Tumor"=tumorSubStatus.polyA, each=TRUE)),
  packageRNAseq     = "edgeR",
  groupColumnName   = rnaseqProject$factorName,
  metadataDF        = rnaseqProject$metaDataDF,
  samplesColumnName = "SAMPLE_ID",
  expressionUnit    = "TMM-RPKM",
  featureType       = "Gene",
  writeFiles        = TRUE
)
DiffExpObj <- dgeObj$performDiffGeneExp()




