rm(list=ls())

## Source all classes and packages ####

source("./utilityPackages.R")
source("./statisticalPackages.R")
source("./class.R")

## Project Title: Expression Analysis for Landscape paper

## Instantiate a new Object of type ProjectSetUp ####
rnaseqProject <- ProjectSetUp$new(
  
  date                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[1]],
  time                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[2]],
  projectName             = "RNASeq.RSEM",
  annotationRDS           = "C:/Users/sindiris/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS",
  pcRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/pc.other.HGNCTableFlat.rds",
  tfRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/TFs_no_epimachines.RDS",
  csRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/CellSurface.RDS",
  cgaRDS                  = "C:/Users/sindiris/R Scribble/Annotation RDS/cancerGermlineAntigens.rds",
  ewsr1Fli1RDS            = "C:/Users/sindiris/R Scribble/Annotation RDS/EWSR1_FL1_DownstreamTargets.RDS",
  pax3Foxo1RDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/PAX3_FOXO1_DownstreamTargets.RDS",
  
  BrainExpRDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Brain.v2.RDS",
  HeartExpRDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Heart.v2.RDS", 
  KidneyExpRDS            = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Kidney.v2.RDS", 
  LiverExpRDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Liver.v2.RDS", 
  LungExpRDS              = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Lung.v2.RDS", 
  
  outputPrefix            = "landscape",
  filterGenes             = TRUE,
  filterGeneMethod        = "bySum",
  factorName              = "DIAGNOSIS.Substatus.Tumor.Normal.Tissue",
  metadataFileRefCol      = "Sample.Biowulf.ID.GeneExp",
  metaDataFileName        = "MetadataMapper.v3.txt",
  outputdirRDSDir         = "GeneRDSOutput",
  outputdirTXTDir         = "GeneTXTOutput",
  gseaDir                 = "GSEA",
  plotsDir                = "Figures",
  plotsDataDir            = "FigureData",
  DiffGeneExpAnaDir       = "DiffExpResults",
  DiffGeneExpRDS          = "DiffGeneExpRDSOutput",
  ## Keep only PolyA
  #factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"),
  #                           "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "PolyA"),
  #                               "Tumors"=list("LIBRARY_TYPE"="Tumor", "LibraryPrep" = "PolyA"))
  ## Keep only Ribozero
  factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), "Normal.ribozero"=list("LibraryPrep" = "Ribozero"))
  ## Remove Celllines
  ##factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"))
)

## Add utility functions to the project ####
corUtilsFuncs <- CoreUtilities$new(  ProjectSetUpObject = rnaseqProject )

## Generate expression matrix ####
rm(mergeObjectsNoDup)
mergeObjectsNoDup <- corUtilsFuncs$getMergedMatrix(dir               = "TPM_Genes.v1",
                                                   fileFormat        = "txt",
                                                   colNameSelect     = "expected_count",
                                                   isRowNames        = TRUE,
                                                   rowNamesColInFile = 1,
                                                   fileSuffix        = ".genes.results",
                                                   primaryID         = "gene_id",
                                                   metadata          = rnaseqProject$metaDataDF,
                                                   metadataFileRefCol=rnaseqProject$metadataFileRefCol
                                                  )
saveRDS(mergeObjectsNoDup, "../RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.PolyA.RDS")

#mergeObjectsNoDup <- readRDS("../RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.excluding celllines.RDS")

## Evaluate presence of duplicate features (genes) and consolidate them ####
setDT(mergeObjectsNoDup, keep.rownames = TRUE)
mergeObjectsNoDup.pre <- mergeObjectsNoDup          %>% 
                         dplyr::rename(GeneID = rn) 
mergeObjectsNoDup.pre <- dplyr::left_join(rnaseqProject$annotationDF[,c("GeneID", "GeneName")], mergeObjectsNoDup.pre, by="GeneID") %>% 
                         data.table()
mergeObjectsConso     <- corUtilsFuncs$consolidateDF(mergeObjectsNoDup.pre[,-c("GeneID")], funcName = "max", featureName = "GeneName")
mergeObjectsConso     <- dplyr::full_join(mergeObjectsConso, rnaseqProject$annotationDF[,c("GeneID", "GeneName")], by="GeneName") %>%  
                         data.table()
mergeObjectsConso     <- subset(mergeObjectsConso,!duplicated(mergeObjectsConso$GeneName))
mergeObjectsConso     <- mergeObjectsConso[complete.cases(mergeObjectsConso), ]; dim(mergeObjectsConso)
mergeObjectsConso     <- mergeObjectsConso[,-c("GeneName")]         %>% 
                         data.frame()                               %>% 
                         tibble::column_to_rownames(var = "GeneID") %>% 
                         as.matrix() ; dim(mergeObjectsConso)
## matching above data frame with the annotationDF
rnaseqProject$annotationDF <- rnaseqProject$annotationDF %>% dplyr::filter(GeneID %in% rownames(mergeObjectsConso)); dim(rnaseqProject$annotationDF)

## Subset metaDataDF by the number of samples in the folder ####
colnamesDF           <- data.frame( "Sample.Biowulf.ID.GeneExp"= colnames(mergeObjectsConso))
corUtilsFuncs$subsetMetaData(colnamesDF=colnamesDF)

## Instantiate a new Object of type GeneExpNormalization ####
expressionObj        <- GeneExpNormalization$new(
  
  countObj          = as.matrix(mergeObjectsConso), 
  featureType       = "Gene", 
  packageRNAseq     = "edgeR", 
  annotationDF      = rnaseqProject$annotationDF, 
  design            = rnaseqProject$metaDataDF[,rnaseqProject$factorName], 
  #design           = newMetaDataDF[,rnaseqProject$factorName],
  proteinCodingOnly = FALSE,
  corUtilsFuncs     = corUtilsFuncs
)

## Get expression in desired units ####
### RawCounts
expressionTMM.Counts          = expressionObj$edgeRMethod("RawCounts")
### RPKM
expressionTMM.RPKM            = expressionObj$edgeRMethod("TMM-RPKM", logtransform = TRUE, zscore = FALSE)

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Biowulf.ID.GeneExp"=colnames(expressionTMM.RPKM)), 
                                                   rnaseqProject$validMetaDataDF[,c("Sample.Biowulf.ID.GeneExp", "Sample.ID.Alias")] )
AliasColnames                 <- c(as.character(AliasNames_df[c(1:7),1]), as.character(AliasNames_df[-c(1:7),2]))

## Perform Sanity Check for the above operations #####
stopifnot( length(colnames(expressionTMM.RPKM)) == length(AliasColnames) )
colnames(expressionTMM.RPKM)  <- AliasColnames

## Save expression (TMM-RPKM/whatwever asked for in the above step) to a file ####
write.table(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirTXTDir,"RPKM",
                                      paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.log2.polyA",rnaseqProject$date,".txt"),sep="/"),
                                      sep="\t", row.names = FALSE, quote = FALSE)
saveRDS(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirRDSDir,"RPKM",
                                      paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.log2.PolyA",rnaseqProject$date,".rds"),sep="/"))

### Performing ssGSEA output analysis. ( Plotting the scores across histology ) ##########

### Prepare input for ssGSEA broad gene pattern
expressionTMM.RPKM.GSEA.Input <- expressionTMM.RPKM[, -c(1:7)]; rownames(expressionTMM.RPKM.GSEA.Input) <- expressionTMM.RPKM[,6]
expressionTMM.RPKM.GSEA.print = corUtilsFuncs$createBroadGCTFile(expressionTMM.RPKM.GSEA.Input)

## Save input for ssGSEA 
write.table(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$gseaDir,
                                                 paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.zscore.",rnaseqProject$date,".txt"),sep="/"),
                                                 sep="\t", row.names = FALSE, quote = FALSE)
saveRDS(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$gseaDir,
                                             paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.zscore.",rnaseqProject$date,".rds"),sep="/"))

## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.zscore2019-01-31.PROJ.gct")

## Add custom expression like cytolytic scre and HLA gene expression to the ssGSEA Outpuut file.
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)
HLA_cytolyticScore      <- rbind(expressionTMM.RPKM.GSEA.Input[c("HLA-A", "HLA-B", "HLA-C"),], cytolyticScore)
ssGSEAScores.HLA.Cyto   <- rbind(ssGSEAScores,HLA_cytolyticScore)

## Plot the one variable plot
## Sanity Check: Checking metadata vs data ##
stopifnot( ncol(ssGSEAScores.HLA.Cyto) == length(as.character(rnaseqProject$metaDataDF$Sample.Biowulf.ID.GeneExp)) )

## Filter specified Diagnosis
factorsToExclude              = paste(c("NS", "YST", "Teratoma"), collapse = "|")
selected.metadata              <- rnaseqProject$metaDataDF  %>% filter_(  .dots = paste0("!grepl(", "'", factorsToExclude , "'" ,",", rnaseqProject$factorName, ")")) %>% 
                                                    dplyr::select_( .dots=c(rnaseqProject$metadataFileRefCol, rnaseqProject$factorName ) )
ssGSEAScores.HLA.Cyto.Selected <- ssGSEAScores.HLA.Cyto %>% dplyr::select_(.dots = selected.metadata[, rnaseqProject$metadataFileRefCol])
dim(ssGSEAScores.HLA.Cyto.Selected)

## sanity check Checking metadata vs data ##
stopifnot( ncol(ssGSEAScores.HLA.Cyto.Selected) == length(as.character(selected.metadata$Sample.Biowulf.ID.GeneExp)) )

## Preparing the expression matrix for string plot, by appending metadata
Scores <- cbind(t(ssGSEAScores.HLA.Cyto.Selected), selected.metadata[,rnaseqProject$factorName, drop=FALSE]) %>% 
          dplyr::rename_(.dots = setNames( list(rnaseqProject$factorName), list("Diagnosis") )) #%>%
#dplyr::mutate(Diagnosis = factor(Diagnosis, ordered = TRUE, levels = orderOfFactor))

### Setting up variables for  string plot
## Set the order of Diagnosis to appear
orderOfFactor    <- unique(Scores$Diagnosis)
## Set the order of signature to appear
orderOfSignature <- colnames(Scores)[-ncol(Scores)]
## Total list of signatures
colList          <- c(1:(ncol(Scores)-1))
## Generate custom colors
customColorDF    <- rnaseqProject$customColorsDF
## Plot the onevariable plot
plotLists        <- corUtilsFuncs$OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =TRUE, logit =FALSE, plotType = "density",
                                 yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF )
## Save the plots
EnrischmentScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName                <- paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir,"TMM-RPKM.ssGSEA.enrichmentScores.all.pc.log.zscore.pdf",sep="/")
## ggsave(SBName, marrangeGrob(EnrischmentScorePlots,ncol=2,nrow=1 ), width = 20, height = 10 )


## Plot to do percent samples enriched across cancer types

## Using  scores
dropSignatures    <- c("Macrophages_M0","Macrophages_M1", "Macrophages_M2","Dendritic_cells_activated")
factorsToExclude  <- paste(c("NS", "YST", "Teratoma"), collapse = "|")
## Read and parse the ssGSEA Output from Broad GenePattern
Scores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.gct")
## Standardizing the raw score to amplify the difference.
ScoresZscore      <- apply(Scores[1:24,],1, corUtilsFuncs$zscore_All)                  
ScoresZscore      %<>%   data.frame()                                                     %<>% 
                         tibble::rownames_to_column(var=rnaseqProject$metadataFileRefCol) %<>% 
                         dplyr::select(-one_of(dropSignatures ))

## Preparing the data and prepare for heatmap
ScoresForGather   <- tidyr::gather(ScoresZscore, key="GeneSet", value="Score", -!!rnaseqProject$metadataFileRefCol )
ScoresForGather   <-  dplyr::left_join(ScoresForGather, 
                                       rnaseqProject$metaDataDF[,c(rnaseqProject$metadataFileRefCol, rnaseqProject$factorName)], 
                                       by=rnaseqProject$metadataFileRefCol)                                                        %>% 
                      dplyr::filter_(  .dots = paste0("!grepl(", "'", factorsToExclude , "'" ,",", rnaseqProject$factorName, ")")) %>% 
                      dplyr::rename_(.dots = setNames(list(rnaseqProject$factorName),c("Diagnosis")) )

## Final check before plotting                        
dim(ScoresForGather);head(ScoresForGather)

ScoresForGatherPercent         <- ScoresForGather                                             %>% 
                                  dplyr::group_by(Diagnosis, GeneSet)                         %>% 
                                  dplyr::mutate(TotalCount = n(), Enriched = sum(Score > 0 )) %>% 
                                  dplyr::mutate(SamplePercent = (Enriched/TotalCount)*100 )   
ScoresForGatherUnique          <- ScoresForGatherPercent[,c(2,4,7)]                           %>% 
                                  ungroup()                                                   %>% 
                                  distinct()
ScoresForSpread                <- tidyr::spread( ScoresForGatherUnique, Diagnosis, SamplePercent ) %>% t() 
colnames(ScoresForSpread)      <- ScoresForSpread[1,]; 
ScoresForSpreadHeat            <- ScoresForSpread[-1,]
ScoresForSpreadHeat            <- t(apply(ScoresForSpreadHeat, 1, as.numeric))
colnames(ScoresForSpreadHeat)  <- colnames(ScoresForSpread)
#write.table(ScoresForSpread, "C:/Users/sindiris/R Scribble/RNASeq/PlotData/ScoresperDiagSigEnrichZscore.txt", sep="\t", quote = F, col.names = T, row.names = F)

## Open the PDF File
pdf( paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir,"PercentSamplesEnrichmentScoreHeatMap.Cibersort.pdf", sep="/"), height=10, width = 20)

### Using two different packages to plot
## Using SuperheatMap package
      # superheat(ScoresForSpreadHeat,
      #           bottom.label.text.angle=90,
      #           title.size = 6,
      #           heat.pal = c( "#4FFC07","#273746", "#F92908"),
      #           pretty.order.rows = T,
      #           pretty.order.cols = T
      # )

## Using fheatmap (presently using)
breaks <- seq(min(ScoresForSpreadHeat),max(ScoresForSpreadHeat), by=0.1)
#Yellow blue grey
#matrix_color_vector <- colorpanel(n=length(breaks)-1,low="#F4D03F",mid="#273746",high="#5DADE2")
# #green black red
# matrix_color_vector <- colorpanel(n=length(breaks)-1,low="#4FFC07",mid="#0B0B0B",high="#F92908")
#black red
matrix_color_vector <- colorpanel(n=length(breaks)-1,low="#4FFC07",mid="#273746",high="#F92908")
fheatmap(ScoresForSpreadHeat, display_tree_col = F,cluster_rows = T, mat_color = matrix_color_vector,
         row_fontsize = 5, col_fontsize = 5, cell_border = T, cell_border_col = "#A6ACAF",seed = 10,
         clustering_method = "complete", title = "Percent samples enriched across cancer types")
dev.off()

## Perform Differential gene expression analysis ####

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

tumorSubStatus.polyA    <- c("RMS.FP" , "RMS.FN", "EWS" ,"ASPS", "DSRCT", "HBL", "ML", "NB.MYCN.NA","NB.MYCN.A", "NB.Unknown", "OS", 
                          "SS", "Teratoma" ,"UDS" ,"YST")
tumorSubStatus.ribozero <-  c("WT" ,"CCSK")
Tumors                  <-  c("ASPS","DSRCT", "EWS" ,"HBL", "ML", "NB" ,"OS", "RMS", "SS", "Teratoma" ,"UDS" ,"YST","WT", "CCSK")


## Testing 
dgeObj  <- DifferentialGeneExp$new(
  countObj          = expressionObj$edgeRMethod("NormFactorDF")$counts,
  group1            = list(list("Normals"=NormalsNoGermLine,each=FALSE)),
  group2            = list(list("Tumor"=tumorSubStatus.polyA, each=TRUE)),
  packageRNAseq     = "edgeR",
  groupColumnName   = rnaseqProject$factorName,
  metadataDF        = rnaseqProject$metaDataDF,
  samplesColumnName = "Sample.Biowulf.ID.GeneExp",
  expressionUnit    = "TMM-RPKM",
  featureType       = "Gene",
  writeFiles        = TRUE,
  fileDirs          = rnaseqProject$fileDirs,
  subsetGenes       = TRUE,
  corUtilsFuncs     = corUtilsFuncs 
)

DiffExpObj <- dgeObj$performDiffGeneExp()

head(DiffExpObj[[1]] %>% dplyr::arrange(-logFC))

### Filtering ####

# Step 0  Define fucnctions ####

getCountObjTXT <- function(fileName, colNumb=1, rowNames=1){
  print(paste(fileName))
  featureCountTxt <- read.csv(fileName, sep="\t", row.names = "GeneName", header = 1);
  return(featureCountTxt[,colNumb, drop=FALSE])
}

mergeDiffTestResults <- function(x, type="", saveDirPath="", extension="", colInterest=1, rowNamesCol =1,
                                 fileSuffix=".txt"){
  
  print(paste(x))
  file_Dir_Gene = x
  fileName <- basename(file_Dir_Gene) 
  dir.create(file.path(paste(saveDirPath,fileName,sep="/")))
  GeneFiles             <- list.files(file_Dir_Gene); GeneFiles <- GeneFiles[grep(fileSuffix, GeneFiles)]
  GeneFilesList         <- paste(file_Dir_Gene, "/", GeneFiles,sep="") ; length(GeneFilesList)
  
  countObj          <- do.call(cbind,lapply(GeneFilesList, getCountObjTXT, colNumb=colInterest, rowNames=rowNamesCol))
  
  write.table(countObj, paste(saveDirPath, paste(fileName, "/", type,"_MergedDiffExpResult.txt",sep=""), sep= "/"), sep="\t",row.names = TRUE, quote = FALSE)
}

# Step 1  Set the filters and annotation ####

PValue = 0.001; FDR = 0.05

#selectedGeneList <- "CancerGermlineAntigen"
#group2FPKM = 1; Zscored.logFC = 0.25 ; Zscore.group2 = 0; 
# group2FPKM = 0 ; group1FPKM = 1;  PValue = 0.01 ; logFC =1 ; FDR = 0.05

selectedGeneList <- "CellSurface"
#group2FPKM = 40; Zscored.logFC = 1 ; Zscore.group2 = 1
#group2FPKM = 40 ; group1FPKM = 1;  PValue = 0.001 ; logFC =2 ; FDR = 0.05
group2FPKM = 2 ; group1FPKM = 2;  PValue = 0.001 ; logFoldDiff =3 ; FDR_value = 0.001 ; vitalFPKM = 1

#selectedGeneList <- "TranscriptionFactor"
#group2FPKM = 1; Zscored.logFC = 1.25 ; Zscore.group2 = 0
#group2FPKM = 1 ; group1FPKM = 5;  PValue = 0.01 ; logFC =1 ; FDR = 0.05

# Step 2  Perform MErging ####

MergedDiffExpResultDir <- paste0("C:/Users/sindiris/R Scribble//RNASeq.RSEM/MergedDiffExpResults/",selectedGeneList)
dir.create(MergedDiffExpResultDir)
#ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })), c("Normals_WT", "Normals_CCSK") )
ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })))
groups <- list.dirs(paste("C:/Users/sindiris/R Scribble//RNASeq.RSEM//DiffExpResults/", sep=""))[-1]; groups[1]
output <- sapply(groups, mergeDiffTestResults, type="Gene", colInterest=c(7,9,10,11,12, 15:28), rowNamesCol = 2,
                 fileSuffix=paste0(selectedGeneList,".txt"),saveDirPath=MergedDiffExpResultDir)

# Step 3  Core Function and save files ####
allTumorStats <- do.call(cbind, lapply(ConditionGroup, function(x){
  tumorData <- read.csv( paste(MergedDiffExpResultDir,"/",x,"/Gene_MergedDiffExpResult.txt",sep=""), sep="\t", header = T, stringsAsFactors = FALSE ) 
 
  ## Actual filtering
  groupsCompare <- unlist(strsplit(x, "_"))
  print(groupsCompare)
  filterDFByColNames <- c("logFC",	groupsCompare[2], groupsCompare[1])
  newColNames <- paste0("Zscored.",c("logFC",	groupsCompare[2], groupsCompare[1]))

  ##Zscoreing matrix
  tumorDataPvalue        <- tumorData %>% tibble::rownames_to_column(var="GeneName") ; print(dim(tumorData))
  tumorDataPvalue_Zscore <- apply(tumorDataPvalue[,c("logFC",	groupsCompare[2], groupsCompare[1])],2,corUtilsFuncs$zscore_All)
  colnames(tumorDataPvalue_Zscore) <- newColNames;
  tumorDataPvalue_Zscore <- cbind(tumorDataPvalue[,c("GeneName"),drop = FALSE], tumorDataPvalue_Zscore ) %>% data.frame()

  # print(colnames(tumorDataPvalue_Zscore))    ; print(head(tumorDataPvalue_Zscore))

  tumorAllData <- left_join(tumorDataPvalue_Zscore, tumorDataPvalue, by="GeneName") ; dim(tumorAllData)
  # tumorAllData.filt <- tumorAllData %>% dplyr::filter_(.dots=paste0(groupsCompare[2]," >= ", group2FPKM ,
  #                               " & ","Zscored.logFC >= ", Zscored.logFC,
  #                               " & ", paste0("Zscored.",groupsCompare[2]), " >= ", Zscore.group2)) %>%
  #                 dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
  
  write.table(tumorAllData, paste(MergedDiffExpResultDir,"/",x,"/",x,".rankFile.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  tumorAllData.filt <- tumorAllData %>% dplyr::filter_(.dots=paste0(       groupsCompare[2], " >= ", group2FPKM ,
                                                                    " & ", groupsCompare[1], " <= ", group1FPKM ,
                                                                    " &   logFC >", logFoldDiff,
                                                                    " &   FDR   <", FDR_value ,
                                                                    " &  Brain.MeanExp  < ", vitalFPKM ,
                                                                    " &  Heart.MeanExp  < ", vitalFPKM ,
                                                                    " &  Kidney.MeanExp < ", vitalFPKM ,
                                                                    " &  Liver.MeanExp  < ", vitalFPKM  ,
                                                                    " &  Lung.MeanExp   < ", vitalFPKM  )) %>%
                                                                    dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
  
  write.table(tumorAllData.filt, paste(MergedDiffExpResultDir,"/",x,"/",x,".filtered.rankFile.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  # ## select genes
  selectGenes <- tumorAllData.filt %>%  dplyr::select(GeneName)

  tumorDataPvalue["status"] <- 0
  statusDF <- tumorDataPvalue %>% mutate(status=ifelse(GeneName %in% selectGenes$GeneName, 1, 0)) %>% dplyr::select(GeneName, status) %>%
                            rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
  
  return(statusDF)
}))
tumorStatusDF <- allTumorStats[, !duplicated(colnames(allTumorStats))]  %>% mutate(RowSum= rowSums(.[-1]))

# Step 4. Perform Clustering using "cluster_data" method from "fheatmap" ####
tumorStatusDF.HM <- tumorStatusDF %>% tibble::column_to_rownames(var="GeneName") %>% dplyr::select(-one_of("RowSum"))
tumorStatusDF.HM.memo <- corUtilsFuncs$memoSort(M=tumorStatusDF.HM)
tumorStatusDF.HM.memo$RowSum <- apply(tumorStatusDF.HM.memo, 1, function(x) sum(x!=0))
tumorStatusDF.HM.memo %<>% tibble::rownames_to_column(var="GeneName")
write.table(tumorStatusDF.HM.memo, paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.DExp.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
# roworder <- unlist( cluster_data(tumorStatusDF.HM, distance = "euclidean", method = "ward.D")["order"])
# colorder <- unlist( cluster_data(t(tumorStatusDF.HM), distance = "euclidean", method = "ward.D")["order"])

# Step 5. Organise the genes and tumors as per the above order . And save the file. ####


# Step 6. Select rows for heatmap ####
allTumorStatsFinal <- read.table( paste(MergedDiffExpResultDir, "/", selectedGeneList,".Summarised.DExp", ".txt" ,sep=""),sep="\t", header = TRUE)
CTA.Filt <- allTumorStatsFinal %>% filter(RowSum>=2) %>% 
  dplyr::arrange(RowSum)
dim(CTA.Filt)
# %>% filter(GeneName %in% c("CD99", "FGFR4", "ALK"))

# Step 8. Plot the heatmap ####
CTA.Filt %<>% dplyr::select(-one_of("RowSum"))
CTA.Filt %<>%  column_to_rownames(var="GeneName") 
colnames(CTA.Filt) <- gsub("NormalsStatus", "", colnames(CTA.Filt))

#pdf( paste("./Plots/",date, ".Differentially Expressed CGAs.v17.pdf", sep=""), height = 19, width = 15)
# superheat(t(CTA.Filt), pretty.order.cols =T,
#           #title = "Differentially Expressed CGAs",
#           #linkage.method = "ward.D2",
#           legend=FALSE,
#           grid.hline = FALSE,
#           grid.vline = FALSE,
#           X.text.size = 20,
#           # grid.hline.size = 0.01,
#           # grid.vline.size = 0.01,
#           # heat.col.scheme = "grey",
#           heat.lim = c(0, 1),
#           #heat.pal = c("#004080",  "#88cc00"),
#           heat.pal = c("#e0e0d1", "#004080"),
#           bottom.label.text.angle=90,
#           title.size = 6)
# dev.off()
pdf( paste("C:/Users/sindiris/R Scribble//RNASeq.RSEM/Figures/", 
  "Differentially Expressed ",  selectedGeneList, ".pdf", sep=""), height = 10, width = 25)
pheatmap(t(CTA.Filt), color =c("#e0e0d1", "#004080"), 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = NA, 
         treeheight_row = 0,
         fontsize = 12,
         legend = FALSE )
dev.off()

###### Step 2 ####
# allTumorStats <- do.call(cbind, lapply(ConditionGroup[1:2], function(x){
#   tumorData <- read.table( paste(MergedDiffExpResultDir,"/",x,"/Gene_MergedDiffExpResult.txt",sep=""), sep="\t",
#                            row.names = 1, header = T, stringsAsFactors = FALSE )
#   groupsCompare <- unlist(strsplit(x, "_"))
#   print(groupsCompare)
#   selectGenes <- tumorData %>% filter_(.dots=paste0(groupsCompare[2]," >= ", group2FPKM ," & ",groupsCompare[1]," <= ", group1FPKM ," & ",
#                                                     "PValue <= ", PValue ," & ", "logFC  >= ", logFC," & ", "FDR  <= ", FDR)) %>%
#     dplyr::arrange_(.dots = c("logFC") ) %>%
#     dplyr::select(GeneName.x)
#   tumorData["status"] <- 0
#   statusDF <- tumorData %>% mutate(status=ifelse(GeneName.x %in% selectGenes$GeneName.x, 1, 0)) %>% dplyr::select(GeneName.x, status) %>%
#     rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
# 
#   print("status DF")
#   print(head(statusDF))
#   
#   tumorStatusDF <- statusDF[, !duplicated(colnames(statusDF))]  %>% mutate(RowSum= rowSums(.[-1]))
#   write.table(tumorStatusDF, paste(MergedDiffExpResultDir,"/",x,"/","SummarisedDExpDF.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
#   returnTumorDF <- tumorStatusDF %>% dplyr::select(GeneName.x, RowSum) %>% rename(c('RowSum'=paste(groupsCompare[2],"StatusSum", sep="")))
#   
#   print("returnTumorDF DF")
#   print(head(returnTumorDF))
# }))













