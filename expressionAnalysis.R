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
  annotationRDS           = "C:/Users/sindiris/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS",
  pcRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/pc.other.HGNCTableFlat.rds",
  tfRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/TFs_no_epimachines.RDS",
  csRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/CellSurface.RDS",
  cgaRDS                  = "C:/Users/sindiris/R Scribble/Annotation RDS/cancerGermlineAntigens.rds",
  ewsr1Fli1RDS            = "C:/Users/sindiris/R Scribble/Annotation RDS/EWSR1_FL1_DownstreamTargets.RDS",
  pax3Foxo1RDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/PAX3_FOXO1_DownstreamTargets.RDS",
  
  BrainExpRDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Brain.RDS",
  HeartExpRDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Heart.RDS", 
  KidneyExpRDS            = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Kidney.RDS", 
  LiverExpRDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Liver.RDS", 
  LungExpRDS              = "C:/Users/sindiris/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Lung.RDS", 
  
  outputPrefix            = "landscape",
  filterGenes             = TRUE,
  filterGeneMethod        = "bySum",
  factorName              = "DIAGNOSIS.Substatus.Tumor.Normal.Tissue",
  metaDataFileName        = "MetadataMapper.v3.txt",
  outputdirRDSDir         = "GeneRDSOutput",
  outputdirTXTDir         = "GeneTXTOutput",
  gseaDir                 = "GSEA",
  plotsDir                = "Figures",
  plotsDataDir            = "FigureData",
  DiffGeneExpAnaDir       = "DiffExpResults",
  DiffGeneExpRDS          = "DiffGeneExpRDSOutput",
  #factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), 
  #                               "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "PolyA"),
  #                               "Tumors"=list("LIBRARY_TYPE"="Tumor", "LibraryPrep" = "PolyA"))
  factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "Ribozero"))                                     
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
                                                   metadataFileRefCol= "Sample.Biowulf.ID")

## Evaluate presence of duplicate features and consolidate them
setDT(mergeObjectsNoDup, keep.rownames = TRUE)
mergeObjectsNoDup.pre <- mergeObjectsNoDup %>% dplyr::rename(GeneID = rn)
mergeObjectsNoDup.pre <- dplyr::left_join(rnaseqProject$annotationDF[,c("GeneID", "GeneName")], mergeObjectsNoDup.pre, by="GeneID") %>% data.table()
mergeObjectsConso <- corUtilsFuncs$consolidateDF(mergeObjectsNoDup.pre[,-c("GeneID")], funcName = "max", 
                                                 featureName = "GeneName")
mergeObjectsConso <- dplyr::full_join(mergeObjectsConso, rnaseqProject$annotationDF[,c("GeneID", "GeneName")], by="GeneName") %>% 
                          data.table()
mergeObjectsConso <-   subset(mergeObjectsConso,!duplicated(mergeObjectsConso$GeneName))
mergeObjectsConso <-   mergeObjectsConso[complete.cases(mergeObjectsConso), ]; dim(mergeObjectsConso)
mergeObjectsConso <- mergeObjectsConso[,-c("GeneName")] %>% data.frame() %>% tibble::column_to_rownames(var = "GeneID") %>% as.matrix()
rnaseqProject$annotationDF <- rnaseqProject$annotationDF %>% dplyr::filter(GeneID %in% rownames(mergeObjectsConso))
  
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
  proteinCodingOnly = FALSE,
  corUtilsFuncs     = corUtilsFuncs
)

## Get expression in desired units
expressionTMM.RPKM = expressionObj$edgeRMethod("TMM-RPKM")
expressionTMM.Counts = expressionObj$edgeRMethod("RawCounts")

## Start here ##
AliasNames_df  <- dplyr::left_join( data.frame("SAMPLE_ID"=colnames(expressionTMM.RPKM)), rnaseqProject$validMetaDataDF[,c("SAMPLE_ID", "SAMPLE_ID.Alias")] )
AliasColnames  <- c(as.character(AliasNames_df[c(1:7),1]), as.character(AliasNames_df[-c(1:7),2]))
## Start here ##
stopifnot( length(colnames(expressionTMM.RPKM)) == length(AliasColnames) )
colnames(expressionTMM.RPKM) <- AliasColnames
write.table(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirTXTDir,"RPKM",
                                      paste0("RPKM_Data_Filt_Consolidated.GeneNames.",rnaseqProject$date,".txt"),sep="/"),
                                      sep="\t", row.names = FALSE, quote = FALSE)
saveRDS(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirRDSDir,"RPKM",
                                      paste0("RPKM_Data_Filt_Consolidated.GeneNames.",rnaseqProject$date,".rds"),sep="/"))

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

tumorSubStatus.polyA <- c("RMS.FP" , "RMS.FN", "EWS" ,"ASPS", "DSRCT", "HBL", "ML", "NB.MYCN.NA","NB.MYCN.A", "NB.Unknown", "OS", 
                          "SS", "Teratoma" ,"UDS" ,"YST")
tumorSubStatus.ribozero <-  c("WT" ,"CCSK")
Tumors         <-  c("ASPS","DSRCT", "EWS" ,"HBL", "ML", "NB" ,"OS", "RMS", "SS", "Teratoma" ,"UDS" ,"YST","WT", "CCSK")


## Testing 
dgeObj  <- DifferentialGeneExp$new(
  countObj          = expressionObj$edgeRMethod("NormFactorDF")$counts,
  group1            = list(list("Normals"=NormalsNoGermLine,each=FALSE)),
  group2            = list(list("Tumor"=tumorSubStatus.polyA, each=TRUE)),
  packageRNAseq     = "edgeR",
  groupColumnName   = rnaseqProject$factorName,
  metadataDF        = rnaseqProject$metaDataDF,
  samplesColumnName = "SAMPLE_ID",
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

getCountObjTXT <- function(fileName, colNumb=1, rowNames=1){
  print(paste(fileName))
  featureCountTxt <- read.csv(fileName, sep="\t", row.names = rowNames, header = 1);
  return(featureCountTxt[,colNumb, drop=FALSE])
}

PValue = 0.001; FDR = 0.05

# Step 1  Set the filters and annotation ####

#selectedGeneList <- "cancergermlineantigen"
#group2FPKM = 1; Zscored.logFC = 0.25 ; Zscore.group2 = 0; 
# group2FPKM = 0 ; group1FPKM = 1;  PValue = 0.01 ; logFC =1 ; FDR = 0.05

#selectedGeneList <- "cellsurface"
#group2FPKM = 40; Zscored.logFC = 1 ; Zscore.group2 = 1
#group2FPKM = 40 ; group1FPKM = 1;  PValue = 0.001 ; logFC =2 ; FDR = 0.05

selectedGeneList <- "transcriptionFactor"
group2FPKM = 1; Zscored.logFC = 1.25 ; Zscore.group2 = 0
#group2FPKM = 1 ; group1FPKM = 5;  PValue = 0.01 ; logFC =1 ; FDR = 0.05

# Step 2  Perform MErging ####

MergedDiffExpResultDir <- paste0("C:/Users/sindiris/R Scribble//RNASeq.RSEM/MergedDiffExp/",selectedGeneList)
dir.create(MergedDiffExpResultDir)
ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })), c("Normals_WT", "Normals_CCSK") )
#ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })))
groups <- list.dirs(paste("C:/Users/sindiris/R Scribble//RNASeq.RSEM//DiffExpResults/", sep=""))[-1]; groups[1]
output <- sapply(groups, mergeDiffTestResults, type="Gene", colInterest=c(5,7,9,10,11,12), rowNamesCol = 5,
                 fileSuffix=paste0(selectedGeneList,".txt"),saveDirPath=MergedDiffExpResultDir)

# Step 3  Core Function and save files ####
allTumorStats <- do.call(cbind, lapply(ConditionGroup, function(x){
  tumorData <- read.table( paste(MergedDiffExpResultDir,"/",x,"/Gene_MergedDiffExpResult.txt",sep=""), sep="\t", 
                           row.names = 1, header = T, stringsAsFactors = FALSE )
 
  ## Actual filtering
  groupsCompare <- unlist(strsplit(x, "_"))
  print(groupsCompare)
  filterDFByColNames <- c("logFC",	groupsCompare[2], groupsCompare[1])
  newColNames <- paste0("Zscored.",c("logFC",	groupsCompare[2], groupsCompare[1]))
  tumorData <- read.table( paste(MergedDiffExpResultDir,"/",x,"/Gene_MergedDiffExpResult.txt",sep=""), 
                           sep="\t", row.names = 1, header = T, stringsAsFactors = FALSE )
  tumorDataPvalue <- tumorData %>% dplyr::filter(PValue <= 0.001); dim(tumorDataPvalue)
  tumorDataPvalue_Zscore <- apply(tumorDataPvalue[,c("logFC",	groupsCompare[2], groupsCompare[1])],2,corUtilsFuncs$zscore_All)
  colnames(tumorDataPvalue_Zscore) <- newColNames
  tumorDataPvalue_Zscore <- cbind(tumorDataPvalue[,c("GeneName.x"),drop = FALSE], tumorDataPvalue_Zscore)
  head(tumorDataPvalue_Zscore)
  tumorAllData <- left_join(tumorDataPvalue_Zscore, tumorData, by="GeneName.x") ; dim(tumorAllData)
  tumorAllData.filt <- tumorAllData %>% dplyr::filter_(.dots=paste0(groupsCompare[2]," >= ", group2FPKM ,
                                " & ","Zscored.logFC >= ", Zscored.logFC,
                                " & ", paste0("Zscored.",groupsCompare[2]), " >= ", Zscore.group2)) %>% 
                  dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
  
  write.table(tumorAllData, paste(MergedDiffExpResultDir,"/",x,"/","rankFile.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  ## select genes
    selectGenes <- tumorAllData.filt %>%  dplyr::select(GeneName.x)
  
  tumorData["status"] <- 0 
  statusDF <- tumorData %>% mutate(status=ifelse(GeneName.x %in% selectGenes$GeneName.x, 1, 0)) %>% dplyr::select(GeneName.x, status) %>% 
                            rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
  
  tumorStatusDF <- statusDF[, !duplicated(colnames(statusDF))]  %>% mutate(RowSum= rowSums(.[-1])) 
  write.table(tumorStatusDF, paste(MergedDiffExpResultDir,"/",x,"/","SummarisedDExpDF.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  returnTumorDF <- tumorStatusDF %>% dplyr::select(GeneName.x, RowSum) %>% rename(c('RowSum'=paste(groupsCompare[2],"StatusSum", sep="")))
  
  return(returnTumorDF)
}))

# Step 4. Perform Clustering using "cluster_data" method from "fheatmap" ####
allTumorStats <- allTumorStats[, !duplicated(colnames(allTumorStats))] ; rownames(allTumorStats) <- allTumorStats[,1]; allTumorStats <- allTumorStats[,-c(1)]
roworder <- unlist( cluster_data(allTumorStats, distance = "euclidean", method = "ward.D")["order"])
colorder <- unlist( cluster_data(t(allTumorStats), distance = "euclidean", method = "ward.D")["order"])

# Step 5. Organise the genes and tumors as per the above order . And save the file. ####
allTumorStatsFinal <- allTumorStats[roworder, colorder]
allTumorStatsFinal$OnesSum <- apply(allTumorStatsFinal, 1, function(x) sum(x!=0))
allTumorStatsFinal[c("CTAG1B", "PRAME"),] ; length(which(allTumorStatsFinal[,"OnesSum"] >=1 ))

# Step 6. Write the Final file ####
write.table(allTumorStatsFinal, paste(MergedDiffExpResultDir, "/", selectedGeneList, ".v9.","PValue-",PValue,".","FDR-",FDR,".",
                                      "NormalRPKM.Less.Than.1", ".txt" ,sep=""),
            sep="\t", row.names = TRUE, quote = FALSE)

# Step 7. Select rows for heatmap ####
allTumorStatsFinal <- read.table( paste(MergedDiffExpResultDir, "/", selectedGeneList, ".v9.","PValue-",PValue,".","FDR-",FDR,".",
                                        "NormalRPKM.Less.Than.1", ".txt" ,sep=""),sep="\t", header = TRUE)
CTA.Filt <- allTumorStatsFinal %>% rownames_to_column(var="CGA") %>% filter(OnesSum>=2) %>% #dplyr::select(-matches("ML|YST|UDS")) %>% 
  dplyr::arrange(OnesSum)
dim(CTA.Filt)
CTA.Filt %>% filter(CGA %in% c("CD99", "FGFR4", "ALK"))

# Step 8. Plot the heatmap ####
CTA.Filt %<>% dplyr::select(-one_of("OnesSum"))
CTA.Filt %<>%  column_to_rownames(var="CGA") 
colnames(CTA.Filt) <- gsub("StatusSum", "", colnames(CTA.Filt))

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
         clustering_method = "ward.D",
         cluster_rows = FALSE, 
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













