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
  #                          "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "PolyA"),
  #                              "Tumors"=list("LIBRARY_TYPE"="Tumor", "LibraryPrep" = "PolyA"))
  ## Keep only Ribozero
  # factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), "Normal.ribozero"=list("LibraryPrep" = "Ribozero"))
  ## Remove Celllines
  factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"))
)

## Add utility functions to the project ####
corUtilsFuncs <- CoreUtilities$new(  ProjectSetUpObject = rnaseqProject )

# Make a Tree Map ####
StatsFinal <-  rnaseqProject$metaDataDF %>% group_by_(.dots= c("DIAGNOSIS.Alias","DIAGNOSIS.Alias.TreeMap",rnaseqProject$factorName, "Color", "LIBRARY_TYPE.TreeMap") ) %>% 
  count_(var=as.name("Sample.ID")) %>% dplyr::summarise(Count=n()) %>% 
  dplyr::group_by_(.dots= c(rnaseqProject$factorName)) %>%  
  dplyr::mutate( SampleSum := sum(Count)) %>% 
  spread_("LIBRARY_TYPE.TreeMap", "Count") %>% 
  mutate_( .dots = setNames( list( interp(~paste(rnaseqProject$factorName ,"(", Sum , ")"), 
                                          factorName=as.name(rnaseqProject$factorName), Sum=as.name("SampleSum") ) ), "LegendSampleSum") ) %>% 
  data.frame() %>% distinct(DIAGNOSIS.Alias.TreeMap,SampleSum, .keep_all = TRUE)

StatsFinal <- StatsFinal %>% dplyr::filter(DIAGNOSIS.Alias != "NS")
StatsFinal[,"LegendSampleSum"] <- paste(StatsFinal[,"DIAGNOSIS.Alias.TreeMap"],"( ",StatsFinal[,"SampleSum"], " )",sep="")
#pdf(file=paste(Plots, date, "Diagnosis Tree Map All",date,"pdf",sep="."), height=8, width= 10)
treemap(dtf=data.frame(StatsFinal), index=c("DIAGNOSIS.Alias", "DIAGNOSIS.Alias.TreeMap"),
        vSize="SampleSum",
        type="categorical",
        vColor="LegendSampleSum",
        palette = as.character(StatsFinal$Color),
        fontcolor.labels=c("black"),
        bg.labels=c("#CCCCCCDC"),
        algorithm = "squarified",
        inflate.labels=F,
        fontsize.labels = 10,
        fontsize.legend = 10,
        border.lwds=0.9,
        title = "Samples Map",
        title.legend = "Histology"
)
dev.off()

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
#saveRDS(mergeObjectsNoDup, "../RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.RiboZeros.RDS")

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
# write.table(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirTXTDir,"RPKM",
#                                       paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.log2.RiboZero",rnaseqProject$date,".txt"),sep="/"),
#             sep="\t", row.names = FALSE, quote = FALSE)
# saveRDS(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirRDSDir,"RPKM",
#                                   paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.log2.RiboZero",rnaseqProject$date,".rds"),sep="/"))
# 

### Performing ssGSEA output analysis. ( Plotting the scores across histology ) ##########

### Prepare input for ssGSEA broad gene pattern
expressionTMM.RPKM.GSEA.Input <- expressionTMM.RPKM[, -c(1:7)]; rownames(expressionTMM.RPKM.GSEA.Input) <- expressionTMM.RPKM[,6]
expressionTMM.RPKM.GSEA.print = corUtilsFuncs$createBroadGCTFile(expressionTMM.RPKM.GSEA.Input)

# ## Save input for ssGSEA 
# write.table(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$gseaDir,
#                                       paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.zscore.",rnaseqProject$date,".txt"),sep="/"),
#                                       sep="\t", row.names = FALSE, quote = FALSE)
# saveRDS(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$gseaDir,
#                                       paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.zscore.",rnaseqProject$date,".rds"),sep="/"))

## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")

## Add custom expression like cytolytic scre and HLA gene expression to the ssGSEA Outpuut file.
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)
HLA_cytolyticScore      <- rbind(expressionTMM.RPKM.GSEA.Input[c("HLA-A", "HLA-B", "HLA-C"),], cytolyticScore)
data.frame(colnames(HLA_cytolyticScore), colnames(ssGSEAScores))
ssGSEAScores.HLA.Cyto   <- rbind(ssGSEAScores,HLA_cytolyticScore)

## Plot the one variable plot
## Sanity Check: Checking metadata vs data ##
stopifnot( ncol(ssGSEAScores.HLA.Cyto) == length(as.character(rnaseqProject$metaDataDF$Sample.Biowulf.ID.GeneExp)) )

## Filter specified Diagnosis
factorsToExclude              = paste(c("NS", "YST", "Teratoma"), collapse = "|")
selected.metadata              <- rnaseqProject$metaDataDF  %>% 
                                  filter_(  .dots = paste0("!grepl(", "'", factorsToExclude , "'" ,",", rnaseqProject$factorName, ")")) %>% 
                                  dplyr::select_( .dots=c(rnaseqProject$metadataFileRefCol, rnaseqProject$factorName ) )
ssGSEAScores.HLA.Cyto.Selected <- ssGSEAScores.HLA.Cyto %>% dplyr::select_(.dots = as.character(selected.metadata[, rnaseqProject$metadataFileRefCol]))
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
ggsave(SBName, marrangeGrob(EnrischmentScorePlots,ncol=2,nrow=1 ), width = 20, height = 10 )


## Plot to do percent samples enriched across cancer types

## Using  scores
dropSignatures    <- c("Macrophages_M0","Macrophages_M1", "Macrophages_M2","Dendritic_cells_activated")
factorsToExclude  <- paste(c("NS", "YST", "Teratoma"), collapse = "|")
## Read and parse the ssGSEA Output from Broad GenePattern
Scores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
## Standardizing the raw score to amplify the difference.
ScoresZscore      <- apply(Scores[c(1:24,43),],1, corUtilsFuncs$zscore_All)                  
ScoresZscore      %<>% data.frame() %<>% tibble::rownames_to_column(var=rnaseqProject$metadataFileRefCol) %<>% dplyr::select(-one_of(dropSignatures ))

## Preparing the data and prepare for heatmap
rnaseqProject$metaDataDF <- as.data.frame( apply(rnaseqProject$metaDataDF, 2, as.character), stringsAsFactors = FALSE )
ScoresForGather   <- tidyr::gather(ScoresZscore, key="GeneSet", value="Score", -!!rnaseqProject$metadataFileRefCol )
ScoresForGather   <- dplyr::left_join(ScoresForGather, 
                                      rnaseqProject$metaDataDF[,c(rnaseqProject$metadataFileRefCol, rnaseqProject$factorName)], 
                                      by=rnaseqProject$metadataFileRefCol) %>% 
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
fheatmap(t(ScoresForSpreadHeat), display_tree_col = F,cluster_rows = T, mat_color = matrix_color_vector,
         row_fontsize = 5, col_fontsize = 5, cell_border = F, cell_border_col = "#A6ACAF",seed = 10,
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
  countObj_print    <- countObj %>% tibble::rownames_to_column(var="GeneName")
  
  write.table(countObj_print, paste(saveDirPath, paste(fileName, "/", type, ".DiffExp.txt",sep=""), sep= "/"), sep="\t",
              row.names = FALSE, quote = FALSE)
}

# Step 1.  Set the filters and annotation ####

## Javed's Filter for all three categories
#group2FPKM.T = 2 ; group1FPKM.T = 2;  PValue.T = 0.001 ; logFoldDiff.T =3 ; FDR_value.T = 0.05 ; vitalFPKM.T = 1

selectedGeneList <- "CancerGermlineAntigen"
# #group2FPKM = 1;
# Zscored.logFC = 0.25 ; Zscore.group2 = 0; group2FPKM = 0 ; group1FPKM = 1;  PValue = 0.01 ; logFC =1 ; FDR = 0.05

# selectedGeneList <- "CellSurface"
# Zscored.logFC = 1 ; Zscore.group2 = 1; group2FPKM = 40 ; group1FPKM = 1;  PValue = 0.001 ; logFC =2 ; FDR = 0.05

## Zscore Filtering
#selectedGeneList <- "TranscriptionFactor"
#Zscored.logFC = 1 ; Zscore.group2 = 0.5; group2FPKM = 2 ; #group1FPKM = 5;  PValue = 0.01 ; logFC =1 ; FDR = 0.05

MergedDiffExpResultDir <- paste0("C:/Users/sindiris/R Scribble//RNASeq.RSEM/MergedDiffExpResults/",selectedGeneList)

# Step 2.  Perform Merging of differential expression file across groups ####

dir.create(MergedDiffExpResultDir)
ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })), c("Normals_WT", "Normals_CCSK") )
#ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })))
groups <- list.dirs(paste("C:/Users/sindiris/R Scribble//RNASeq.RSEM//DiffExpResults/", sep=""))[-1]; groups[1]
output <- sapply(groups, mergeDiffTestResults, type="Gene", colInterest=c(7,9,10,11,12, 15:28), rowNamesCol = 2,
                 fileSuffix=paste0(selectedGeneList,".txt"),saveDirPath=MergedDiffExpResultDir)

# Step 3.  Core Function and save files ####
allTumorStats <- do.call(cbind, lapply(ConditionGroup, function(x){
  tumorData <- read.csv( paste(MergedDiffExpResultDir,"/",x,"/Gene.DiffExp.txt",sep=""), sep="\t", header = T, stringsAsFactors = FALSE ) 
  
  ## Actual filtering
  groupsCompare <- unlist(strsplit(x, "_"))
  print(groupsCompare)
  filterDFByColNames <- c("logFC",	groupsCompare[2], groupsCompare[1])
  newColNames <- paste0("Zscored.",c("logFC",	groupsCompare[2], groupsCompare[1]))
  
  ##Zscoreing matrix
  tumorDataPvalue        <- tumorData ; print(dim(tumorData))
  tumorDataPvalue_Zscore <- apply(tumorDataPvalue[,c("logFC",	groupsCompare[2], groupsCompare[1])],2,corUtilsFuncs$zscore_All)
  colnames(tumorDataPvalue_Zscore) <- newColNames;
  tumorDataPvalue_Zscore <- cbind(tumorDataPvalue[,c("GeneName"),drop = FALSE], tumorDataPvalue_Zscore ) %>% data.frame()
  tumorAllData <- left_join(tumorDataPvalue_Zscore, tumorDataPvalue, by="GeneName") ; dim(tumorAllData)
  
  write.table(tumorAllData, paste(MergedDiffExpResultDir,"/",x,"/",x,".allgenes.DiffExp.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  ## Zscore Ranking filter
  tumorAllData.zscore <- tumorAllData %>% 
    dplyr::filter_(.dots=paste0( 
      groupsCompare[2]," >= ", group2FPKM ,
      " &  Zscored.logFC   >= ", Zscored.logFC,
      " & ", paste0("Zscored.",groupsCompare[2]), " >= ", Zscore.group2)) %>%
    dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
  
  ## Complete filtered gene List with zscoring filter
  write.table(tumorAllData.zscore, paste(MergedDiffExpResultDir,"/",x,"/",x,".filteredgenes.ZscoringRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  ## Complete filtered gene List (Binary) with zscoring filter
  selectGenes <- tumorAllData.zscore %>%  dplyr::select(GeneName)
  tumorData.zscoreRanking <- tumorDataPvalue
  tumorData.zscoreRanking["status"] <- 0
  statusDF.zscoreRanking <- tumorData.zscoreRanking %>% mutate(status=ifelse(GeneName %in% selectGenes$GeneName, 1, 0)) %>% dplyr::select(GeneName, status) %>%
    rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
  statusDF.zscoreRanking <- statusDF.zscoreRanking[, !duplicated(colnames(statusDF.zscoreRanking))] 
  
  # write.table(statusDF.zscoreRanking, paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.zscoreRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  
  ## Javed's Filtering
  tumorAllData.filt <- tumorAllData %>% dplyr::filter_(.dots=paste0(       groupsCompare[2], " >= ", group2FPKM.T ,
                                                                           " & ", groupsCompare[1], " <= ", group1FPKM.T ,
                                                                           " &   logFC >", logFoldDiff.T,
                                                                           " &   FDR   <", FDR_value.T ,
                                                                           " &  Brain.MeanExp  < ", vitalFPKM.T ,
                                                                           " &  Heart.MeanExp  < ", vitalFPKM.T ,
                                                                           " &  Kidney.MeanExp < ", vitalFPKM.T ,
                                                                           " &  Liver.MeanExp  < ", vitalFPKM.T  ,
                                                                           " &  Lung.MeanExp   < ", vitalFPKM.T  )) %>%
    dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
  
  ## Complete filtered gene List with traditional filter
  write.table(tumorAllData.filt, paste(MergedDiffExpResultDir,"/",x,"/",x,".filteredgene.traditionalRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  ## Complete filtered gene List (Binary) with traditional filter 
  selectGenes <- tumorAllData.filt %>%  dplyr::select(GeneName)
  tumorData.traditionalRanking <- tumorDataPvalue
  tumorData.traditionalRanking["status"] <- 0
  statusDF.traditionalRanking <- tumorData.traditionalRanking %>% mutate(status=ifelse(GeneName %in% selectGenes$GeneName, 1, 0)) %>% dplyr::select(GeneName, status) %>%
    rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
  statusDF.traditionalRanking <- statusDF.traditionalRanking[, !duplicated(colnames(statusDF.traditionalRanking))] 
  
  # write.table(statusDF.traditionalRanking, paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  
  return(list(statusDF.zscoreRanking,statusDF.traditionalRanking))
}))

allInOneFile <- do.call(rbind, lapply(ConditionGroup, function(x){
  tumorData <- read.csv( paste(MergedDiffExpResultDir,"/",x,"/Gene.DiffExp.txt",sep=""), sep="\t", header = T, stringsAsFactors = FALSE ) 
  colnames(tumorData)[6] <- "Tumor"
  tumorData$Group <- x
  return(tumorData)
})); dim(allInOneFile)
write.table(allInOneFile, paste(MergedDiffExpResultDir,"/",selectedGeneList,".allSamples.txt",  sep=""),
            sep="\t", row.names = FALSE, quote = FALSE)

#tumorStatusDF <- allTumorStats[, !duplicated(colnames(allTumorStats))]  %>% mutate(RowSum= rowSums(.[-1]))

# Step 4.  Merge multiple DFs, memo-sort each DF and plot ####
allTumorMergedStats <- lapply(1:nrow(allTumorStats), function(x){
  mergedDF <- do.call(cbind, allTumorStats[x,])
  mergedDF <- mergedDF[, !duplicated(colnames(mergedDF))]        %>% 
    tibble::column_to_rownames(var="GeneName")
  mergedDF.Memo        <- corUtilsFuncs$memoSort(M=mergedDF)
  mergedDF.Memo$RowSum <- apply(mergedDF.Memo, 1, function(x) sum(x!=0))
  mergedDF.Memo %<>% tibble::rownames_to_column(var="GeneName")
  return(mergedDF.Memo)
})

# Step 5.  Save the files ####
write.table(allTumorMergedStats[[1]], paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.ZscoreRank.Dexp.txt",  sep=""),
            sep="\t", row.names = FALSE, quote = FALSE)
write.table(allTumorMergedStats[[2]], paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.Dexp.txt",sep=""),
            sep="\t", row.names = FALSE, quote = FALSE)

# Step 6.  Select rows for heatmap ####
allTumorStatsFinal <- read.table(paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.ZscoreRank.Dexp.txt",sep=""),sep="\t", header = TRUE)
CTA.Filt <- allTumorStatsFinal %>% filter(RowSum>=3) %>% 
  dplyr::arrange(RowSum)
dim(CTA.Filt)
# %>% filter(GeneName %in% c("CD99", "FGFR4", "ALK"))

# Step 8.  Plot the heatmap ####
CTA.Filt %<>% dplyr::select(-one_of("RowSum"))
CTA.Filt %<>%  column_to_rownames(var="GeneName") 
colnames(CTA.Filt) <- gsub("NormalsStatus", "", colnames(CTA.Filt))

pdf( paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir, "zscore.Differentially Expressed TF.v18.pdf", sep="/"), height = 10, width = 25)
superheat(t(CTA.Filt), pretty.order.cols =T,
          #title = "Differentially Expressed CGAs",
          #linkage.method = "ward.D2",
          legend=FALSE,
          grid.hline = FALSE,
          grid.vline = FALSE,
          X.text.size = 20,
          # grid.hline.size = 0.01,
          # grid.vline.size = 0.01,
          # heat.col.scheme = "grey",
          heat.lim = c(0, 1),
          #heat.pal = c("#004080",  "#88cc00"),
          heat.pal = c("#e0e0d1", "#004080"),
          bottom.label.text.angle=90,
          title.size = 6)
dev.off()

pdf( paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir, "zscore.Differentially Expressed cs.v19.pdf", sep="/"), height = 15, width = 25)
pheatmap(t(CTA.Filt), color =c("#e0e0d1", "#004080"), 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = NA, 
         treeheight_row = 0,
         fontsize = 12,
         legend = FALSE )
dev.off()

### Performing TCR analysis

## Perfroming TCR analysis ####

### Placeholder DF ####
emptyDF <- data.frame(count=c(0), freq=c(0), cdr3nt=c("None"),cdr3aa=c("NF"),v=c("NF"),d=c("NF"),j=c("NF"),VEnd=c(0),DStart=c(0),
                      DEnd=c(0),JStart=c(0),SampleName=c(0))
emptyDFEntropy <- data.frame(VJcombo=c(), Counts =c(), Vcassette=c(), Jcassette=c(), aaCDR3_filtered = c(), ntCDR3= c())

emptyDFEntropyResults <- data.frame(FileName=c(), Hcdr3 =c(), Htot=c(), CLcdr3=c(), CLHvj= c(), CLtot= c(),
                                    Hcdr3_max=c(), Hvj_max =c(), Htot_max=c(), CLcdr3_max=c(), Num_CDR3= c(), Num_VJ= c(),
                                    Num_totCDR3 =c())

readCloneFiles <- function(x, cloneType=NA){
  print(x)
  exomeData <- read.csv( paste(TCRDir, x, sep=""), sep="\t", header = TRUE )
  if(nrow(exomeData)>0){
    exomeData$SampleName <- x
    if( !cloneType %in% unique(substr(exomeData$v,1,3)) ) {
      emptyDF$SampleName <- c(x)
      exomeData <- rbind(exomeData, emptyDF)
    }
    
  } else {
    emptyDF$SampleName <- c(x)
    exomeData <- emptyDF
  }
  
  return(exomeData)
}

## Start analysis for clone type: Choose clone type ####
## cloneType = "IGH";
cloneType = "TRB";
correlationPlots <- function(varName="", constName="", df=NA, customColorDF=NA, xlab="Log Total Clones"){
  
  print(paste(varName))
  customColorsVector <- setNames( as.character(customColorDF$Color), as.character(customColorDF$Diagnosis))
  corrTest <- cor.test(df[,constName], df[,varName], method = "spearman")
  if  ( corrTest$p.value < 2.2e-16 ) { corrTest$p.value = 2.2e-16 }
  plot <- ggplot(df, aes_string(x=constName, y=varName)) + 
    geom_smooth(method=lm,  fill="grey") +
    geom_point(aes(colour = factor(Diagnosis)), show.legend = T, size=3, shape=16) + 
    scale_colour_manual(values=customColorsVector) +
    theme_bw() +
    theme(axis.text=element_text(size=13)
          ,axis.title=element_text(size=13,face="bold")) +
    xlab(xlab) +
    ylab(paste("Standardised Enrichment Score", sep=" "))+
    ggtitle(paste("Corr.Coeff = ", signif(corrTest$estimate[[1]],5), "\np-value = ", signif(corrTest$p.value,5), "", varName,sep=""))
  
  return(list(plot))
}   
customColorDF    <- rnaseqProject$customColorsDF
TCRResultsDir <- paste0(rnaseqProject$workDir,rnaseqProject$projectName,"/TCR.Results/")

### List files and read data into a single data matrix ####
TCRDir <- paste0(rnaseqProject$workDir,rnaseqProject$projectName,"/TCR.clones.files/")
fileList <- list.files(TCRDir) ; length(fileList)
AllClonesData             <- rbindlist( lapply(fileList, readCloneFiles, cloneType=cloneType) ) ; dim(AllClonesData)

### Filter Clones by clone types ####
cloneObj               <- corUtilsFuncs$filterSpecificCloneTypes(cloneData = AllClonesData, cloneType = cloneType) %>% 
                          dplyr::rename(Sample.ID=SampleName) %>% 
                          dplyr::mutate(Sample.ID = gsub("Sample_|convert.|.clones.txt","", Sample.ID)) %>% 
                          data.frame() ; dim(cloneObj); head(cloneObj)
cloneObj.Expansion.GE3 <- cloneObj %>%  dplyr::filter(grepl(cloneType,v) & count >= 3); dim(cloneObj); head(cloneObj)

### Attach metadata and generate countObj ####
#countObj <- countObj %>% dplyr::rename(Sample.ID=SampleName); 
#countObj$Sample.ID <- gsub("Sample_|convert.|.clones.txt","", countObj$Sample.ID)
#countObj$SAMPLE_ID <- gsub("-","_", countObj$SAMPLE_ID)
countObj.Annot        <- dplyr::full_join(cloneObj, rnaseqProject$metaDataDF, by="Sample.ID") %>% 
                         dplyr::select_(.dots=c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", 
                                                "DStart", "DEnd", "JStart", "Sample.ID", "Sample.ID.Alias",
                                                "LIBRARY_TYPE",rnaseqProject$factorName)) %>% 
                          dplyr::filter(complete.cases(.))
dim(countObj.Annot);head(countObj.Annot)

### Section 2 Aggregating and summarizing ####

### Summarize By Samples 
countObj.Annot.gb.Samples  <- countObj.Annot  %>% dplyr::group_by(Sample.ID) %>% 
                                                  dplyr::summarise(
                                                  TotalCloneSum=sum(count),
                                                  TotalClones=ifelse(n()==1 & TotalCloneSum==0, 0, n()) )  %>% 
                                                  dplyr::rename_(.dots=setNames(list("TotalClones"),c(cloneType)))
dim(countObj.Annot.gb.Samples); head(countObj.Annot.gb.Samples); tbl_df(countObj.Annot.gb.Samples)

countObj.Annot.gb.Samples.Annotate  <- left_join(countObj.Annot.gb.Samples, 
                                                 rnaseqProject$metaDataDF[,c("Sample.ID","Sample.ID.Alias","LIBRARY_TYPE", rnaseqProject$factorName)], 
                                                 by="Sample.ID") %>% dplyr::filter(complete.cases(.))
###Replace NA by 0
countObj.Annot.gb.Samples.Annotate[which(is.na(countObj.Annot.gb.Samples.Annotate$TotalCloneSum)), c(cloneType,"TotalCloneSum")] <- 0
countObj.Annot.gb.Samples.Annotate[which(countObj.Annot.gb.Samples.Annotate$TotalCloneSum == 0), cloneType] <- 0

### Saving files
# saveRDS(countObj.Annot.gb.Samples, paste0(TCRResultsDir,"countObj.Annot.gb.Samples.Annotate",".", cloneType,".RDS" , sep="") )
# write.table(countObj.Annot.gb.Samples.Annotate,paste0(TCRResultsDir,"countObj.Annot.gb.Samples",".", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)

### Summarize By Diagnosis
countObj.Annot.gb.Diagnosis  <- countObj.Annot  %>% dplyr::group_by_(.dots= list("cdr3aa", "v", "d", "j", rnaseqProject$factorName)) %>% 
                                                    summarise(
                                                      TotalSamples=n(), 
                                                      Samples = paste(Sample.ID, collapse = ','), 
                                                      CloneCount = paste(count, collapse = ','),
                                                      MedianCloneCount = median(count)
                                                    ) # %>% 
                                                      # dplyr::filter(complete.cases(.))
# saveRDS(countObj.Annot.gb.Diagnosis, paste0(TCRResultsDir,"countObj.Annot.gb.Diagnosis.Annotate",".", cloneType,".RDS" , sep=""))
# write.table(countObj.Annot.gb.Diagnosis, paste0(TCRResultsDir,"countObj.Annot.gb.Diagnosis",".", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)

### Section 3 Filter normal TCR and Samples ####
countObj.Annot.CL.Normal <- countObj.Annot %>% dplyr::filter(LIBRARY_TYPE %in% c("Normal", "CellLine")); dim(countObj.Annot.CL.Normal)
cdr3_Normal_CellLine     <- countObj.Annot.CL.Normal %>% filter(cdr3aa != c("NF")) %>% distinct(cdr3aa)
countObj.Annot.Tumor     <- countObj.Annot %>% 
                                    dplyr::filter( !LIBRARY_TYPE %in% c("Normal", "CellLine") ) %>% 
                                    dplyr::filter( !cdr3aa %in% c(cdr3_Normal_CellLine$cdr3aa)); dim(countObj.Annot.Tumor) 
## Summarise
countObj.gb.Samples.Annotate.NoNS  <- countObj.Annot.Tumor  %>% dplyr::group_by(Sample.ID) %>% 
                                                  dplyr::summarise(
                                                  TotalCloneSum=sum(count),
                                                  TotalClones=ifelse(n()==1 & TotalCloneSum==0, 0, n()) )  %>% 
                                                  dplyr::rename_(.dots=setNames(list("TotalClones"),c(cloneType)))
dim(countObj.gb.Samples.Annotate.NoNS); head(countObj.gb.Samples.Annotate.NoNS); tbl_df(countObj.gb.Samples.Annotate.NoNS)

countObj.gb.Samples.Annotate.NoNS  <- left_join(countObj.gb.Samples.Annotate.NoNS, 
                                                 rnaseqProject$metaDataDF[,c("Sample.ID","Sample.ID.Alias","LIBRARY_TYPE", rnaseqProject$factorName)], 
                                                 by="Sample.ID") %>% dplyr::filter(complete.cases(.))
## Change to factors
countObj.gb.Samples.Annotate.NoNS     <- countObj.gb.Samples.Annotate.NoNS %>% 
                                          dplyr::mutate(LIBRARY_TYPE=factor(countObj.gb.Samples.Annotate.NoNS$LIBRARY_TYPE, 
                                                               levels = unique(countObj.gb.Samples.Annotate.NoNS$LIBRARY_TYPE))) %>%
                                          dplyr::mutate(DIAGNOSIS.Substatus.Tumor.Normal.Tissue=factor(countObj.gb.Samples.Annotate.NoNS$DIAGNOSIS.Substatus.Tumor.Normal.Tissue, 
                                                               levels = unique(countObj.gb.Samples.Annotate.NoNS$DIAGNOSIS.Substatus.Tumor.Normal.Tissue)))
# saveRDS(countObj.gb.Samples.Annotate.NoNS, paste0(TCRResultsDir,"countObj.gb.Samples.Annotate.NoNS",".", cloneType,".RDS" , sep=""))
# write.table(countObj.gb.Samples.Annotate.NoNS, paste0(TCRResultsDir,"countObj.gb.Samples.Annotate.NoNS",".", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)

### plot for TCR COunt Bean plot

### Prepare data for one variable plot ####

selectCol="TotalCloneSum" ; StatsFinalCol=rnaseqProject$factorName ; SampleNames <- "Sample.ID.Alias"
tcrcloneCountPre          <- countObj.gb.Samples.Annotate.NoNS %>% 
                                dplyr::select_(.dots=c(paste0("selectCol"), paste0("StatsFinalCol"), paste0("SampleNames")))
  
tcrcloneCountPre.Diag   <- tcrcloneCountPre %>%  dplyr::rename_(.dots = setNames(list(SampleNames,StatsFinalCol),c("Samples","Diagnosis"))) 
ScoresPre               <- tcrcloneCountPre.Diag[,!(colnames(tcrcloneCountPre.Diag) %in% c("Samples")), drop=FALSE]
orderOfFactor           <- as.character( unique(ScoresPre$Diagnosis) )
orderOfSignature        <- colnames(ScoresPre)[-ncol(ScoresPre)]
colList                 <- c(1:(ncol(ScoresPre)-1)) ; Scores <- ScoresPre
## Generate custom colors
customColorDF    <- rnaseqProject$customColorsDF

## Filter for diagnosis
Scores <- Scores %>% filter(!Diagnosis %in% c("Teratoma", "YST"))

### Plot and Save ###
plotLists <- corUtilsFuncs$OneVariablePlotSort( colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =FALSE, logit =TRUE, logBase=10,
                                                yLab = "log(IGH counts)", legendDisplay = FALSE, customColorDF = customColorDF, 
                                                plotType = "StringBean", sizeOfDots = 0.6  )
plotLists
tcrcloneCountPlots <- lapply(plotLists, function(l) l[[1]])
tcrcloneCountData  <- lapply(plotLists, function(l) l[[2]]) %>% data.frame(check.names = FALSE) %>% bind_rows()
  
SBName = paste0(TCRResultsDir,"/",cloneType, ".BeanPlot.v3.pdf")
ggsave(SBName, marrangeGrob(tcrcloneCountPlots, ncol=1, nrow=1), width = 10, height = 5)
dev.off()


### Coorelation with Immune Signature

### Prepare data for correlation between immunescore and clone count Read the enrichment score data ####
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.gct")
ssGSEA.zscore <- apply(ssGSEAScores, 1, corUtilsFuncs$zscore_All) ; 
ssGSEA.t <- ssGSEA.zscore %>% data.frame() %>% tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp"); 
dim(ssGSEA.t)
ssGSEA.t <- left_join(ssGSEA.t, 
                      rnaseqProject$metaDataDF[,c("Sample.ID","Sample.Biowulf.ID.GeneExp","Sample.ID.Alias","LIBRARY_TYPE", rnaseqProject$factorName)], 
                      by="Sample.Biowulf.ID.GeneExp") %>% dplyr::filter(complete.cases(.)); dim(ssGSEA.t)

immuneScore.Clones <- left_join(countObj.gb.Samples.Annotate.NoNS[,c("TotalCloneSum",cloneType, "Sample.ID")], ssGSEA.t, by="Sample.ID") %>% 
                      dplyr::select(-one_of(c("Sample.Biowulf.ID.GeneExp","Sample.ID"))) %>% 
                      data.frame()
#immuneScore.Clones %<>% tibble::column_to_rownames("Sample.ID")
immuneScore.Clones <- immuneScore.Clones %>% dplyr::rename_(.dots = setNames(list(rnaseqProject$factorName),c("Diagnosis"))) %>% dplyr::mutate(TotalCloneSum = log10(TotalCloneSum+1))

varNames <- colnames(immuneScore.Clones[,3:44])  
plotLists <- lapply(varNames, correlationPlots, constName="TotalCloneSum",  df= data.frame(immuneScore.Clones), customColorDF=customColorDF)
ImmuneScorePlots <- lapply(plotLists, function(l) l[[1]] )

SBName =paste0(TCRResultsDir,"/ImmuneScore.vs.TotalCloneSum.",cloneType,".pdf")
ggsave(SBName, marrangeGrob(ImmuneScorePlots, ncol=1, nrow=1), width = 15, height = 10)
dev.off()

### Prepare data to compare TCR with public databases

### For Venn plots Public dataBases ####
## Group by TCR
countObj.Annot.gb <- countObj.Annot  %>% dplyr::group_by(cdr3aa, v, d, j) %>% 
  dplyr::summarise(
    TotalSamples=n(), 
    Samples = paste(Sample.ID, collapse = ','), 
    CloneCount = paste(count, collapse = ','),
    MedianCloneCount = median(count),
    Diagnosis= paste(unique(DIAGNOSIS.Substatus.Tumor.Normal.Tissue), collapse = ',' )
  )

## Using VDJdb
vdjdb <- read.csv("../RNASeq.RSEM/TCR.Public.DB/vdjdb-2018-01-17/vdjdb_all_CDR3aa.txt", sep="\t")
vdjdbHealthy <- vdjdb %>% filter(grepl("healthy",meta.subject.cohort))
vdjdb.CDR3Beta <- unique(as.character(vdjdbHealthy$cdr3.beta))
length(vdjdb.CDR3Beta)
vdjdb.CDR3Beta.Complete <- unlist(stringr::str_extract_all(vdjdb.CDR3Beta, "^C.*F")); 
length(vdjdb.CDR3Beta.Complete)

## Using Waren
warenetal <- read.csv("../RNASeq.RSEM/TCR.Public.DB/waren_et_al/waren_et_al_all_CDR3aa.txt", sep="\t", header = F)
warenetal.CDR3Beta <- unique(as.character(warenetal$V1))
length(warenetal.CDR3Beta)
warenetal.CDR3Beta.Complete <- unlist(stringr::str_extract_all(warenetal.CDR3Beta, "^C.*F")); 
length(warenetal.CDR3Beta.Complete)


## All TCGA ## Using Bo et al/TCGA 
TCGA.Tumor.Normal <- read.csv("../RNASeq.RSEM/TCR.Public.DB/TCGA/Bo.et.al.CDR3.fa", sep="\t", header = F); dim(TCGA.Tumor.Normal)
TCGATumorNormalTab <- data.frame(Names=gsub(">","",TCGA.Tumor.Normal[seq(1, 1366836,2),]), Sequence=TCGA.Tumor.Normal[seq(2, 1366836,2),])
head(TCGATumorNormalTab); dim(TCGATumorNormalTab)
#Normal
TCGANormal <- TCGATumorNormalTab %>% dplyr::filter(grepl('^N', Names)) ; dim(TCGANormal)
#Tumor
TCGATumor  <- TCGATumorNormalTab %>% dplyr::filter(grepl('^T', Names)) ; dim(TCGATumor)
## Getting Only Shared TCGA Tumor CDR3aa
TCGATumor.Shared <- TCGATumor %>% dplyr::filter( grepl(".*C[A|S]+.*F", Sequence)) %>%
                                  dplyr::group_by(Sequence) %>% 
                                  dplyr::mutate(Samples = paste(unique(Names), collapse = ','), Count = n(),
                                                      DistinctSamples = length(unique(Names))) %>% dplyr::distinct()
TCGATumor.Shared.Multi <- TCGATumor.Shared %>% dplyr::filter(DistinctSamples > 1)
TCGATumor.Shared.Multi$Length <- sapply(as.character(TCGATumor.Shared.Multi$Sequence), nchar)
TCGATumor.Shared <- TCGATumor %>% dplyr::filter( grepl(".*C[A|S]+.*F", Sequence))
TCGATumorDistinct <- unique(TCGATumor.Shared$Sequence)

## Getting Only Shared TCGA Normal CDR3aa
TCGANormal.Shared <- TCGANormal %>% dplyr::filter( grepl(".*C[A|S]+.*F", Sequence)) %>%
  dplyr::group_by(Sequence) %>% 
  dplyr::mutate(Samples = paste(unique(Names), collapse = ','), Count = n(),
                DistinctSamples = length(unlist(strsplit(Samples, ",")))) %>% dplyr::distinct()
TCGANormal.Shared.Multi <- TCGANormal.Shared %>% dplyr::filter(DistinctSamples > 1)
TCGANormal.Shared.Multi$Length <- sapply(as.character(TCGANormal.Shared.Multi$Sequence), nchar)
TCGANormal.Shared <- TCGANormal %>% dplyr::filter( grepl(".*C[A|S]+.*F", Sequence))
TCGANormalDistinct <- unique(TCGANormal.Shared$Sequence)

## countObj.Annot.gb
countObj.Annot.gb.Cancer <- countObj.Annot.gb %>% filter(!grepl("NS",Diagnosis))
countObj.Annot.gb.Cancer$length <- sapply(as.character(countObj.Annot.gb.Cancer$cdr3aa), nchar)
cancerShared             <- countObj.Annot.gb.Cancer %>% filter(TotalSamples > 1); length(cancerShared$cdr3aa)
cancerShared.Complete    <- unlist(stringr::str_extract_all(cancerShared$cdr3aa, "^C.*F")); length(cancerShared.Complete)
cancerPrivate            <- countObj.Annot.gb.Cancer %>% filter(TotalSamples == 1); length(cancerPrivate$cdr3aa)
cancerPrivate.Complete    <- unlist(stringr::str_extract_all(cancerPrivate$cdr3aa, "^C.*F")); length(cancerPrivate.Complete)

countObj.Annot.gb.NS <- countObj.Annot.gb %>% filter(grepl("NS",Diagnosis)); length(countObj.Annot.gb.NS$cdr3aa)
inHouseNormal.Complete    <- unlist(stringr::str_extract_all(countObj.Annot.gb.NS$cdr3aa, "^C.*F")); length(inHouseNormal.Complete)

vennCDR3aaList <- list("Tumor, Private" = unique(cancerPrivate$cdr3aa),
                       "Tumor, Shared" = unique(cancerShared$cdr3aa),
                       "Tumor, Public" = TCGATumorDistinct,
                       "TCGA,  Normal"    =   TCGANormalDistinct,
                       "Warren et al (healthy)" = warenetal.CDR3Beta, 
                       "Chudakov et al (healthy)"=vdjdb.CDR3Beta,
                       "In House Normal" = unique(countObj.Annot.gb.NS$cdr3aa) )

# vennCDR3aaList <- list("Tumor, Private" = unique(cancerPrivate.Complete),
#                        "Tumor, Shared" = unique(cancerShared.Complete),
#                        "Warreb et al (healthy)" = unique(warenetal.CDR3Beta.Complete), 
#                        "Chudakov et al (healthy)"=unique(vdjdb.CDR3Beta.Complete),
#                        "In House Normal" = unique(inHouseNormal.Complete) 
# )

library(venn)
v.table <- venn::venn(vennCDR3aaList, ilab=TRUE, zcolor = "style", size = 15, cexil = 1, cexsn = 1)

### Plot Inforgraphics ####
TumorPrivate <- attr(v.table,"intersections")[["Tumor, Private"]]; length(TumorPrivate)
TumorShared  <- attr(v.table,"intersections")[["Tumor, Shared"]]; length(TumorShared)
InHouseNormals  <- attr(v.table,"intersections")[["In House Normal"]]; length(InHouseNormals)
warenetalNormals  <- attr(v.table,"intersections")[["Warren et al (healthy)"]]; length(warenetalNormals)


TumorPrivateFinal <- data.frame("AA"=TumorPrivate, "length"=sapply(TumorPrivate,nchar))
colnames(TumorPrivateFinal) <- c("AA","Len")
TumorPrivateFinal_15 <- TumorPrivateFinal %>% filter(Len == 15)
TumorPrivatePlot <- ggseqlogo(as.character(TumorPrivateFinal_15$AA), seq_type='aa',  method = "probability")

warenetalNormalsFinal <- data.frame("AA"=warenetalNormals, "length"=sapply(warenetalNormals,nchar))
colnames(warenetalNormalsFinal) <- c("AA","Len")
warenetalNormalsFinal <- warenetalNormalsFinal %>% filter(Len == 15)
warenetalNormalsPlot <- ggseqlogo(as.character(warenetalNormalsFinal$AA), seq_type='aa',  method = "probability")


InHouseNormalsFinal <- data.frame("AA"=InHouseNormals, "length"=sapply(InHouseNormals,nchar))
colnames(InHouseNormalsFinal) <- c("AA","Len")
InHouseNormalsFinal <- InHouseNormalsFinal %>% filter(Len == 15)
InHouseNormalsPlot <- ggseqlogo(as.character(InHouseNormalsFinal$AA), seq_type='aa',  method = "probability")

pdf("ven.probability.pdf", height = 10, width = 10)
ggarrange(plotlist = list(TumorPrivatePlot, InHouseNormalsPlot, warenetalNormalsPlot),common.legend=TRUE, nrow = 3)
dev.off()

### Clonality #####

ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
ssGSEA.zscore <- apply(ssGSEAScores, 1, corUtilsFuncs$zscore_All) ; 
ssGSEA.t <- ssGSEA.zscore %>% data.frame() %>% tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp"); 
dim(ssGSEA.t)
ssGSEA.t <- left_join(ssGSEA.t, 
                      rnaseqProject$metaDataDF[,c("Sample.Data.ID", "Sample.Biowulf.ID.GeneExp", rnaseqProject$factorName)], 
                      by="Sample.Biowulf.ID.GeneExp") %>% dplyr::filter(complete.cases(.)); dim(ssGSEA.t)

entropy <- read.table(paste0("../RNASeq.RSEM/TCR.Clones.Entropy/AllEntropyData_H_CL_JS.landscape.",cloneType,".v3.txt"), sep="\t", header = T)
entropyMeta <- dplyr::left_join(entropy, rnaseqProject$metaDataDF[, c("Sample.Data.ID", "LIBRARY_TYPE")], by="Sample.Data.ID")
entropyMeta.Filt <- entropyMeta %>% dplyr::filter(!LIBRARY_TYPE %in% c("Normal", "CellLine") )

entropyMetassGSEA <- dplyr::left_join(entropyMeta.Filt, ssGSEA.t, by="Sample.Data.ID") 
entropyMetassGSEA <- entropyMetassGSEA[complete.cases(entropyMetassGSEA), ] %>% dplyr::rename_(.dots = setNames(list(rnaseqProject$factorName),c("Diagnosis")))

## Correlation Plot
varNames <- colnames(entropyMetassGSEA[,16:58]) 
plotLists <- lapply(varNames, correlationPlots, constName="Htot..Entropy.", xlab="Entropy", df= data.frame(entropyMetassGSEA), customColorDF=customColorDF)
ImmuneScorePlots <- lapply(plotLists, function(l) l[[1]] )
SBName =paste0(TCRResultsDir,"/ImmuneScore.vs.Htot..Entropy",cloneType,".pdf")
ggsave(SBName, marrangeGrob(ImmuneScorePlots, ncol=1, nrow=1), width = 15, height = 10)
dev.off()

## Bean plot
entropyScores <- entropyMetassGSEA[,c("Htot..Entropy.", "Diagnosis")] 
orderOfFactor           <- as.character( unique(entropyScores$Diagnosis) )
orderOfSignature        <- colnames(entropyScores)[-ncol(entropyScores)]
colList                 <- c(1:(ncol(entropyScores)-1)) ; Scores <- entropyScores
## Generate custom colors
customColorDF    <- rnaseqProject$customColorsDF

## Filter for diagnosis
Scores <- entropyScores %>% filter(!Diagnosis %in% c("Teratoma", "YST")) %>% dplyr::filter(complete.cases(.))

### Plot and Save ###
plotLists <- corUtilsFuncs$OneVariablePlotSort( colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =FALSE, logit =TRUE, logBase=10,
                                                yLab = "log( Entropy )", legendDisplay = FALSE, customColorDF = customColorDF, 
                                                plotType = "StringBean", sizeOfDots = 0.6  )

