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
  annotationRDS           = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS",
  pcRDS                   = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/pc.other.HGNCTableFlat.rds",
  emRDS                   = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/EMGenes.RDS",
  tfRDS                   = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/TFs_no_epimachines.RDS",
  csRDS                   = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/CellSurface.RDS",
  cgaRDS                  = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/cancerGermlineAntigens.rds",
  ewsr1Fli1RDS            = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/EWSR1_FL1_DownstreamTargets.RDS",
  pax3Foxo1RDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/PAX3_FOXO1_DownstreamTargets.RDS",
  
  BrainExpRDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Brain.v2.RDS",
  HeartExpRDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Heart.v2.RDS", 
  KidneyExpRDS            = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Kidney.v2.RDS", 
  LiverExpRDS             = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Liver.v2.RDS", 
  LungExpRDS              = "T:/Sivasish_Sindiri/R Scribble/Annotation RDS/VitalExpression/expressionTMM.RPKM.Lung.v2.RDS", 
  
  outputPrefix            = "landscape",
  filterGenes             = TRUE,
  filterGeneMethod        = "bySum",
  factorName              = "DIAGNOSIS.Substatus.Tumor.Normal.Tissue",
  #factorName              = "DIAGNOSIS.Substatus.Tumor.Tissue",
  metadataFileRefCol      = "Sample.Biowulf.ID.GeneExp",
  metaDataFileName        = "MetadataMapper.v3.txt",
  outputdirRDSDir         = "GeneRDSOutput",
  outputdirTXTDir         = "GeneTXTOutput",
  gseaDir                 = "GSEA",
  plotsDir                = "Figures",
  plotsDataDir            = "FigureData",
  DiffGeneExpAnaDir       = "DiffExpResults",
  DiffGeneExpRDS          = "DiffGeneExpRDSOutput",
  ## Keep only Ribozero
  # factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"),
  #                          "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "PolyA"),
  #                              "Tumors"=list("LIBRARY_TYPE"="Tumor", "LibraryPrep" = "PolyA"))
  #Keep only PolyA
  # factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), "Normal.ribozero"=list("LibraryPrep" = "Ribozero"))
  ## Remove Celllines
  factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"))
  # factorsToExclude          = list('None'=list("LIBRARY_TYPE"=""))
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
#dev.off()

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
                                                   metadataFileRefCol=rnaseqProject$metadataFileRefCol )

# saveRDS(mergeObjectsNoDup, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.RDS")
# saveRDS(mergeObjectsNoDup, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.Cellline.RDS")

## Tumor Normal and no cell line
mergeObjectsNoDup <- readRDS("../RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.RDS")
## Tumor Normal and Cellline
## mergeObjectsNoDup <- readRDS("../RNASeq.RSEM/GeneRDSOutput/RawCount/All.samples.Tumor.Normal.Cellline.RDS")

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
  proteinCodingOnly = TRUE,
  corUtilsFuncs     = corUtilsFuncs
)


## Get expression in desired units ####
### RawCounts
#expressionTMM.Counts          = expressionObj$edgeRMethod("RawCounts")
## Normalised counts
#expressionTMM.NormDF         = expressionObj$edgeRMethod("NormFactorDF")
### RPKM
expressionTMM.RPKM            = expressionObj$edgeRMethod("TMM-RPKM", logtransform = FALSE, zscore = FALSE)

## Arrange data by histology and Library type
arrange_metadata <- rnaseqProject$validMetaDataDF %>% arrange(DIAGNOSIS.Substatus.Tumor.Normal.Tissue, desc(LIBRARY_TYPE))
expressionTMM.RPKM.arr <- expressionTMM.RPKM %>% dplyr::select(one_of("Chr","Start","End","Strand","GeneID","GeneName","Length",
                                                                      as.character(factor(arrange_metadata$Sample.Biowulf.ID.GeneExp, 
                                                                      ordered = TRUE, 
                                                                      levels = arrange_metadata$Sample.Biowulf.ID.GeneExp))))

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Biowulf.ID.GeneExp"=colnames(expressionTMM.RPKM)), 
                                                   rnaseqProject$validMetaDataDF[,c("Sample.Biowulf.ID.GeneExp", "Sample.ID.Alias", "Sample.Data.ID", 
                                                                                    "DIAGNOSIS.Alias",
                                                                                    rnaseqProject$factorName)] )
AliasColnames                 <- c(as.character(AliasNames_df[c(1:7),1]), as.character(AliasNames_df[-c(1:7),2]))


## Perform Sanity Check for the above operations #####
stopifnot( length(colnames(expressionTMM.RPKM)) == length(AliasColnames) )
colnames(expressionTMM.RPKM)  <- AliasColnames

### Save expression (TMM-RPKM/whatwever asked for in the above step) to a file ####
write.table(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirTXTDir,"RPKM",
                                      paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.log2.",rnaseqProject$date,".txt"),sep="/"),
            sep="\t", row.names = FALSE, quote = FALSE)
saveRDS(expressionTMM.RPKM, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$outputdirRDSDir,"RPKM",
                                  paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.log2.",rnaseqProject$date,".rds"),sep="/"))



### Performing ssGSEA output analysis. ( Plotting the scores across histology ) ##########

### Prepare input for ssGSEA broad gene pattern
expressionTMM.RPKM.GSEA.Input <- expressionTMM.RPKM[, -c(1:7)]; rownames(expressionTMM.RPKM.GSEA.Input) <- expressionTMM.RPKM[,6]
expressionTMM.RPKM.GSEA.print = corUtilsFuncs$createBroadGCTFile(expressionTMM.RPKM.GSEA.Input)

## Only For TCGA+Khanlab dataSet 
# khanlab.TCGA.geneExp <- readRDS("../RNASeq.RSEM/GeneRDSOutput/RPKM/RPKM_Data_Filt_Consolidated.GeneNames.all.TCGA.Khanlab.pc.log22019-03-19.rds")
# expressionTMM.RPKM.GSEA.Input <- khanlab.TCGA.geneExp[, -c(1:7)]; rownames(expressionTMM.RPKM.GSEA.Input) <- khanlab.TCGA.geneExp[,6]

# ## Save input for ssGSEA 
# write.table(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$gseaDir,
#                                       paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.",rnaseqProject$date,".txt"),sep="/"),
#                                       sep="\t", row.names = FALSE, quote = FALSE)
# saveRDS(expressionTMM.RPKM.GSEA.print, paste(rnaseqProject$workDir,rnaseqProject$projectName,rnaseqProject$gseaDir,
#                                       paste0("RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.",rnaseqProject$date,".rds"),sep="/"))

## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/New analysis/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.gct")
#ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
#ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM.TMM.RPKM.GSEA.Input.All.TCGA.Khanlab.2019-03-19.PROJ.gct")

## Add custom expression like cytolytic scre and HLA gene expression to the ssGSEA Outpuut file.
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)
HLA_cytolyticScore      <- rbind(expressionTMM.RPKM.GSEA.Input[c("HLA-A", "HLA-B", "HLA-C"),], cytolyticScore)
View(data.frame(colnames(HLA_cytolyticScore), colnames(ssGSEAScores)))
ssGSEAScores.HLA.Cyto   <- rbind(ssGSEAScores,HLA_cytolyticScore)

## Plot the one variable plot
## Sanity Check: Checking metadata vs data ##
stopifnot( ncol(ssGSEAScores.HLA.Cyto) == length(as.character(rnaseqProject$validMetaDataDF$Sample.Biowulf.ID.GeneExp)) )

## Filter specified Diagnosis
factorsToExclude              = paste(c("NS.", "YST", "Teratoma"), collapse = "|")
selected.metadata              <- rnaseqProject$validMetaDataDF  %>% 
                                  filter_(  .dots = paste0("!grepl(", "'", factorsToExclude , "'" ,",", rnaseqProject$factorName, ")")) %>% 
                                  dplyr::select_( .dots=c(rnaseqProject$metadataFileRefCol, rnaseqProject$factorName ) )

ssGSEAScores.HLA.Cyto.Selected <- ssGSEAScores.HLA.Cyto %>% dplyr::select(one_of(as.character(selected.metadata[, rnaseqProject$metadataFileRefCol])))
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
customColorDF    <- rnaseqProject$customColorsDFAll

### Filter the score matrix for diffent categories in the plot

### For all nothing to be changed
### For only khanlab filter the tidy score matrix & color matrix

Scores <- Scores %>% filter(!grepl("TCGA.",Diagnosis)) %>% filter(!grepl(".CellLine",Diagnosis))
customColorDF <- customColorDF %>% filter(!grepl("TCGA.",Diagnosis)) %>% filter(!grepl(".CellLine",Diagnosis))

## Plot the onevariable plot
plotLists        <- corUtilsFuncs$OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =FALSE, logit =FALSE, plotType = "density",
                                                      yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF )
## Save the plots
EnrischmentScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName                <- paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir,"TMM-RPKM.ssGSEA.enrichmentScores.all.pc.Khanlab.log.pdf",sep="/")
ggsave(SBName, marrangeGrob(EnrischmentScorePlots,ncol=2,nrow=1 ), width = 20, height = 15 )

### For Everything except Cellline

Scores <- Scores %>%  filter(!grepl(".CellLine",Diagnosis))
customColorDF <- customColorDF  %>% filter(!grepl(".CellLine",Diagnosis))

## Plot the onevariable plot
plotLists        <- corUtilsFuncs$OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =TRUE, logit =FALSE, plotType = "density",
                                                      yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF )
## Save the plots
EnrischmentScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName                <- paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir,"TMM-RPKM.ssGSEA.enrichmentScores.all.pc.Khanlab.no.Celline.log.pdf",sep="/")
ggsave(SBName, marrangeGrob(EnrischmentScorePlots,ncol=2,nrow=1 ), width = 20, height = 15 )

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

## Order of genesets
genesets <- c("ImmuneSignature",
              "StromalSignature",
              "Antigen_processing_and_presentation",
              "T.cells_CD8",
              "T.cells_CD4_naive",
              "T.cells_CD4_memory_resting",
              "T.cells_CD4_memory_activated",
              "T.cells_follicular_helper",
              "T.cells_regulatory",
              "T.cells_gamma_delta",
              "NK.cells_activated",
              "NK.cells_resting",
              "B.cells_naive",
              "B.cells_memory",
              "Plasma_cells",
              "Monocytes",
              "Dendritic_cells_resting",
              "M1 Macrophages",
              "M2 Macrophages",
              "Neutrophils",
              "Eosinophils",
              "Mast_cells_resting",
              "Mast_cells_activated"
)

Diagnosis <- c("WT", "SS", "CCSK", "EWS",  "RMS.FN", "RMS.FP", "NB.MYCN.A","NB.Unknown", "DSRCT", "NB.MYCN.NA",  "OS", "UDS", 
                "ML", "HBL", "ASPS")
ScoresForGather$GeneSet <- factor(ScoresForGather$GeneSet, levels = genesets, ordered = TRUE)
ScoresForGather$Diagnosis <- factor(ScoresForGather$Diagnosis, levels = Diagnosis, ordered = TRUE)

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
#write.table(ScoresForSpread, "T:/Sivasish_Sindiri/R Scribble/RNASeq/PlotData/ScoresperDiagSigEnrichZscore.txt", sep="\t", quote = F, col.names = T, row.names = F)

## Open the PDF File
#pdf( paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir,"PercentSamplesEnrichmentScoreHeatMap.Cibersort.pdf", sep="/"), height=10, width = 20)

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

fheatmap(t(ScoresForSpreadHeat), display_tree_col = F,cluster_rows = F, cluster_cols = F, mat_color = matrix_color_vector,
         row_fontsize = 5, col_fontsize = 5, cell_border = F, cell_border_col = "#A6ACAF",seed = 10,
         title = "Percent samples enriched across cancer types")

dev.off()

### Perform correlation analysis  for DDR genes ####
library("Hmisc")
library("psych")
## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
## Read the DDR genes file
DDRGenes <- read.table("T:/Sivasish_Sindiri/R Scribble/Annotation RDS/276_DDR_Genes.v2.txt", sep="\t", header = T, stringsAsFactors = FALSE)
## sanity check if any DDR genes are missing form the expression matrix
excludeGenesIndx <- which(! DDRGenes$DDRGenes %in% rownames(expressionTMM.RPKM.GSEA.Input))
DDRGenesPresent <- DDRGenes$DDRGenes[-excludeGenesIndx]; length(DDRGenesPresent)
## Get gene expression the present genes
DDRGenesGeneExp     <- expressionTMM.RPKM.GSEA.Input[DDRGenesPresent,]; dim(DDRGenesGeneExp)
## Sanity check for NA or Inf
indx <- apply(DDRGenesGeneExp, 1, function(x) any(is.na(x) | is.infinite(x)))
rownames(DDRGenesGeneExp)[indx];
DDRGenesGeneExpFinal <-  DDRGenesGeneExp[complete.cases(DDRGenesGeneExp), ]; dim(DDRGenesGeneExpFinal)
## bind geneexpression matrix to immunescore matrix and perform correlation
ssGSEAScores.DDRGenes  <- rbind(ssGSEAScores,DDRGenesGeneExpFinal)
ssGSEAScores.DDRGenes.corr <- rcorr(t(ssGSEAScores.DDRGenes),  type = "spearman");View(ssGSEAScores.DDRGenes.corr)
## Plot the heatmap
col<- colorRampPalette(c("blue", "white", "red"))(20)
pdf("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/Figures/ImmuneScoreVSDDRgeneExp.pdf", height = 55, width = 55)
heatmap(x = ssGSEAScores.DDRGenes.corr$r, col = col, symm = TRUE, keep.dendro = FALSE)
dev.off()
## slicing out the interesting part
ssGSEAScores.DDRGenes  <- rbind(ssGSEAScores[c("T.regulatory_PMID_30127393_neg", "T.regulatory_PMID_30127393_pos"),],DDRGenesGeneExpFinal)
ssGSEAScores.DDRGenes.corr <- rcorr(t(ssGSEAScores.DDRGenes),  type = "spearman"); View(ssGSEAScores.DDRGenes.corr)
pdf("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/Figures/ImmuneScoreVSDDRgeneExpSelected.pdf", height = 55, width = 55)
heatmap(x = ssGSEAScores.DDRGenes.corr$r, col = col, symm = TRUE, keep.dendro = FALSE)
dev.off()
## Writing data to files
write.table(ssGSEAScores.DDRGenes, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/FigureData/ssGSEAScores.DDRGenes.txt", quote = FALSE, sep = "\t")
write.table(ssGSEAScores.DDRGenes.corr$p, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/FigureData/ImmuneSoreVSDDRgexp.p.txt", quote = FALSE)

### Perform correlation analysis  for CTA genes ####
## Read the ssGSEA output
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.kegg.immune.rm.dups.gct")
## Read the DDR genes file
rm(CGAGenes)
CGAGenes <- as.character(rnaseqProject$csDF$GeneName)
## sanity check if any DDR genes are missing form the expression matrix
excludeGenesIndx <- which(! CGAGenes %in% rownames(expressionTMM.RPKM.GSEA.Input)); length(excludeGenesIndx)
DDRGenesabsent <- CGAGenes[excludeGenesIndx]; as.character(DDRGenesabsent)
DDRGenesPresent <- CGAGenes[-excludeGenesIndx]; length(DDRGenesPresent)
## Get gene expression the present genes
DDRGenesGeneExp     <- expressionTMM.RPKM.GSEA.Input[DDRGenesPresent,]; dim(DDRGenesGeneExp)
## Sanity check for NA or Inf
indx <- apply(DDRGenesGeneExp, 1, function(x) any(is.na(x) | is.infinite(x)))
rownames(DDRGenesGeneExp)[indx];
DDRGenesGeneExpFinal <-  DDRGenesGeneExp[complete.cases(DDRGenesGeneExp), ]; dim(DDRGenesGeneExpFinal)
## bind geneexpression matrix to immunescore matrix and perform correlation
cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)
ssGSEAScores.CGAGenes  <- rbind(ssGSEAScores[4,],cytolyticScore,DDRGenesGeneExpFinal) %>% data.frame()

## Remove Normal samples
Tumor_samples_annot <- AliasNames_df %>% dplyr::filter(!grepl('Teratoma',DIAGNOSIS.Substatus.Tumor.Normal.Tissue))
Tumor_samples_annot <- Tumor_samples_annot[complete.cases(Tumor_samples_annot),]
ssGSEAScores.CGAGenes.Tumor <- ssGSEAScores.CGAGenes %>% dplyr::select(one_of(as.character(Tumor_samples_annot$Sample.Biowulf.ID.GeneExp)))
dim(ssGSEAScores.CGAGenes.Tumor)

## For all tumor samples
# ssGSEAScores.CGAGenes.corr <- rcorr(t(ssGSEAScores.CGAGenes.Tumor), type = "spearman");
# ssGSEAScores.CGAGenes.cytolytic <- flattenCorrMatrix(ssGSEAScores.CGAGenes.corr$r, ssGSEAScores.CGAGenes.corr$P) %>% filter(grepl('CytolyticScore|T.cells_CD8', row))
# write.table(ssGSEAScores.CGAGenes.cytolytic, 
#             "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/CD8.cytolytic.tf.all.Tumor.corr.txt", quote = FALSE, sep="\t")

Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue <- gsub('NS.*','NS',Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
## Perform rcorr for each tumor groups 
correlationDF <- rcorr_groups(Tumor_samples_annot, ssGSEAScores.CGAGenes.Tumor, rnaseqProject$factorName, "Sample.Biowulf.ID.GeneExp")
correlationDF_cytolytic <- correlationDF %>% filter(grepl('CytolyticScore|T.cells_CD8', row))
write.table(correlationDF_cytolytic, 
            "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/CD8.cytolytic.cs.each.Tumor.corr.txt", quote = FALSE, sep="\t")


### Perform annotation for the above
table_diff_Exp_annot <- read.table("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/CD8/CD8.annotation/TranscriptionFactor.Summarised.traditionalRank.Dexp.txt",
                                   sep="\t", header = T, row.names = 1)
colnames(table_diff_Exp_annot) <- gsub("NormalsStatus", '', colnames(table_diff_Exp_annot))

cor.data <- read.table("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/CD8/CD8.annotation/CD8.cytolytic.tf.each.Tumor.corr.txt",
                       sep="\t", header = T, row.names = 1)

diff_exp_annotate <- function(x, lookupDF=NA){
  value = lookupDF[x["column"],x["group"]]
  if(is.null(value)){
    value= NA
  }
  x$DiffExp <- value
  return(data.frame(x))
}

annotatedList <- apply(cor.data, 1, diff_exp_annotate, lookupDF=table_diff_Exp_annot)
annotatedDF <- do.call(rbind.data.frame,annotatedList)
write.table(annotatedDF, 
            "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/CD8.cytolytic.tf.each.Tumor.corr.diff.exp.txt", quote = FALSE, sep="\t")



#### Perform correlation analysis between Exhaustion markers and immune signatures ####
## Read the ssGSEA output
#ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.kegg.immune.rm.dups.gct")
## Read the DDR genes file
EMGenes <- read.table("T:/Sivasish_Sindiri/R Scribble/Annotation RDS/Exhaustion_markers_genes.txt", sep="\t", header = T, stringsAsFactors = FALSE)
## sanity check if any EM genes are missing form the expression matrix
excludeGenesIndx <- which(!EMGenes$Exhaustion_Marker_Genes %in% rownames(expressionTMM.RPKM.GSEA.Input))
if( length(EMGenesPresent) != 0) {
  print("Some Genes are not found !!")
  EMGenesPresent <- EMGenes$Exhaustion_Marker_Genes[-excludeGenesIndx]; length(EMGenesPresent)
} else {
  print("All Genes are found !!")
  EMGenesPresent <- EMGenes$Exhaustion_Marker_Genes; length(EMGenesPresent)
}
## Get gene expression the present genes
EMGenesGeneExp     <- expressionTMM.RPKM.GSEA.Input[EMGenesPresent,]; dim(EMGenesGeneExp)
## Sanity check for NA or Inf
index <- apply(EMGenesGeneExp, 1, function(x) any(is.na(x) | is.infinite(x)))
rownames(EMGenesGeneExp)[index];
EMGenesGeneExpFinal <-  EMGenesGeneExp[complete.cases(EMGenesGeneExp), ]; dim(EMGenesGeneExpFinal)
## bind geneexpression matrix to immunescore matrix and perform correlation
ssGSEAScores.EMGenes  <- rbind(ssGSEAScores, EMGenesGeneExpFinal); dim(ssGSEAScores.EMGenes)

# ## For making heatmap
# ssGSEAScores.EMGenes  <- rbind(ssGSEAScores[c("T.cells_CD4_memory_activated", "T.cells_CD8"),], EMGenesGeneExpFinal); dim(ssGSEAScores.EMGenes)

## Get correlation between Gene expression and Immune signature irrespective of diagnosis
ssGSEAScores.EMGenes.T <- t(ssGSEAScores.EMGenes)
ssGSEAScores.EMGenes.corr <- rcorr(ssGSEAScores.EMGenes.T,  type = "spearman");View(ssGSEAScores.EMGenes.corr)

## Get correlation Immune gene signatures irrespective of diagnosis
library(psych)
library(Hmisc)
ssGSEAScores.T <- t(ssGSEAScores)
ssGSEAScores.corr <- rcorr(ssGSEAScores.T,  type = "spearman");View(ssGSEAScores.corr$r)

## coorelation plot
ssGSEAScores.T.df <- data.frame(ssGSEAScores.T)
ssGSEAScores.T.df[,"diff"] <- abs(ssGSEAScores.T[,"Kegg_Antigen_processing_and_presentation"]-ssGSEAScores.T[,"ImmuneSignature"])
ggplot(ssGSEAScores.T.df, aes(x = Kegg_Antigen_processing_and_presentation, y = ImmuneSignature), color="darkblue") +
  geom_point() +
  theme_bw() + geom_smooth(method = "lm")

## Plot heatmap for the above
col<- colorRampPalette(c("blue", "white", "red"))(100)
pdf("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/Figures/ImmuneScoreVSDDRgeneExp.pdf", height = 55, width = 55)
heatmap(x = ssGSEAScores.EMGenes.corr$r, col = col, symm = TRUE, keep.dendro = FALSE)
pheatmap(ssGSEAScores.corr$r, col = col)
dev.off()

## Add annotation
ssGSEAScores.EMGenes.T <- t(ssGSEAScores.EMGenes)
ssGSEAScores.EMGenes.T <- data.frame(ssGSEAScores.EMGenes.T) %>%  tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")
ssGSEAScores.EMGenes.T.annot.all <- dplyr::left_join(ssGSEAScores.EMGenes.T, AliasNames_df[,c(1,3)], by="Sample.Biowulf.ID.GeneExp") %>%
                                dplyr::rename(Diagnosis=DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
ssGSEAScores.EMGenes.T.annot <- ssGSEAScores.EMGenes.T.annot.all %>% filter( !grepl("NS", Diagnosis) )  %>%
                                tibble::column_to_rownames(var="Sample.Biowulf.ID.GeneExp")
## Make Correlation matrix
rownames <- rep( EMGenes$Exhaustion_Marker_Genes , length( unique(ssGSEAScores.EMGenes.T.annot$Diagnosis)) ); length(rownames)
gp = dplyr::group_by(ssGSEAScores.EMGenes.T.annot, Diagnosis)
gpTowrite <- gp ; rownames(gpTowrite) <- rownames(ssGSEAScores.EMGenes.T.annot)

CorrDF.Out.R <- dplyr::do(gp, data.frame(Cor=t(corr.test(.[,1:43], .[,44:65], method = "spearman")$r))) %>% data.frame() %>%  
                mutate_all( funs_( interp( ~replace(., is.na(.),0) ) ) ); dim(CorrDF.Out.R); head(CorrDF.Out.R)
CorrDF.Out.R <- tibble::add_column(CorrDF.Out.R, GeneNames= rownames, .after=1) 
CorrDF.Out.R$GeneNames <- factor(CorrDF.Out.R$GeneNames, levels = sort(unique(CorrDF.Out.R$GeneNames)))
#write.table(CorrDF.Out.R, paste("./PlotData/", date, "CorrDF.Out.R.EM.spearman.txt", sep=""), sep="\t", row.names = F, quote = FALSE )
CorrDF.Out.P <- dplyr::do(gp, data.frame(Cor=t(corr.test(.[,1:24], .[,44:65], method = "spearman")$p))) %>% data.frame()
CorrDF.Out.P <- tibble::add_column(CorrDF.Out.P, GeneNames= rownames, .after=1)
#write.table(CorrDF.Out.P, paste("./PlotData/", date, "CorrDF.Out.P.EM.spearnman.txt", sep=""), sep="\t", row.names = F, quote = FALSE )      

## Post Analysis; Merge coefficient and p values
CorrDF.Out.R.TC8.TC4 <- CorrDF.Out.R %>% dplyr::select(matches("Diagnosis|GeneNames|CD8|CD4"))
CorrDF.Out.P.TC8.TC4 <- CorrDF.Out.P %>% dplyr::select(matches("Diagnosis|GeneNames|CD8|CD4"))      
CorrDF.Out.R.P.TC8.TC4 <- merge(CorrDF.Out.R.TC8.TC4, CorrDF.Out.P.TC8.TC4, by=c("Diagnosis", "GeneNames"))
CorrDF.Out.R.P.TC8.TC4$Diagnosis <- factor(as.character(CorrDF.Out.R.P.TC8.TC4$Diagnosis),ordered = TRUE,
                                           #levels = sort(unique( CorrDF.Out.R.P.TC8.TC4$Diagnosis )) )
                                           levels = c("NB.MYCN.NA","NB.MYCN.A","EWS","DSRCT","OS","RMS.FP","RMS.FN",
                                                      "Teratoma","SS","CCSK","NB.Unknown","ASPS","HBL","WT","ML","UDS","YST") )
## Filter correlation based on the filter
CorrDF.Out.R.P.CD8.CD4_MemoryAct <- CorrDF.Out.R.P.TC8.TC4 %>% 
                                    mutate(TC8.TC4.MemoryAct=ifelse( abs( Cor.T.cells_CD8.x >= 0.3 & Cor.T.cells_CD8.y <= 0.05 ) |
                                                                    abs( Cor.T.cells_CD4_memory_activated.x >= 0.3 & Cor.T.cells_CD4_memory_activated.y <= 0.05 ), 1, 0 ) ) %>%
                                    dplyr::select(Diagnosis, GeneNames, TC8.TC4.MemoryAct) 
CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread <- CorrDF.Out.R.P.CD8.CD4_MemoryAct %>% tidyr::spread(Diagnosis,TC8.TC4.MemoryAct)
CD8.CD4_MemoryAct.Spread <- CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread  %>% dplyr::mutate(Count=rowSums(.[2:ncol(CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread)]))
CD8.CD4_MemoryAct.Spread <- tibble::add_column( CD8.CD4_MemoryAct.Spread, Legend=paste(CD8.CD4_MemoryAct.Spread$GeneNames, 
                                                "(", 
                                                CD8.CD4_MemoryAct.Spread$Count, ")"
                                              ), .after=1) %>% data.frame() %>%  
                                              dplyr::select(-contains("GeneNames"))

## Filter correlation based on the filter plot correlation values only.
CorrDF.Out.R.P.CD8.heatmap <- CorrDF.Out.R.P.TC8.TC4 %>% 
  mutate(TC8.TC4.MemoryAct=ifelse( abs( Cor.T.cells_CD8.y <= 0.05 ), Cor.T.cells_CD8.x, 0 ) ) %>%
  dplyr::select(Diagnosis, GeneNames, TC8.TC4.MemoryAct) 
CorrDF.Out.R.P.CD8.heatmap.Spread <- CorrDF.Out.R.P.CD8.heatmap %>% tidyr::spread(Diagnosis,TC8.TC4.MemoryAct)
CorrDF.Out.R.P.CD8.heatmap.Spread <- CorrDF.Out.R.P.CD8.heatmap.Spread  %>% dplyr::mutate(Count=rowSums(.[2:ncol(CorrDF.Out.R.P.CD8.heatmap.Spread)]))
CorrDF.Out.R.P.CD8.heatmap.Spread <- tibble::add_column( CorrDF.Out.R.P.CD8.heatmap.Spread, Legend=paste(CorrDF.Out.R.P.CD8.heatmap.Spread$GeneNames), .after=1) %>% 
                                                data.frame() %>% dplyr::select(-contains("GeneNames"))

## Plot the heatmap
CorrDF <- CorrDF.Out.R.P.CD8.heatmap.Spread
#CorrDF <- CorrDF.Out.R.P.CD8.CD4_MemoryAct
CorrDF %<>% 
  dplyr::arrange(Count) %<>% 
  dplyr::filter( Count > 0 ) %<>% 
  dplyr::select(-one_of("Count")) %<>% 
  dplyr::select(-matches("ML|UDS|YST")) %<>%
  tibble::column_to_rownames(var="Legend")
# CD8.CD4_MemoryAct.Spread %<>%   dplyr::arrange(desc(Legend)) %<>% tibble::column_to_rownames(var="Legend")
# CD8.CD4_MemoryAct.Spread %<>% column_to_rownames(var="Legend")

pdf(paste("./Plots/", date, "CD8 significant correlation with EM.spearman.TMM.RPKM.GP.log2.Output.pdf", sep=""), height = 10, width = 15)
superheat(CorrDF,
          #title = "CD4-Memory.Activated/CD8 significant correlation with EM",
          title = "CD8 significant correlation with EM",
          legend=FALSE,
          grid.hline = FALSE,
          grid.vline = FALSE,
          pretty.order.rows = FALSE,
          pretty.order.cols = FALSE,
          # grid.hline.size = 0.01,
          # grid.vline.size = 0.01,
          # heat.col.scheme = "grey",
          heat.lim = c(0, 1),
          #heat.pal = c("#004080",  "#88cc00"),
          heat.pal = c("#e0e0d1", "#004080"),
          bottom.label.text.angle=90,
          title.size = 6)
dev.off()


#### Perform correlation analysis between Exhaustion markers and immune signatures ####
## Read the ssGSEA output
#ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.pc.log2.2019-01-31.PROJ.KeggSig.gct")
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.kegg.immune.rm.dups.gct")
## Read the DDR genes file
EMGenes <- read.table("C:/Users/sindiris/R Scribble/Annotation RDS/Exhaustion_markers_genes.txt", sep="\t", header = T, stringsAsFactors = FALSE)
## sanity check if any EM genes are missing form the expression matrix
excludeGenesIndx <- which(!EMGenes$Exhaustion_Marker_Genes %in% rownames(expressionTMM.RPKM.GSEA.Input))
if( length(EMGenesPresent) != 0) {
  print("Some Genes are not found !!")
  EMGenesPresent <- EMGenes$Exhaustion_Marker_Genes[-excludeGenesIndx]; length(EMGenesPresent)
} else {
  print("All Genes are found !!")
  EMGenesPresent <- EMGenes$Exhaustion_Marker_Genes; length(EMGenesPresent)
}
## Get gene expression the present genes
EMGenesGeneExp     <- expressionTMM.RPKM.GSEA.Input[EMGenesPresent,]; dim(EMGenesGeneExp)
## Sanity check for NA or Inf
index <- apply(EMGenesGeneExp, 1, function(x) any(is.na(x) | is.infinite(x)))
rownames(EMGenesGeneExp)[index];
EMGenesGeneExpFinal <-  EMGenesGeneExp[complete.cases(EMGenesGeneExp), ]; dim(EMGenesGeneExpFinal)
## bind geneexpression matrix to immunescore matrix and perform correlation
ssGSEAScores.EMGenes  <- rbind(ssGSEAScores, EMGenesGeneExpFinal); dim(ssGSEAScores.EMGenes)

# ## For making heatmap
# ssGSEAScores.EMGenes  <- rbind(ssGSEAScores[c("T.cells_CD4_memory_activated", "T.cells_CD8"),], EMGenesGeneExpFinal); dim(ssGSEAScores.EMGenes)

## Get correlation between Gene expression and Immune signature irrespective of diagnosis
ssGSEAScores.EMGenes.T <- t(ssGSEAScores.EMGenes)
ssGSEAScores.EMGenes.corr <- rcorr(ssGSEAScores.EMGenes.T,  type = "spearman");View(ssGSEAScores.EMGenes.corr)

## Get correlation Immune gene signatures irrespective of diagnosis
library(psych)
library(Hmisc)
ssGSEAScores.T <- t(ssGSEAScores)
ssGSEAScores.corr <- rcorr(ssGSEAScores.T,  type = "spearman");View(ssGSEAScores.corr$r)

## coorelation plot
ssGSEAScores.T.df <- data.frame(ssGSEAScores.T)
ssGSEAScores.T.df[,"diff"] <- abs(ssGSEAScores.T[,"Kegg_Antigen_processing_and_presentation"]-ssGSEAScores.T[,"ImmuneSignature"])
ggplot(ssGSEAScores.T.df, aes(x = Kegg_Antigen_processing_and_presentation, y = ImmuneSignature), color="darkblue") +
  geom_point() +
  theme_bw() + geom_smooth(method = "lm")

## Plot heatmap for the above
col<- colorRampPalette(c("blue", "white", "red"))(100)
pdf("C:/Users/sindiris/R Scribble/RNASeq.RSEM/Figures/ImmuneScoreVSDDRgeneExp.pdf", height = 55, width = 55)
heatmap(x = ssGSEAScores.EMGenes.corr$r, col = col, symm = TRUE, keep.dendro = FALSE)
pheatmap(ssGSEAScores.corr$r, col = col)
dev.off()

## Add annotation
ssGSEAScores.EMGenes.T <- t(ssGSEAScores.EMGenes)
ssGSEAScores.EMGenes.T <- data.frame(ssGSEAScores.EMGenes.T) %>%  tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")
ssGSEAScores.EMGenes.T.annot.all <- dplyr::left_join(ssGSEAScores.EMGenes.T, AliasNames_df[,c(1,3)], by="Sample.Biowulf.ID.GeneExp") %>%
                                dplyr::rename(Diagnosis=DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
ssGSEAScores.EMGenes.T.annot <- ssGSEAScores.EMGenes.T.annot.all %>% filter( !grepl("NS", Diagnosis) )  %>%
                                tibble::column_to_rownames(var="Sample.Biowulf.ID.GeneExp")
## Make Correlation matrix
rownames <- rep( EMGenes$Exhaustion_Marker_Genes , length( unique(ssGSEAScores.EMGenes.T.annot$Diagnosis)) ); length(rownames)
gp = dplyr::group_by(ssGSEAScores.EMGenes.T.annot, Diagnosis)
gpTowrite <- gp ; rownames(gpTowrite) <- rownames(ssGSEAScores.EMGenes.T.annot)

CorrDF.Out.R <- dplyr::do(gp, data.frame(Cor=t(corr.test(.[,1:43], .[,44:65], method = "spearman")$r))) %>% data.frame() %>%  
                mutate_all( funs_( interp( ~replace(., is.na(.),0) ) ) ); dim(CorrDF.Out.R); head(CorrDF.Out.R)
CorrDF.Out.R <- tibble::add_column(CorrDF.Out.R, GeneNames= rownames, .after=1) 
CorrDF.Out.R$GeneNames <- factor(CorrDF.Out.R$GeneNames, levels = sort(unique(CorrDF.Out.R$GeneNames)))
#write.table(CorrDF.Out.R, paste("./PlotData/", date, "CorrDF.Out.R.EM.spearman.txt", sep=""), sep="\t", row.names = F, quote = FALSE )
CorrDF.Out.P <- dplyr::do(gp, data.frame(Cor=t(corr.test(.[,1:24], .[,44:65], method = "spearman")$p))) %>% data.frame()
CorrDF.Out.P <- tibble::add_column(CorrDF.Out.P, GeneNames= rownames, .after=1)
#write.table(CorrDF.Out.P, paste("./PlotData/", date, "CorrDF.Out.P.EM.spearnman.txt", sep=""), sep="\t", row.names = F, quote = FALSE )      

## Post Analysis; Merge coefficient and p values
CorrDF.Out.R.TC8.TC4 <- CorrDF.Out.R %>% dplyr::select(matches("Diagnosis|GeneNames|CD8|CD4"))
CorrDF.Out.P.TC8.TC4 <- CorrDF.Out.P %>% dplyr::select(matches("Diagnosis|GeneNames|CD8|CD4"))      
CorrDF.Out.R.P.TC8.TC4 <- merge(CorrDF.Out.R.TC8.TC4, CorrDF.Out.P.TC8.TC4, by=c("Diagnosis", "GeneNames"))
CorrDF.Out.R.P.TC8.TC4$Diagnosis <- factor(as.character(CorrDF.Out.R.P.TC8.TC4$Diagnosis),ordered = TRUE,
                                           #levels = sort(unique( CorrDF.Out.R.P.TC8.TC4$Diagnosis )) )
                                           levels = c("NB.MYCN.NA","NB.MYCN.A","EWS","DSRCT","OS","RMS.FP","RMS.FN",
                                                      "Teratoma","SS","CCSK","NB.Unknown","ASPS","HBL","WT","ML","UDS","YST") )
## Filter correlation based on the filter
CorrDF.Out.R.P.CD8.CD4_MemoryAct <- CorrDF.Out.R.P.TC8.TC4 %>% 
                                    mutate(TC8.TC4.MemoryAct=ifelse( abs( Cor.T.cells_CD8.x >= 0.3 & Cor.T.cells_CD8.y <= 0.05 ) |
                                                                    abs( Cor.T.cells_CD4_memory_activated.x >= 0.3 & Cor.T.cells_CD4_memory_activated.y <= 0.05 ), 1, 0 ) ) %>%
                                    dplyr::select(Diagnosis, GeneNames, TC8.TC4.MemoryAct) 
CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread <- CorrDF.Out.R.P.CD8.CD4_MemoryAct %>% tidyr::spread(Diagnosis,TC8.TC4.MemoryAct)
CD8.CD4_MemoryAct.Spread <- CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread  %>% dplyr::mutate(Count=rowSums(.[2:ncol(CorrDF.Out.R.P.CD8.CD4_MemoryAct.Spread)]))
CD8.CD4_MemoryAct.Spread <- tibble::add_column( CD8.CD4_MemoryAct.Spread, Legend=paste(CD8.CD4_MemoryAct.Spread$GeneNames, 
                                                "(", 
                                                CD8.CD4_MemoryAct.Spread$Count, ")"
                                              ), .after=1) %>% data.frame() %>%  
                                              dplyr::select(-contains("GeneNames"))

## Filter correlation based on the filter plot correlation values only.
CorrDF.Out.R.P.CD8.heatmap <- CorrDF.Out.R.P.TC8.TC4 %>% 
  mutate(TC8.TC4.MemoryAct=ifelse( abs( Cor.T.cells_CD8.y <= 0.05 ), Cor.T.cells_CD8.x, 0 ) ) %>%
  dplyr::select(Diagnosis, GeneNames, TC8.TC4.MemoryAct) 
CorrDF.Out.R.P.CD8.heatmap.Spread <- CorrDF.Out.R.P.CD8.heatmap %>% tidyr::spread(Diagnosis,TC8.TC4.MemoryAct)
CorrDF.Out.R.P.CD8.heatmap.Spread <- CorrDF.Out.R.P.CD8.heatmap.Spread  %>% dplyr::mutate(Count=rowSums(.[2:ncol(CorrDF.Out.R.P.CD8.heatmap.Spread)]))
CorrDF.Out.R.P.CD8.heatmap.Spread <- tibble::add_column( CorrDF.Out.R.P.CD8.heatmap.Spread, Legend=paste(CorrDF.Out.R.P.CD8.heatmap.Spread$GeneNames), .after=1) %>% 
                                                data.frame() %>% dplyr::select(-contains("GeneNames"))

## Plot the heatmap
CorrDF <- CorrDF.Out.R.P.CD8.heatmap.Spread
#CorrDF <- CorrDF.Out.R.P.CD8.CD4_MemoryAct
CorrDF %<>% 
  dplyr::arrange(Count) %<>% 
  dplyr::filter( Count > 0 ) %<>% 
  dplyr::select(-one_of("Count")) %<>% 
  dplyr::select(-matches("ML|UDS|YST")) %<>%
  tibble::column_to_rownames(var="Legend")
# CD8.CD4_MemoryAct.Spread %<>%   dplyr::arrange(desc(Legend)) %<>% tibble::column_to_rownames(var="Legend")
# CD8.CD4_MemoryAct.Spread %<>% column_to_rownames(var="Legend")

pdf(paste("./Plots/", date, "CD8 significant correlation with EM.spearman.TMM.RPKM.GP.log2.Output.pdf", sep=""), height = 10, width = 15)
superheat(CorrDF,
          #title = "CD4-Memory.Activated/CD8 significant correlation with EM",
          title = "CD8 significant correlation with EM",
          legend=FALSE,
          grid.hline = FALSE,
          grid.vline = FALSE,
          pretty.order.rows = FALSE,
          pretty.order.cols = FALSE,
          # grid.hline.size = 0.01,
          # grid.vline.size = 0.01,
          # heat.col.scheme = "grey",
          heat.lim = c(0, 1),
          #heat.pal = c("#004080",  "#88cc00"),
          heat.pal = c("#e0e0d1", "#004080"),
          bottom.label.text.angle=90,
          title.size = 6)
dev.off()

colfunc <- colorRampPalette(c("#f2f2f2", "gold","firebrick3"))
pheatmap(CorrDF, 
         #color =c("#e0e0d1", "#004080"), 
         #color =c("#f2f2f2", "#004080"), 
         color = colfunc(10),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = "grey", 
         #border_color = NA,
         treeheight_row = 0,
         fontsize = 12,
         legend = TRUE,
         main="CD8 correlation with EM")

### Perform Differential gene expression analysis ####

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


### Intantiate a new Differential Gene Expression Object ####
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
# 
DiffExpObj <- dgeObj$performDiffGeneExp()
# head(DiffExpObj[[1]] %>% dplyr::arrange(-logFC))
# 

### Filtering ####
# Step 0   Define fucnctions ####

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
group2FPKM.T = 1 ; group1FPKM.T = 1;  PValue.T = 0.05 ; logFoldDiff.T = 0 ; FDR_value.T = 0.05 ; vitalFPKM.T = 0

# selectedGeneList <- "CancerGermlineAntigen"
# Zscored.logFC = 0.25; Zscore.group2 = 0.5; group2FPKM = 5; group1FPKM = 1;  PValue = 0.001; logFC = 4; FDR = 0.05
# 
# selectedGeneList <- "CellSurface"
# Zscored.logFC = 1 ; Zscore.group2 = 0.5; group2FPKM = 5 ; group1FPKM = 1;  PValue = 0.001 ; logFC = 4 ; FDR = 0.05
#
# selectedGeneList <- "TranscriptionFactor"
# Zscored.logFC = 0.25 ; Zscore.group2 = 0.5; group2FPKM = 5 ; group1FPKM = 1;  PValue = 0.001 ; logFC = 4 ; FDR = 0.05
#
selectedGeneList <- "ExhaustionMarkers"
Zscored.logFC = 0.25 ; Zscore.group2 = 0.5; group2FPKM = 5 ; group1FPKM = 1;  PValue = 0.001 ; logFC = 4 ; FDR = 0.05

MergedDiffExpResultDir <- paste0("T:/Sivasish_Sindiri/R Scribble//RNASeq.RSEM/MergedDiffExpResults/",selectedGeneList)

# Step 2.  Perform Merging of differential expression file across groups ####
dir.create(MergedDiffExpResultDir)
ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })), c("Normals_WT", "Normals_CCSK") )
#ConditionGroup <- c(unique(sapply(dgeObj$pairedList, function(x){ return(paste(x[1],x[2],sep = "_"))  })))
groups <- list.dirs(paste("T:/Sivasish_Sindiri/R Scribble//RNASeq.RSEM//DiffExpResults/", sep=""))[-1]; groups[1]
output <- sapply(groups, mergeDiffTestResults, type="Gene", colInterest=c(7,9,10,11,12, 15:28), rowNamesCol = 2,
                 fileSuffix=paste0(selectedGeneList,".txt"), saveDirPath=MergedDiffExpResultDir)

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
  
  ## write.table(tumorAllData, paste(MergedDiffExpResultDir,"/",x,"/",x,".allgenes.DiffExp.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  ## Zscore Ranking filter
  tumorAllData.zscore <- tumorAllData %>% 
    dplyr::filter_(.dots=paste0( 
      groupsCompare[2]," >= ", group2FPKM ,
      " &  Zscored.logFC   >= ", Zscored.logFC,
      " &   logFC >", logFC,
      " & ", paste0("Zscored.",groupsCompare[2]), " >= ", Zscore.group2 )) %>%
    dplyr::arrange_(.dots = paste0("desc(","Zscored.",groupsCompare[2], ")" ) )
  
  ## Complete filtered gene List with zscoring filter
  ## write.table(tumorAllData.zscore, paste(MergedDiffExpResultDir,"/",x,"/",x,".filteredgenes.ZscoringRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  ## Complete filtered gene List (Binary) with zscoring filter
  selectGenes <- tumorAllData.zscore %>%  dplyr::select(GeneName)
  tumorData.zscoreRanking <- tumorDataPvalue
  tumorData.zscoreRanking["status"] <- 0
  statusDF.zscoreRanking <- tumorData.zscoreRanking %>% mutate(status=ifelse(GeneName %in% selectGenes$GeneName, 1, 0)) %>% dplyr::select(GeneName, status) %>%
    rename(c('status'=paste(groupsCompare[2],groupsCompare[1],"Status", sep="")))
  statusDF.zscoreRanking <- statusDF.zscoreRanking[, !duplicated(colnames(statusDF.zscoreRanking))] 
  
  # write.table(statusDF.zscoreRanking, paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.zscoreRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  
  ## Javed's Filtering
  print(colnames(tumorAllData))
  tumorAllData.filt <- tumorAllData %>% dplyr::filter_(.dots=paste0( groupsCompare[2], " >= ", group2FPKM.T ,
                                                                           " &  logFC >=", logFoldDiff.T,
                                                                           " &  PValue   <=", PValue.T ,
                                                                           " &  Brain.MeanExp  < ", vitalFPKM.T ,
                                                                           " &  Heart.MeanExp  < ", vitalFPKM.T   )) %>%
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
  
  write.table(statusDF.traditionalRanking, paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.txt",sep=""), sep="\t", row.names = FALSE, quote = FALSE)
  
  
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
# write.table(allTumorMergedStats[[1]], paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.ZscoreRank.Dexp.txt",  sep=""),
#             sep="\t", row.names = FALSE, quote = FALSE)
write.table(allTumorMergedStats[[2]], paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.Dexp.txt",sep=""),
            sep="\t", row.names = FALSE, quote = FALSE)

# Step 6.  Select rows for heatmap ####
#allTumorStatsFinal <- read.table(paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.ZscoreRank.Dexp.txt",sep=""),sep="\t", header = TRUE)
allTumorStatsFinal <- read.table(paste(MergedDiffExpResultDir,"/",selectedGeneList,".Summarised.traditionalRank.Dexp.txt",sep=""),sep="\t", header = TRUE)
CTA.Filt <- allTumorStatsFinal %>% filter(RowSum >= 1) %>% 
  dplyr::arrange(-RowSum) %>% t() %>% data.frame()
colnames(CTA.Filt) <- as.character(unlist(CTA.Filt[c("GeneName"),]))
CTA.Filt.sorted <- CTA.Filt[-1,]
CTA.Filt.sorted <- CTA.Filt.sorted %>% tibble::rownames_to_column("Diagnosis")
CTA.Filt.sorted <- CTA.Filt.sorted %>% dplyr::arrange_(.dots = list(paste0("desc(",colnames(CTA.Filt.sorted)[2], ")")))
rownames(CTA.Filt.sorted) <- CTA.Filt.sorted[,1]
CTA.Filt.sorted <- CTA.Filt.sorted[-1,]
CTA.Filt.sorted <- CTA.Filt.sorted[,-1]
#filter(, GeneName %in% c("CD99", "FGFR4", "ALK", "GPC2", "MYCN", "MYOG", "MYOD1", "IGF2", "CTAG1B"))
View(CTA.Filt.sorted);dim(CTA.Filt.sorted)

# Step 7.  Plot the heatmap ####
rownames(CTA.Filt.sorted) <- gsub("NormalsStatus", "", rownames(CTA.Filt.sorted))
indx <- sapply(CTA.Filt.sorted, is.factor)
CTA.Filt.sorted[indx] <- lapply(CTA.Filt.sorted[indx], function(x) as.numeric(as.character(x)))

pdf( paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir, 
           paste0("Differentially Expressed .", selectedGeneList, ".v19.pdf"),  sep="/"), height = 10, width = 25)
pheatmap(CTA.Filt.sorted, 
         #color =c("#e0e0d1", "#004080"), 
         color =c("#f2f2f2", "#004080"), 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         #border_color = "grey", 
         border_color = NA,
         treeheight_row = 0,
         fontsize = 12,
         legend = FALSE )
dev.off()

## Alternate heatmap
#pdf( paste(rnaseqProject$workDir, rnaseqProject$projectName, rnaseqProject$plotsDir, "zscore.Differentially Expressed TF.v18.pdf", sep="/"), height = 10, width = 25)
superheat(CTA.Filt.sorted, pretty.order.cols =F,pretty.order.rows=F,
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
#dev.off()


### Performing TCR analysis

#### Perfroming TCR analysis ####

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

### Start analysis for clone type: Choose clone type ####
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
#### Adding Jun's Color Scheme data obtained from TCRv2.R script
### Unfortunately, change of color decision 
### was taken at the end of project, and its very difficult me to change multiple things 
### Team decided to keep the previous color scheme as it is.
customColorsVector <- data.frame(Color=unique(as.character(toPlotDF$Color.Jun)), Diagnosis= unique(as.character(toPlotDF$Diagnosis)) )
#### 
varNames <- colnames(entropyMetassGSEA[,16:58]) 
plotLists <- lapply(varNames, correlationPlots, constName="Htot..Entropy.", xlab="Entropy", df= data.frame(entropyMetassGSEA), 
                    customColorDF=customColorsVector)
ImmuneScorePlots <- lapply(plotLists, function(l) l[[1]] )
SBName =paste0(TCRResultsDir,"ImmuneScore.vs.Htot..Entropy",cloneType,"v2.jun.pdf")
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



### Neoantigen analysis

#### Neoantigen Post-Processing #####

## neoantigen from variants

neoantigenFromVariants <- read.csv("T:/Sivasish_Sindiri/R Scribble/RNASeq.Mutation.data/NeoantigenCountFromVariants.txt", sep = "\t", header = T) %>% data.table()
neoantigenFromVariantsAnnot <- dplyr::full_join( rnaseqProject$metaDataDF, neoantigenFromVariants, by="Sample.Biowulf.ID") %>% 
                                dplyr::filter(!is.na(Patient.ID)) %>%
                                dplyr::filter( ! LIBRARY_TYPE %in% c("CellLine","Normal")) %>%
                                dplyr::mutate(VariantNeoAntigenCount = ifelse(is.na(VariantNeoAntigenCount),0,VariantNeoAntigenCount)) %>%
                                dplyr::select_(.dots=c("Sample.Biowulf.ID","DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Sample.ID.Alias", "VariantNeoAntigenCount"))
dim(neoantigenFromVariantsAnnot); 
#View(neoantigenFromVariantsAnnot)
## neoantigen from fusions
neoantigenFromFusions <- read.csv("T:/Sivasish_Sindiri/R Scribble/RNASeq.Mutation.data/NeoantigenCountFromFusions.txt", sep = "\t", header = T) %>% data.table() %>% 
                               dplyr::select_(.dots=c("Sample.Biowulf.ID","DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Sample.ID.Alias", "FusionNeoAntigenCount"))
dim(neoantigenFromFusions); 
#View(neoantigenFromFusions)                                

##Merge the above two tables
neoantigenFromSamples <- dplyr::full_join( neoantigenFromVariantsAnnot, neoantigenFromFusions,  
                                                          by=c("Sample.Biowulf.ID","DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Sample.ID.Alias") ); 
neoantigenFromSamplesFinal <- neoantigenFromSamples %>% 
                              dplyr::mutate(FusionNeoAntigenCount = ifelse(is.na(FusionNeoAntigenCount),0,FusionNeoAntigenCount)) %>% 
                              dplyr::mutate(VariantNeoAntigenCount = ifelse(is.na(VariantNeoAntigenCount),0,VariantNeoAntigenCount)) %>% 
                              dplyr::mutate(TotalNeoantigenCount = FusionNeoAntigenCount +  VariantNeoAntigenCount ) %>% 
                              data.table()
dim(neoantigenFromSamplesFinal);
#View(neoantigenFromSamplesFinal)    

## Bean Plot
neoantigenBurden <- neoantigenFromSamplesFinal[,c("TotalNeoantigenCount", "DIAGNOSIS.Substatus.Tumor.Normal.Tissue")] %>% 
                                                                  dplyr::rename(Diagnosis=DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
orderOfFactor           <- as.character( unique(neoantigenBurden$Diagnosis) )
orderOfSignature        <- colnames(neoantigenBurden)[-ncol(neoantigenBurden)]
colList                 <- c(1:(ncol(neoantigenBurden)-1)) ; Scores <- neoantigenBurden
## Generate custom colors
customColorDF    <- rnaseqProject$customColorsDF

## Filter for diagnosis
Scores <- neoantigenBurden %>% filter(!Diagnosis %in% c("Teratoma", "YST")) %>% dplyr::filter(complete.cases(.))

### Plot and Save ###
plotLists <- corUtilsFuncs$OneVariablePlotSort( colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =FALSE, logit =TRUE, logBase=10,
                                                yLab = "log10( NeoantigenBurden )", legendDisplay = FALSE, customColorDF = customColorDF, 
                                                plotType = "StringBean", sizeOfDots = 0.8)

### Get quartiles
## summary of TotalNeoantigenCount
summary_TotalNeoantigenCount <- summary(neoantigenFromSamplesFinal$TotalNeoantigenCount)
FirststQu <- as.numeric(summary_TotalNeoantigenCount["1st Qu."]); paste("FirststQu:", FirststQu)
Median    <- as.numeric(summary_TotalNeoantigenCount["Median"]); paste("Median ", Median)
thirdQu   <- as.numeric(summary_TotalNeoantigenCount["3rd Qu."]); paste("thirdQu ", thirdQu)

# neoantigenFromSamplesFinal[TotalNeoantigenCount <= FirststQu, Sample.ID.Alias ] %>% length()
# neoantigenFromSamplesFinal[TotalNeoantigenCount > FirststQu & TotalNeoantigenCount < Median, Sample.ID.Alias ] %>% length()
neoantigenFromSamplesFinal[TotalNeoantigenCount < Median, Sample.ID.Alias ] %>% length()
neoantigenFromSamplesFinal[TotalNeoantigenCount >= Median, Sample.ID.Alias ] %>% length()

#### Split data matrix into High, Intermediate and Low expression matrices
### Expression Matrix
expressionTMM.RPKM.Neoantigen <- expressionTMM.RPKM %>% dplyr::select(-one_of(c("Chr","Start","End","GeneID", "Length", "Strand")))

## Following Samples removed because neoantigen burden for variant calling failed.
samplesToRemove <- c("OS.PALWWX", "OS.PANPUM", "OS.PAUUML", "OS.PAUYTT", "OS.PALZGU", "OS.PAMHYN", "OS.PANSEN")
expressionTMM.RPKM.Neoantigen.SR <- expressionTMM.RPKM.Neoantigen %>% dplyr::select(-one_of(samplesToRemove))
dim(expressionTMM.RPKM.Neoantigen.SR)
### Get Mean expression of each category ####

############################################### Neoantigen analysis for the cohort #######################################
### Get mean expression of low neoantgen burden samples.
FPKM.Data.NeoAntiBurd.Low    <- expressionTMM.RPKM.Neoantigen.SR %>% dplyr::select( one_of( 
  c("GeneName",as.character(neoantigenFromSamplesFinal[TotalNeoantigenCount <= FirststQu, Sample.ID.Alias ]) ) ))
FPKM.Data.NeoAntiBurd.Low %<>% tibble::column_to_rownames(var="GeneName")
dim(FPKM.Data.NeoAntiBurd.Low); FPKM.Data.NeoAntiBurd.Low[1:5,1:5]

#FPKM.Data.NeoAntiBurd.Low.mean     <- as.data.frame(apply(log2(FPKM.Data.NeoAntiBurd.Low + 1), 1, mean)); 
FPKM.Data.NeoAntiBurd.Low.mean     <- as.data.frame(apply(FPKM.Data.NeoAntiBurd.Low, 1, mean)); 
colnames(FPKM.Data.NeoAntiBurd.Low.mean) <- "FPKM.Data.NeoAntiBurd.Low.mean"
head(FPKM.Data.NeoAntiBurd.Low.mean); dim(FPKM.Data.NeoAntiBurd.Low);

### Get mean expression of intermediate neoantgen burden samples.
FPKM.Data.NeoAntiBurd.Intermediate    <- expressionTMM.RPKM.Neoantigen.SR %>% dplyr::select( one_of( 
  c("GeneName",as.character(neoantigenFromSamplesFinal[TotalNeoantigenCount > FirststQu & TotalNeoantigenCount < Median, Sample.ID.Alias ]) ) ))
FPKM.Data.NeoAntiBurd.Intermediate %<>% tibble::column_to_rownames(var="GeneName")
dim(FPKM.Data.NeoAntiBurd.Intermediate); FPKM.Data.NeoAntiBurd.Intermediate[1:5,1:5]

############################################ New Addition Start ###############

### Get mean expression of intermediate neoantgen burden samples.
FPKM.Data.NeoAntiBurd.Intermediate    <- expressionTMM.RPKM.Neoantigen.SR %>% dplyr::select( one_of( 
  c("GeneName",as.character(neoantigenFromSamplesFinal[TotalNeoantigenCount < Median, Sample.ID.Alias ]) ) ))
FPKM.Data.NeoAntiBurd.Intermediate %<>% tibble::column_to_rownames(var="GeneName")
dim(FPKM.Data.NeoAntiBurd.Intermediate); FPKM.Data.NeoAntiBurd.Intermediate[1:5,1:5]

#FPKM.Data.NeoAntiBurd.Intermediate.mean     <- as.data.frame(apply(log2(FPKM.Data.NeoAntiBurd.Intermediate + 1), 1, mean)); 
FPKM.Data.NeoAntiBurd.Intermediate.mean     <- as.data.frame(apply(FPKM.Data.NeoAntiBurd.Intermediate, 1, mean))
colnames(FPKM.Data.NeoAntiBurd.Intermediate.mean) <- "FPKM.Data.NeoAntiBurd.Intermediate.mean"
head(FPKM.Data.NeoAntiBurd.Intermediate.mean); dim(FPKM.Data.NeoAntiBurd.Intermediate);

############################################ New Addition End #################

### Get mean expression of high neoantgen burden samples.
FPKM.Data.NeoAntiBurd.High    <- expressionTMM.RPKM.Neoantigen.SR %>% dplyr::select( one_of( 
  c("GeneName",as.character(neoantigenFromSamplesFinal[TotalNeoantigenCount >= Median, Sample.ID.Alias ]) ) ))
FPKM.Data.NeoAntiBurd.High %<>% tibble::column_to_rownames(var="GeneName")
dim(FPKM.Data.NeoAntiBurd.High); FPKM.Data.NeoAntiBurd.High[1:5,1:5]

#FPKM.Data.NeoAntiBurd.High.mean     <- as.data.frame(apply(log2(FPKM.Data.NeoAntiBurd.High + 1), 1, mean)); 
FPKM.Data.NeoAntiBurd.High.mean     <- as.data.frame(apply(FPKM.Data.NeoAntiBurd.High, 1, mean)); 
colnames(FPKM.Data.NeoAntiBurd.High.mean) <- "FPKM.Data.NeoAntiBurd.High.mean"
head(FPKM.Data.NeoAntiBurd.High.mean); dim(FPKM.Data.NeoAntiBurd.High);

### Get expression of 2 categories for comparision ####
IntermediateLow <- cbind(FPKM.Data.NeoAntiBurd.Intermediate, FPKM.Data.NeoAntiBurd.Low)
#IntermediateLow.Mean <- as.data.frame(apply(log2(IntermediateLow+1), 1, mean)); colnames(IntermediateLow.Mean) <- "IntermediateLow"
IntermediateLow.Mean <- as.data.frame(apply(IntermediateLow, 1, mean)); colnames(IntermediateLow.Mean) <- "IntermediateLow"
dim(IntermediateLow.Mean); head(IntermediateLow.Mean)

HighLow <- cbind(FPKM.Data.NeoAntiBurd.High, FPKM.Data.NeoAntiBurd.Low)
#HighLow.Mean <- as.data.frame(apply(log2(HighLow+1), 1, mean)); colnames(HighLow.Mean) <- "HighLow"
HighLow.Mean <- as.data.frame(apply(HighLow, 1, mean)); colnames(HighLow.Mean) <- "HighLow"
dim(HighLow.Mean); head(HighLow.Mean)

HighIntermediate <- cbind(FPKM.Data.NeoAntiBurd.High, FPKM.Data.NeoAntiBurd.Intermediate)
#HighIntermediate.Mean <- as.data.frame(apply(log2(HighIntermediate+1), 1, mean)); colnames(HighIntermediate.Mean) <- "HighIntermediate"
HighIntermediate.Mean <- as.data.frame(apply(HighIntermediate, 1, mean)); colnames(HighIntermediate.Mean) <- "HighIntermediate"
dim(HighIntermediate.Mean); head(HighIntermediate.Mean)

### Find ratio between one vs rest of the groups ####
FC.High.IntermediateLow <-  log2(FPKM.Data.NeoAntiBurd.High.mean + 1) - log2(IntermediateLow.Mean + 1)
FC.High.IntermediateLow <- cbind(expressionTMM.RPKM[,c("Chr","Start","End","GeneName","GeneID")], FC.High.IntermediateLow)
FC.High.IntermediateLow.rnk <- FC.High.IntermediateLow %>% dplyr::arrange(-FPKM.Data.NeoAntiBurd.High.mean)
head(FC.High.IntermediateLow.rnk); dim(FC.High.IntermediateLow.rnk)
write.table(FC.High.IntermediateLow.rnk[,c(4,6)], paste0("../RNASeq.RSEM/GSEA/rnk/FC.High.IntermediateLow.pc.v6.rnk"), sep="\t", quote = FALSE, row.names = FALSE )

FC.Intermediate.HighLow <-  log2(FPKM.Data.NeoAntiBurd.Intermediate.mean + 1) - log2(HighLow.Mean + 1)
FC.Intermediate.HighLow <- cbind(expressionTMM.RPKM[,c("Chr","Start","End","GeneName","GeneID")], FC.Intermediate.HighLow)
FC.Intermediate.HighLow.rnk <- FC.Intermediate.HighLow %>% dplyr::arrange(-FPKM.Data.NeoAntiBurd.Intermediate.mean)
head(FC.Intermediate.HighLow.rnk); dim(FC.Intermediate.HighLow.rnk)
write.table(FC.Intermediate.HighLow.rnk[,c(4,6)], paste0("../RNASeq.RSEM/GSEA/rnk/FC.Intermediate.HighLow.pc.v6.rnk"), sep="\t", quote = FALSE, row.names = FALSE )

FC.Low.HighIntermediate <-  log2(FPKM.Data.NeoAntiBurd.Low.mean + 1) - log2(HighIntermediate.Mean + 1)
FC.Low.HighIntermediate <- cbind(expressionTMM.RPKM[,c("Chr","Start","End","GeneName","GeneID")], FC.Low.HighIntermediate)
FC.Low.HighIntermediate.rnk <- FC.Low.HighIntermediate %>% dplyr::arrange(-FPKM.Data.NeoAntiBurd.Low.mean)
head(FC.Low.HighIntermediate.rnk); dim(FC.Low.HighIntermediate.rnk)
write.table(FC.Low.HighIntermediate.rnk[,c(4,6)], paste0("../RNASeq.RSEM/GSEA/rnk/FC.Low.HighIntermediate.pc.v6.rnk"), sep="\t", quote = FALSE, row.names = FALSE )

############################################ New Addition Start #################
FPKM.intermediate.rnk <- FPKM.Data.NeoAntiBurd.Intermediate.mean %>% tibble::rownames_to_column(var="GeneName") %>% 
                                            dplyr::arrange(-FPKM.Data.NeoAntiBurd.Intermediate.mean) %>% 
                                            dplyr::mutate(FPKM.Data.NeoAntiBurd.Intermediate.mean = log2(FPKM.Data.NeoAntiBurd.Intermediate.mean+1))
head(FPKM.intermediate.rnk); dim(FPKM.intermediate.rnk)
write.table(FPKM.intermediate.rnk, paste0("../RNASeq.RSEM/GSEA/rnk/FPKM.intermediate.v7.rnk"), sep="\t", quote = FALSE, row.names = FALSE )

FPKM.High.rnk <- FPKM.Data.NeoAntiBurd.High.mean %>% tibble::rownames_to_column(var="GeneName") %>% 
                                      dplyr::arrange(-FPKM.Data.NeoAntiBurd.High.mean) %>%
                                      dplyr::mutate(FPKM.Data.NeoAntiBurd.High.mean = log2(FPKM.Data.NeoAntiBurd.High.mean+1))
head(FPKM.High.rnk); dim(FPKM.High.rnk)
write.table(FPKM.High.rnk, paste0("../RNASeq.RSEM/GSEA/rnk/FPKM.High.v7.rnk"), sep="\t", quote = FALSE, row.names = FALSE )
############################################ New Addition End #################

### Extract data from GSEA files ####
folder = "../RNASeq.RSEM/GSEA/results/NeoantigenVsImmueScore.PC.Cibersort.v6/"
allDirs                <- list.dirs(folder, recursive=FALSE);allDirs
preRankedGSEA.DF.NES  <- cbind(as.data.frame(lapply(allDirs, corUtilsFuncs$NESorPvalGSEAPrerank, colNumb=5))) 
preRankedGSEA.DF.Pval <- cbind(as.data.frame(lapply(allDirs, corUtilsFuncs$NESorPvalGSEAPrerank, colNumb=7)))
#preRankedGSEA.DF.Pval.Zeros <- as.data.frame(apply(preRankedGSEA.DF.Pval, 2, function(x){ x[x==0]<-0.001; return(x)}))
preRankedGSEA.DF.Pval.Zeros <- as.data.frame(apply(preRankedGSEA.DF.Pval, 2, function(x){ x = x+0.001; return(x)}))

High <- as.data.frame(cbind(NES=preRankedGSEA.DF.NES[,1], 
                            Pval=-log10(preRankedGSEA.DF.Pval.Zeros[,1]) ))
Intermediate <- as.data.frame(cbind(NES=preRankedGSEA.DF.NES[,2], 
                              Pval=-log10(preRankedGSEA.DF.Pval.Zeros[,2]) ))
Low    <- as.data.frame(cbind(NES=preRankedGSEA.DF.NES[,3], 
                              Pval=-log10(preRankedGSEA.DF.Pval.Zeros[,3]) ))

rownames(High) <- rownames(preRankedGSEA.DF.NES)
rownames(Intermediate) <- rownames(preRankedGSEA.DF.NES)
rownames(Low) <- rownames(preRankedGSEA.DF.NES)


### Plot the volcano plot ####
HighPlot <- ggplot(High, aes(x=NES, y=Pval, colour=NES>0 )) + 
  scale_colour_manual(name = 'NES > 0', values = setNames(c('red','darkgreen'),c(T, F))) +
  geom_point(size=3) +
  geom_text_repel(aes(x=NES, y=Pval, colour=NES>0, label=gsub("-CELLS_", ".",rownames(High))),  
                  size=3.5, force = 3) +
  geom_hline(yintercept=1.30103, size=1) +
  geom_vline(xintercept = 0, size=1) +
  xlim(-2,2.5)+
  ylim(0,3.5) +
  xlab("Normalised enrichment Score (NES)") +
  ylab("-log10(Pval)")+
  ggtitle("High Neoantigen Burden ")+
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        legend.position="none")

IntermediatePlot <- ggplot(Intermediate, aes(x=NES, y=Pval, colour=NES>0)) + 
  scale_colour_manual(name = 'NES > 0', values = setNames(c('red','darkgreen'),c(T, F))) +
  geom_point(size=3) +
  geom_text_repel(aes(x=NES, y=Pval, colour=NES>0, label=gsub("-CELLS_", ".",rownames(Intermediate))),  
                  size=3.5, force = 3) +
  geom_hline(yintercept=1.30103, size=1) +
  geom_vline(xintercept = 0, size=1) +
  xlim(-2,2)+
  ylim(0,3.5) +
  xlab("Normalised enrichment Score (NES)") +
  ylab("-log10(Pval)")+
  ggtitle("Intermediate Neoantigen Burden")+
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        legend.position="none")

LowPlot <- ggplot(Low, aes(x=NES, y=Pval, colour=NES>0)) + 
  scale_colour_manual(name = 'NES > 0', values = setNames(c('red','darkgreen'),c(T, F))) +
  geom_point(size=3) +
  geom_text_repel(aes(x=NES, y=Pval, colour=NES>0, label=gsub("-CELLS_", ".",rownames(Low))),  
                  size=3.5, force = 2) +
  geom_hline(yintercept=1.30103, size=1) +
  geom_vline(xintercept = 0, size=1) +
  xlim(-2,2)+
  ylim(0,3.5) +
  xlab("Normalised enrichment Score (NES)") +
  ylab("-log10(Pval)")+
  ggtitle("Low Neoantigen Burden") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        legend.position="none") 

pdf(paste0("../RNASeq.RSEM/Figures/Figure 4b.PC.v6",".pdf"), height=10, width = 20)
ggarrange(HighPlot, IntermediatePlot,LowPlot,  
          labels = c("A.", "B.", "C."),
          ncol = 3, nrow = 1,
          legend = NULL)
dev.off()

############################################### Neoantigen analysis for the OS #######################################
## Perform the same for OS samples only
neoantigenFromSamplesFinal %<>% filter(DIAGNOSIS.Substatus.Tumor.Normal.Tissue == "OS") %<>% data.table()

## Perform the same for OS samples only
## Following Samples removed because neoantigen burden for variant calling failed.
samplesToRemove <- c("OS.PALWWX", "OS.PANPUM", "OS.PAUUML", "OS.PAUYTT", "OS.PALZGU", "OS.PAMHYN", "OS.PANSEN")
expressionTMM.RPKM.Neoantigen.SR <- expressionTMM.RPKM.Neoantigen %>% 
                                dplyr::select(one_of("GeneName",neoantigenFromSamplesFinal$Sample.ID.Alias)) %>% 
                                dplyr::select(-one_of(samplesToRemove))
dim(expressionTMM.RPKM.Neoantigen.SR)


### Get mean expression of low neoantgen burden samples.
FPKM.Data.NeoAntiBurd.Low    <- expressionTMM.RPKM.Neoantigen.SR %>% dplyr::select( one_of( 
  c("GeneName",as.character(neoantigenFromSamplesFinal[TotalNeoantigenCount <= FirststQu, Sample.ID.Alias ]) ) ))
FPKM.Data.NeoAntiBurd.Low %<>% tibble::column_to_rownames(var="GeneName")
dim(FPKM.Data.NeoAntiBurd.Low); FPKM.Data.NeoAntiBurd.Low[1:5,1:5]

#FPKM.Data.NeoAntiBurd.Low.mean     <- as.data.frame(apply(log2(FPKM.Data.NeoAntiBurd.Low + 1), 1, mean)); 
FPKM.Data.NeoAntiBurd.Low.mean     <- as.data.frame(apply(FPKM.Data.NeoAntiBurd.Low, 1, mean)); 
colnames(FPKM.Data.NeoAntiBurd.Low.mean) <- "FPKM.Data.NeoAntiBurd.Low.mean"
head(FPKM.Data.NeoAntiBurd.Low.mean); dim(FPKM.Data.NeoAntiBurd.Low);

### Get mean expression of intermediate neoantgen burden samples.
FPKM.Data.NeoAntiBurd.Intermediate    <- expressionTMM.RPKM.Neoantigen.SR %>% dplyr::select( one_of( 
  c("GeneName",as.character(neoantigenFromSamplesFinal[TotalNeoantigenCount > FirststQu & TotalNeoantigenCount < thirdQu, Sample.ID.Alias ]) ) ))
FPKM.Data.NeoAntiBurd.Intermediate %<>% tibble::column_to_rownames(var="GeneName")
dim(FPKM.Data.NeoAntiBurd.Intermediate); FPKM.Data.NeoAntiBurd.Intermediate[1:5,1:5]

#FPKM.Data.NeoAntiBurd.Intermediate.mean     <- as.data.frame(apply(log2(FPKM.Data.NeoAntiBurd.Intermediate + 1), 1, mean)); 
FPKM.Data.NeoAntiBurd.Intermediate.mean     <- as.data.frame(apply(FPKM.Data.NeoAntiBurd.Intermediate, 1, mean))
colnames(FPKM.Data.NeoAntiBurd.Intermediate.mean) <- "FPKM.Data.NeoAntiBurd.Intermediate.mean"
head(FPKM.Data.NeoAntiBurd.Intermediate.mean); dim(FPKM.Data.NeoAntiBurd.Intermediate);

### Get mean expression of high neoantgen burden samples.
FPKM.Data.NeoAntiBurd.High    <- expressionTMM.RPKM.Neoantigen.SR %>% dplyr::select( one_of( 
  c("GeneName",as.character(neoantigenFromSamplesFinal[TotalNeoantigenCount >= thirdQu, Sample.ID.Alias ]) ) ))
FPKM.Data.NeoAntiBurd.High %<>% tibble::column_to_rownames(var="GeneName")
dim(FPKM.Data.NeoAntiBurd.High); FPKM.Data.NeoAntiBurd.High[1:5,1:5]

#FPKM.Data.NeoAntiBurd.High.mean     <- as.data.frame(apply(log2(FPKM.Data.NeoAntiBurd.High + 1), 1, mean)); 
FPKM.Data.NeoAntiBurd.High.mean     <- as.data.frame(apply(FPKM.Data.NeoAntiBurd.High, 1, mean)); 
colnames(FPKM.Data.NeoAntiBurd.High.mean) <- "FPKM.Data.NeoAntiBurd.High.mean"
head(FPKM.Data.NeoAntiBurd.High.mean); dim(FPKM.Data.NeoAntiBurd.High);

### Get expression of 2 categories for comparision ###
IntermediateLow <- cbind(FPKM.Data.NeoAntiBurd.Intermediate, FPKM.Data.NeoAntiBurd.Low)
#IntermediateLow.Mean <- as.data.frame(apply(log2(IntermediateLow+1), 1, mean)); colnames(IntermediateLow.Mean) <- "IntermediateLow"
IntermediateLow.Mean <- as.data.frame(apply(IntermediateLow, 1, mean)); colnames(IntermediateLow.Mean) <- "IntermediateLow"
dim(IntermediateLow.Mean); head(IntermediateLow.Mean)

HighLow <- cbind(FPKM.Data.NeoAntiBurd.High, FPKM.Data.NeoAntiBurd.Low)
#HighLow.Mean <- as.data.frame(apply(log2(HighLow+1), 1, mean)); colnames(HighLow.Mean) <- "HighLow"
HighLow.Mean <- as.data.frame(apply(HighLow, 1, mean)); colnames(HighLow.Mean) <- "HighLow"
dim(HighLow.Mean); head(HighLow.Mean)

HighIntermediate <- cbind(FPKM.Data.NeoAntiBurd.High, FPKM.Data.NeoAntiBurd.Intermediate)
#HighIntermediate.Mean <- as.data.frame(apply(log2(HighIntermediate+1), 1, mean)); colnames(HighIntermediate.Mean) <- "HighIntermediate"
HighIntermediate.Mean <- as.data.frame(apply(HighIntermediate, 1, mean)); colnames(HighIntermediate.Mean) <- "HighIntermediate"
dim(HighIntermediate.Mean); head(HighIntermediate.Mean)

### Find ratio between one vs rest of the groups ###
FC.High.IntermediateLow <-  log2(FPKM.Data.NeoAntiBurd.High.mean + 1) - log2(IntermediateLow.Mean + 1)
FC.High.IntermediateLow <- cbind(expressionTMM.RPKM[,c("Chr","Start","End","GeneName","GeneID")], FC.High.IntermediateLow)
FC.High.IntermediateLow.rnk <- FC.High.IntermediateLow %>% dplyr::arrange(-FPKM.Data.NeoAntiBurd.High.mean)
head(FC.High.IntermediateLow.rnk); dim(FC.High.IntermediateLow.rnk)
write.table(FC.High.IntermediateLow.rnk[,c(4,6)], paste0("../RNASeq.RSEM/GSEA/rnk/FC.High.IntermediateLow.pc.v6.OS.thirdQ.rnk"), sep="\t", quote = FALSE, row.names = FALSE )

FC.Intermediate.HighLow <-  log2(FPKM.Data.NeoAntiBurd.Intermediate.mean + 1) - log2(HighLow.Mean + 1)
FC.Intermediate.HighLow <- cbind(expressionTMM.RPKM[,c("Chr","Start","End","GeneName","GeneID")], FC.Intermediate.HighLow)
FC.Intermediate.HighLow.rnk <- FC.Intermediate.HighLow %>% dplyr::arrange(-FPKM.Data.NeoAntiBurd.Intermediate.mean)
head(FC.Intermediate.HighLow.rnk); dim(FC.Intermediate.HighLow.rnk)
write.table(FC.Intermediate.HighLow.rnk[,c(4,6)], paste0("../RNASeq.RSEM/GSEA/rnk/FC.Intermediate.HighLow.pc.v6.OS.thirdQ.rnk"), sep="\t", quote = FALSE, row.names = FALSE )

FC.Low.HighIntermediate <-  log2(FPKM.Data.NeoAntiBurd.Low.mean + 1) - log2(HighIntermediate.Mean + 1)
FC.Low.HighIntermediate <- cbind(expressionTMM.RPKM[,c("Chr","Start","End","GeneName","GeneID")], FC.Low.HighIntermediate)
FC.Low.HighIntermediate.rnk <- FC.Low.HighIntermediate %>% dplyr::arrange(-FPKM.Data.NeoAntiBurd.Low.mean)
head(FC.Low.HighIntermediate.rnk); dim(FC.Low.HighIntermediate.rnk)
write.table(FC.Low.HighIntermediate.rnk[,c(4,6)], paste0("../RNASeq.RSEM/GSEA/rnk/FC.Low.HighIntermediate.pc.v6.OS.thirdQ.rnk"), sep="\t", quote = FALSE, row.names = FALSE )

### Extract data from GSEA files ###
folder = "../RNASeq.RSEM/GSEA/results/NeoantigenVsImmueScore.PC.Cibersort.v6.OS.thirdQ/"
allDirs                <- list.dirs(folder, recursive=FALSE);allDirs
preRankedGSEA.DF.NES  <- cbind(as.data.frame(lapply(allDirs, corUtilsFuncs$NESorPvalGSEAPrerank, colNumb=5))) 
preRankedGSEA.DF.Pval <- cbind(as.data.frame(lapply(allDirs, corUtilsFuncs$NESorPvalGSEAPrerank, colNumb=7)))
#preRankedGSEA.DF.Pval.Zeros <- as.data.frame(apply(preRankedGSEA.DF.Pval, 2, function(x){ x[x==0]<-0.001; return(x)}))
preRankedGSEA.DF.Pval.Zeros <- as.data.frame(apply(preRankedGSEA.DF.Pval, 2, function(x){ x = x+0.001; return(x)}))

High <- as.data.frame(cbind(NES=preRankedGSEA.DF.NES[,1], 
                            Pval=-log10(preRankedGSEA.DF.Pval.Zeros[,1]) ))
Intermediate <- as.data.frame(cbind(NES=preRankedGSEA.DF.NES[,2], 
                                    Pval=-log10(preRankedGSEA.DF.Pval.Zeros[,2]) ))
Low    <- as.data.frame(cbind(NES=preRankedGSEA.DF.NES[,3], 
                              Pval=-log10(preRankedGSEA.DF.Pval.Zeros[,3]) ))

rownames(High) <- rownames(preRankedGSEA.DF.NES)
rownames(Intermediate) <- rownames(preRankedGSEA.DF.NES)
rownames(Low) <- rownames(preRankedGSEA.DF.NES)


### Plot the volcano plot ###
HighPlot <- ggplot(High, aes(x=NES, y=Pval, colour=NES>0 )) + 
  scale_colour_manual(name = 'NES > 0', values = setNames(c('red','darkgreen'),c(T, F))) +
  geom_point(size=3) +
  geom_text_repel(aes(x=NES, y=Pval, colour=NES>0, label=gsub("-CELLS_", ".",rownames(High))),  
                  size=3.5, force = 3) +
  geom_hline(yintercept=1.30103, size=1) +
  geom_vline(xintercept = 0, size=1) +
  xlim(-2,2.5)+
  ylim(0,3.5) +
  xlab("Normalised enrichment Score (NES)") +
  ylab("-log10(Pval)")+
  ggtitle("High Neoantigen Burden ")+
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        legend.position="none")

IntermediatePlot <- ggplot(Intermediate, aes(x=NES, y=Pval, colour=NES>0)) + 
  scale_colour_manual(name = 'NES > 0', values = setNames(c('red','darkgreen'),c(T, F))) +
  geom_point(size=3) +
  geom_text_repel(aes(x=NES, y=Pval, colour=NES>0, label=gsub("-CELLS_", ".",rownames(Intermediate))),  
                  size=3.5, force = 3) +
  geom_hline(yintercept=1.30103, size=1) +
  geom_vline(xintercept = 0, size=1) +
  xlim(-2,2)+
  ylim(0,3.5) +
  xlab("Normalised enrichment Score (NES)") +
  ylab("-log10(Pval)")+
  ggtitle("Intermediate Neoantigen Burden")+
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        legend.position="none")

LowPlot <- ggplot(Low, aes(x=NES, y=Pval, colour=NES>0)) + 
  scale_colour_manual(name = 'NES > 0', values = setNames(c('red','darkgreen'),c(T, F))) +
  geom_point(size=3) +
  geom_text_repel(aes(x=NES, y=Pval, colour=NES>0, label=gsub("-CELLS_", ".",rownames(Low))),  
                  size=3.5, force = 2) +
  geom_hline(yintercept=1.30103, size=1) +
  geom_vline(xintercept = 0, size=1) +
  xlim(-2,2)+
  ylim(0,3.5) +
  xlab("Normalised enrichment Score (NES)") +
  ylab("-log10(Pval)")+
  ggtitle("Low Neoantigen Burden") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        legend.position="none") 

pdf(paste0("../RNASeq.RSEM/Figures/Figure 4b.PC.v6.OS.median",".pdf"), height=10, width = 20)
ggarrange(HighPlot, IntermediatePlot,LowPlot,  
          labels = c("A.", "B.", "C."),
          ncol = 3, nrow = 1,
          legend = NULL)
dev.off()











### Violin plot for exhaustion markers ####

RPKM.Data.Exhaustion   <- expressionTMM.RPKM %>% dplyr::filter(GeneName %in% as.character(rnaseqProject$emDF$GeneName)) %>% dplyr::arrange(GeneName)
Exhaustion.Transpose           <- as.data.frame(t(RPKM.Data.Exhaustion[,-c(1:7)]))
colnames(Exhaustion.Transpose) <- RPKM.Data.Exhaustion$GeneName
Exhaustion.Transpose <- Exhaustion.Transpose %>% tibble::rownames_to_column(var="Sample.ID.Alias")
Exhaustion.Transpose.diag <- dplyr::full_join(Exhaustion.Transpose, rnaseqProject$metaDataDF[,c("Sample.ID.Alias", "DIAGNOSIS.Substatus.Tumor.Normal.Tissue")], 
                                              by="Sample.ID.Alias") %>% dplyr::rename(Diagnosis=DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
Rm.Normal.Exhaustion.Transpose <- Exhaustion.Transpose.diag %>% filter(!grepl("^NS.*", Exhaustion.Transpose.diag$Diagnosis))
finalExhaustionMatrix <-  melt(Rm.Normal.Exhaustion.Transpose[,-1], id.var = "Diagnosis")

finalExhaustionMatrix.tidy <- finalExhaustionMatrix %>% dplyr::group_by(Diagnosis, variable) %>% 
  dplyr::mutate(Med=median(value)) %>% arrange(Diagnosis, variable, value) %>% 
  arrange(desc(Med)) %>% 
  ungroup() %>% 
  mutate( Diagnosis.Marker = factor(paste(Diagnosis,variable,sep="."), levels= unique(paste(Diagnosis,variable,sep=".")),
                                    order = TRUE) ) %>% 
  arrange(Diagnosis.Marker)

pdf("ExhautionMarkers.Variable.RPKM.v5.pdf",height=25,width=20)
ggplot(finalExhaustionMatrix.tidy, aes(x=Diagnosis.Marker, y=value)) + 
  ##ggplot(data, aes(x=Group, y=log2(ENSG00000182752))) + 
  #geom_boxplot(varwidth = TRUE,notch = FALSE) + 
  geom_violin(scale = "width",trim = FALSE, draw_quantiles = c(0.5)) + 
  #geom_jitter(width=0.1) +
  theme_bw() + 
  ylab( paste("log2(TPM)") ) +
  #ylim(0,10) +
  xlab( "Diagnosis" ) +
  #geom_hline(yintercept=0, size=0.1) + 
  theme( title = element_text(size=13, face="bold")
         ,axis.title.x = element_text(size=13, face="bold")
         ,axis.title.y = element_text(size=13, face="bold")
         ,axis.text.x = element_text(size=10, face="bold", angle=90, vjust=1)
         ,axis.text.y = element_text(size=10, face="bold")
         ,axis.ticks.x =element_blank()
         ,strip.text.y= element_blank()
         ,strip.text.x=element_text(size=13,face="bold")
         ,strip.background=element_blank()
         ,panel.grid.major.x=element_blank()
         ,panel.grid.minor.x=element_blank()
         ,panel.border = element_rect(colour = "black", fill=NA, size=0.0000000002, linetype = 2)
         ,panel.spacing = unit(0, "cm")
         ,strip.switch.pad.grid = unit(0, "cm")
  ) + facet_wrap( ~ variable, scales="free", nrow = 7) +
  scale_x_discrete(labels=setNames(as.character(finalExhaustionMatrix.tidy$Diagnosis), finalExhaustionMatrix.tidy$Diagnosis.Marker))
dev.off()

############ TP53 Analysis ###################################


TP53.mutant.samples <- read.table("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/TP53.Unique.txt",header = T)
TP53.mutant.samples <- left_join(TP53.mutant.samples, rnaseqProject$validMetaDataDF[,c("Sample.Biowulf.ID.GeneExp","Sample.ID","DIAGNOSIS.Alias")], 
                                 by="Sample.ID") %>% arrange(Sample.ID)
HLA_geneexp      <- rbind(expressionTMM.RPKM.GSEA.Input[c("HLA-A", "HLA-B", "HLA-C"),])
HLA_geneexp.t <- t(HLA_geneexp) %>% data.frame() %>% tibble::rownames_to_column(var = "Sample.Biowulf.ID.GeneExp") 
HLA_geneexp.t <- HLA_geneexp.t %>% dplyr::mutate(Group = ifelse(Sample.Biowulf.ID.GeneExp %in% 
                                                                  TP53.mutant.samples$Sample.Biowulf.ID.GeneExp, "TP53-mutant", "TP53-wildtype"))
HLA_geneexp.t <- left_join(HLA_geneexp.t, rnaseqProject$validMetaDataDF[,c("Sample.Biowulf.ID.GeneExp", "DIAGNOSIS.Alias", "LIBRARY_TYPE")], 
                           by="Sample.Biowulf.ID.GeneExp")
HLA_geneexp.t <- HLA_geneexp.t %>% dplyr::mutate(Group = ifelse(DIAGNOSIS.Alias == "NS", "Normal", Group))

HLA_geneexp.t <- HLA_geneexp.t %>% dplyr::filter(!grepl("CellLine", LIBRARY_TYPE))

my_comparisons=list(c("TP53-mutant","TP53-wildtype"), c("TP53-mutant", "Normal"), c("TP53-wildtype","Normal"))


pdf("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/HLA-Tumor-samples.pdf", height = 20, width = 15)
HLA.A <- ggplot(HLA_geneexp.t, aes(x=Group, y=HLA.A , fill=Group)) + geom_boxplot() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
                        geom_jitter(shape=16, position=position_jitter(0.2)) +  theme(legend.position = "none") +
                        stat_compare_means( paired = FALSE,comparisons = my_comparisons) +
                        stat_compare_means(method = "anova", label.y = 16) + 
                        ylab("HLA.A RPKM")
HLA.B <- ggplot(HLA_geneexp.t, aes(x=Group, y=HLA.B, fill=Group)) + geom_boxplot()  + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
                        geom_jitter(shape=16, position=position_jitter(0.2)) +  theme(legend.position = "none") +
                        stat_compare_means( paired = FALSE,comparisons = my_comparisons) +
                        stat_compare_means(method = "anova", label.y = 12.5) + 
                        ylab("HLA.B RPKM")
                       
HLA.C <- ggplot(HLA_geneexp.t, aes(x=Group, y=HLA.C , fill=Group)) + geom_boxplot()  + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
                        geom_jitter(shape=16, position=position_jitter(0.2)) +  theme(legend.position = "none") +
                        stat_compare_means( paired = FALSE,comparisons = my_comparisons) +
                        stat_compare_means(method = "anova", label.y = 12.5) + 
                        ylab("HLA.C RPKM")
grid.arrange(HLA.A,HLA.B,HLA.C,nrow = 3)
dev.off()

################################################################# ERV Analysis #################################################

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

## coorelation plots
correlationPlots <- function(varName="", constName="", df=NA, customColorDF=NA, xlab="Log Total Clones"){
  
  print(paste(varName))
  df <- df %>% dplyr::filter_(.dots = paste0(varName, ">=", 5))
  customColorsVector <- setNames( as.character(customColorDF$Color), as.character(customColorDF$Diagnosis))
  corrTest <- cor.test(df[,constName], df[,varName], method = "spearman")
  if  ( corrTest$p.value < 2.2e-16 ) { corrTest$p.value = 2.2e-16 }
  plot <- ggplot(df, aes_string(x=constName, y=varName)) + 
    geom_smooth(method=lm,  fill="grey") +
    geom_point(aes(colour = factor(Diagnosis)), show.legend = T, size=2, shape=16) + 
    scale_colour_manual(values=customColorsVector) +
    theme_bw() +
    theme(axis.text=element_text(size=13)
          ,axis.title=element_text(size=13,face="bold")) +
    xlab(xlab) +
    ylab(paste("Erv expression", sep=" "))+
    # theme(legend.position = "none") +
    ggtitle(paste("Corr.Coeff = ", signif(corrTest$estimate[[1]],5), "\np-value = ", signif(corrTest$p.value,5), "", varName,sep=""))
  
  return(list(plot))
} 

## Read ERV expression file

## Need to improve the below code
## Log the expression and saved  ## NO need to use everytime.
# erv.data.max <- read.csv("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/data.geneNames.merged.max.txt", sep="\t", header = T)
# colnames(erv.data.max) <- gsub(".rsem.ENS.FPKM","",colnames(erv.data.max))
# erv.data.max <- erv.data.max %>% tibble::column_to_rownames(var="geneName")
# erv.data.max <- log(erv.data.max+1,2)
# # write.table(erv.data.max, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/data.geneNames.merged.max.log2.txt", sep="\t", quote = F)
# ## Get cytolytic score and save
# cytolyticScore          <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM.GSEA.Input)
# write.table(cytolyticScore, "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/cytolyticScore.txt", sep="\t", quote = F)

## Read the final expression file
erv.data.max.cytolytic <- read.csv("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/ERV_final_files/data.geneNames.merged.max.cytolytic.v2.log2.txt",
                                     sep="\t", header = T)
# erv.data.max.cytolytic.complete <- erv.data.max.cytolytic[complete.cases(erv.data.max.cytolytic),] %>% data.frame()
erv.data.max.cytolytic <- tibble::column_to_rownames(erv.data.max.cytolytic, var = "GeneID")

## Add additional annotations (sample Id alias) ####
AliasNames_df                 <- dplyr::left_join( data.frame("Sample.Data.ID"=colnames(erv.data.max.cytolytic)), 
                                          rnaseqProject$validMetaDataDF[,c("Sample.Biowulf.ID.GeneExp",
                                          "Sample.ID.Alias", "Sample.Data.ID", "DIAGNOSIS.Alias" ,"Color",
                                          rnaseqProject$factorName)] )
## Remove Normal samples
Tumor_samples_annot <- AliasNames_df %>% dplyr::filter(!grepl('Teratoma',DIAGNOSIS.Substatus.Tumor.Normal.Tissue))
Tumor_samples_annot <- Tumor_samples_annot[complete.cases(Tumor_samples_annot),]
erv.data.max.cytolytic.Tumor <- erv.data.max.cytolytic %>% dplyr::select(one_of(Tumor_samples_annot$Sample.Data.ID)); dim(erv.data.max.cytolytic.Tumor)

## For all tumor samples
erv.data.max.cytolytic.Tumor.corr <- rcorr(t(erv.data.max.cytolytic.Tumor), type = "spearman");
erv.data.max.cytolytic.Tumor.corr.cytolytic <- flattenCorrMatrix(erv.data.max.cytolytic.Tumor.corr$r, erv.data.max.cytolytic.Tumor.corr$P) %>% filter(grepl('CytolyticScore', column))
write.table(erv.data.max.cytolytic.Tumor.corr.cytolytic, 
             "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/erv.data.max.cytolytic.Tumor.all.corr.txt", quote = FALSE, sep="\t")

Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue <- gsub('NS.*','NS',Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue)
## Perform rcorr for each tumor groups 
correlationDF <- rcorr_groups(Tumor_samples_annot, erv.data.max.cytolytic.Tumor, "DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Sample.Data.ID")
correlationDF_cytolytic <- correlationDF %>% filter(grepl('CytolyticScore', column))
write.table(correlationDF_cytolytic, 
            "T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/erv.data.max.cytolytic.Tumor.each.corr.txt", quote = FALSE, sep="\t")


### Plot correlation plot for the following ERVs ###
sig.Ervs <- c("LTR12E_LTR/ERV1","LTR16E2_LTR/ERVL","LTR23_LTR/ERV1","LTR70_LTR/ERV1",
              "LTR86B2_LTR/ERVL","MER11B_LTR/ERVK","MER52A_LTR/ERV1","MER61C_LTR/ERV1",
              "MER83C_LTR/ERV1","MLT2B5_LTR/ERVL")

subset_data_annot <- cbind(t(erv.data.max.cytolytic.Tumor[c(sig.Ervs, "CytolyticScore"),]),
                           Diagnosis=as.character(Tumor_samples_annot$DIAGNOSIS.Substatus.Tumor.Normal.Tissue)) %>% data.table()

#colnames(subset_data_annot) <- gsub("/",".", colnames(subset_data_annot))
## Generate custom colors
customColorDF    <- rbind(rnaseqProject$customColorsDF, data.frame(Diagnosis=c("NS"), Color=c("#000000")))
sig.Ervs.valid <- sub("/",".", sig.Ervs) ; cols = sig.Ervs.valid
dataDF <- data.frame(subset_data_annot)
dataDF[c(sig.Ervs.valid, "CytolyticScore")] <- sapply(dataDF[c(sig.Ervs.valid, "CytolyticScore")],as.numeric)
dataDF$Diagnosis <- gsub('^NS.*','NS',dataDF$Diagnosis)
dataDF %<>% dplyr::filter(!grepl("NS",Diagnosis))
#dataDF.RMS.FP <- dataDF %>% dplyr::filter(Diagnosis %in% c("RMS.FN"))

## Plot all of them
plotLists <- lapply(sig.Ervs.valid, correlationPlots, constName="CytolyticScore", xlab="CTL Score", df=dataDF, 
                          customColorDF=customColorDF)
ERVCytolyticScorePlots <- lapply(plotLists, function(l) l[[1]] )
ggsave("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/Figures/erv.correlation_plots.pdf", marrangeGrob(ERVCytolyticScorePlots, ncol=1, nrow=1),
       width = 15, height = 10)

## plot individual ERV specific to diagnosis
sig.Ervs.specific <- "MLT2B5_LTR.ERVL"
dataDF.diagnosis <- dataDF %>% dplyr::filter(grepl("NB.Unknown",Diagnosis)); dim(dataDF.diagnosis)
plotLists <- lapply(sig.Ervs.specific, correlationPlots, constName="CytolyticScore", xlab="CTL Score", df=dataDF.diagnosis, 
                    customColorDF=customColorDF)
ggsave(paste0("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/Figures/",sig.Ervs.specific,".pdf"),width = 10, height = 10)



### Analyse Nanostring datasets
## Function declaration
plotViolins <- function(x, df = NA, group = NA , ylab = NA, xlab = NA, title=NA, fill = NA){
  genemarker <- x
  # genemarker <- "4.1BB"
  # group <- "Risk"
  # ylab <- "Standardised ROI count"
  # xlab <- "Risk"
  customColorsVector <- c('0' = "#C7C2B8", '1'='#E69F00' )
  if(is.na(title)){
    title = genemarker
  } else {
    title = paste(title,genemarker,sep="-")
  }
  if(is.na(xlab)){
    xlab = group
  } else {
    xlab = xlab
  }
  if(is.na(fill)){
    fill = group
  } else {
    fill = fill
  }

  p <- ggplot(df, aes_string(x = group, y = genemarker, fill= group)) + 
    #geom_boxplot(aes_string(fill = group)) +
    geom_violin(aes_string(fill = group)) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4, position=position_dodge(1),
                 fill=c("black")) +
    scale_fill_manual(values=customColorsVector) +
    theme_bw() +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title) +
    #ylim(0,34)+
    theme(axis.title=element_text(face="bold",size="17"), text = element_text(size=10, face="bold"),
          plot.title = element_text(hjust=0.5),
          legend.text= element_text(face="bold",size="12"),
          strip.text.x=element_text(size=10,face="bold", angle=90),
          legend.position="none") +
    stat_compare_means( method = "t.test", paired = FALSE, comparisons = list(c("0","1")) )
  return(p)
}

## Read nanostring dataset
nanostringData <- read.csv("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/NanoString/ROI_Map_ClinicalInfo_with Zscored ImmuneMarker_shortlist_send_For Sivasish.txt",
                            sep="\t", header = TRUE)
## Factor variables to use as Groups/categories
nanostringData$Risk <- as.factor(nanostringData$Risk)
nanostringData$NMYC <- as.factor(nanostringData$NMYC)

## Get the average of duplicate samples
genemarkers <- colnames(nanostringData)[-c(1:8)]
nanostringData.group <- nanostringData %>% group_by(Samples..TMA.14.001.) %>% 
                        summarize_at(vars(one_of(genemarkers)), mean)
nanostringData.group.final <- dplyr::left_join(nanostringData[,c("Samples..TMA.14.001.", "NMYC", "Risk", "Post.Treatment")], nanostringData.group, 
                                               by="Samples..TMA.14.001.")

## Plot Violin plots for all markers ## Individual plot for each marker
plotLists.all <- lapply(genemarkers, plotViolins, df = nanostringData.group.final, group="NMYC", title="All", ylab="Standardised marker counts")

nanostringData.post.treat <- nanostringData.group.final %>% dplyr::filter(grepl('post', Post.Treatment));dim(nanostringData.post.treat)
plotLists.post <- lapply(genemarkers, plotViolins, df = nanostringData.post.treat, group="NMYC", title="Post", ylab="Standardised marker counts")

nanostringData.no.treat <- nanostringData.group.final %>% dplyr::filter(!grepl('post', Post.Treatment));dim(nanostringData.no.treat)
plotLists.no.treat <- lapply(genemarkers, plotViolins, df = nanostringData.no.treat, group="NMYC", title="Pre", ylab="Standardised marker counts")

all_plotlists <- c(rbind(plotLists.no.treat, plotLists.post, plotLists.all))
all_plotlists <- plotLists.all

ggsave(paste0("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/Figures/Nanostring.NMYC",".pdf"), marrangeGrob(all_plotlists, ncol=3, nrow = 1),
       width = 15, height = 10)

## Plot violin plots for all markers under same scale
nanostringData.group.final.tidy   <- tidyr::gather_(nanostringData.group.final, key="GeneMarker", value="Score", genemarkers )
immuneMarkers <- c("beta.2.microglobulin","CD45", "HLA.DR", "CD3", "CD4", "CD8", "CD45RO", "GZMB", "X4.1BB", "CD40", "CD40L",
                   "CD86", "CD11c", "CD14", "CD163", "CD68", "CD20", "CD25", "CD27", "CD34", "CD44", "CD56", "CD80", "FoxP3",
                   "GITR", "CTLA4", "ICOS", "IDO.1", "LAG3", "PD.1", "PD.L1", "STING.TMEM173", "TGF.beta.1", "TIM.3", "VISTA")

nanostringData.group.tidy.immune <- nanostringData.group.final.tidy %>% filter(GeneMarker %in% immuneMarkers); dim(nanostringData.group.tidy.immune)
plot <- plotViolins(x = "Score", df = nanostringData.group.tidy.immune, group = "GeneMarker", fill="NMYC", ylab="Standardised marker counts")

### Box plot
plot <- ggplot(nanostringData.group.tidy.immune, aes_string(x="GeneMarker", y = "Score")) + 
            geom_boxplot(aes_string(fill = "NMYC"),
                         position = position_dodge(0.9)) +
          scale_fill_manual(values = c("#C7C2B8", "#E69F00")) +
          theme(axis.title=element_text(face="bold",size="17"), 
                text = element_text(size=10, face="bold"),
                axis.text.x = element_text(angle = 90, hjust = 0.5, size =8, face = "bold"),
                plot.title = element_text(hjust=0.5),
                legend.text= element_text(face="bold",size="12"),
                strip.text.x=element_text(size=10,face="bold", angle=90),
                legend.position="none") +
          stat_compare_means(method = "t.test", paired = FALSE, comparisons = list(c("0","1"))) +
          ggtitle("Nanostring immunemarker")
ggsave(paste0("T:/Sivasish_Sindiri/R Scribble/RNASeq.RSEM/Figures/Nanostring.NMYC.singlePlot",".pdf"), plot,
                 width = 45, height = 10)
