setwd("T:/Sivasish_Sindiri/R_workspace/MiXCR/")

## Khanlab meta data
metaData <- read.csv("MetadataMapper.v3.txt", sep="\t")

############################################ Section 0 Declare functions ##########################################################

### Required Custom Functions ####

filterSpecificCloneTypes  <- function(cloneData, cloneType){
  cloneDataFilt           <- cloneData %>% dplyr::filter(grepl(cloneType,v) | grepl("NF",v))
  return(cloneDataFilt)
}

correlationPlots <- function(varName="", constName="", df=NA){
  
  print(paste(varName))
  # plot <- ggscatter(df, x = constName, y =varName, 
  #                   add = "reg.line", conf.int = TRUE, 
  #                   cor.coef = TRUE, cor.method = "spearman",
  #                   xlab = constName, ylab = varName)
  
  corrTest <- cor.test(df[,constName], df[,varName], method = "spearman")
  if  ( corrTest$p.value < 2.2e-16 ) { corrTest$p.value = 2.2e-16 }
  plot <- ggplot(df, aes_string(x=constName, y=varName)) + 
    geom_smooth(method=lm,  fill="grey") +
    geom_point(aes(colour = factor(DIAGNOSIS.Substatus.Tumor.Normal.Tissue)), show.legend = T, size=4) + 
    scale_colour_manual(values=setNames( StatsFinal$Color, StatsFinal$DIAGNOSIS.Substatus.Tumor.Normal.Tissue)  ) +
    theme_bw() +
    theme(axis.text=element_text(size=13)
          ,axis.title=element_text(size=13,face="bold")) +
    xlab("Log Total Clones") +
    ylab(paste("Standardised Enrichment Score", sep=" "))+
    ggtitle(paste("Corr.Coeff = ", signif(corrTest$estimate[[1]],5), "\np-value = ", signif(corrTest$p.value,5), "                                            ", varName,sep=""))
  
  return(list(plot))
}   

makeEntropyInput <- function(filename, inputDir="", outputDir="", cloneType="") {
  outfileName <- paste0(outputDir, gsub("convert.|.clones.txt","",filename ), ".Entropy.txt")
  exomeData <- read.csv( paste(inputDir, filename, sep=""), sep="\t", header = TRUE )
  
  if(cloneType != "") {
    exomeData               <- filterSpecificCloneTypes(cloneData = exomeData, cloneType = cloneType)
  }
  
  if(nrow(exomeData)>0){
      exomeDataEntropy <- data.frame(VJcombo=paste(exomeData$v,exomeData$j,sep="."), Counts=exomeData$count, Vcassette=exomeData$v, 
                                 Jcassette=exomeData$j, aaCDR3_filtered = exomeData$cdr3aa, ntCDR3= exomeData$cdr3nt)
  } else {
    exomeDataEntropy <- emptyDFEntropy
  }
  write.table(exomeDataEntropy, outfileName, sep = "\t", row.names = FALSE, quote = FALSE)
}

immunoseqv2 <- function(x) {
  filename = x
  outfileName <- paste0("./immunoseqv2EntropyNitin/", gsub("convert.|.clones.txt","",filename ), ".Entropy.txt")
  exomeData <- read.csv( paste("./immunoseqv2/", filename, sep=""), sep="\t", header = TRUE, stringsAsFactors = FALSE )
  exomeDataFilt <- exomeData[,1:47] %>% dplyr::filter(sequenceStatus %in% c("In"))
  colnames(exomeDataFilt)[3] <- "count..templates.reads."
  
  print(paste("./immunoseqv2/", filename, sep=""))
  print(dim(exomeDataFilt))
  
  if(nrow(exomeDataFilt)>0){
    
    emptyVGeneName <- which(exomeDataFilt$vGeneName == "")
    emptyJGeneName <- which(exomeDataFilt$jGeneName == "")
    exomeDataFilt[emptyVGeneName, c("vGeneName")] <- sapply(exomeDataFilt[emptyVGeneName,c("vGeneNameTies")], function(x){  unlist(strsplit(as.character(x), ","))[1] })
    exomeDataFilt[emptyJGeneName, c("jGeneName")] <- sapply(exomeDataFilt[emptyJGeneName,c("jGeneNameTies")], function(x){  unlist(strsplit(as.character(x), ","))[1] })
    
    
    
    exomeDataEntropy <- data.frame(VJcombo=paste(exomeDataFilt$vGeneName,exomeDataFilt$jGeneName,sep="."), 
                                   Counts=exomeDataFilt$count..templates.reads., 
                                   Vcassette=exomeDataFilt$vGeneName, 
                                   Jcassette=exomeDataFilt$jGeneName, 
                                   aaCDR3_filtered = exomeDataFilt$aminoAcid, 
                                   ntCDR3= exomeDataFilt$nucleotide)
  } else {
    exomeDataEntropy <- emptyDFEntropy
  }
  write.table(exomeDataEntropy, outfileName, sep = "\t", row.names = FALSE, quote = FALSE)
  
}

############################################ Section 1 Read data  #################################################################

### Placeholder DF ####
emptyDF <- data.frame(count=c(0), freq=c(0), cdr3nt=c("NA"),cdr3aa=c("NF"),v=c("NF"),d=c("NF"),j=c("NF"),VEnd=c(0),DStart=c(0),
                      DEnd=c(0),JStart=c(0),SampleName=c(0))
emptyDFEntropy <- data.frame(VJcombo=c(), Counts =c(), Vcassette=c(), Jcassette=c(), aaCDR3_filtered = c(), ntCDR3= c())

emptyDFEntropyResults <- data.frame(FileName=c(), Hcdr3 =c(), Htot=c(), CLcdr3=c(), CLHvj= c(), CLtot= c(),
                                    Hcdr3_max=c(), Hvj_max =c(), Htot_max=c(), CLcdr3_max=c(), Num_CDR3= c(), Num_VJ= c(),
                                    Num_totCDR3 =c())

### List files and read data into a single data matrix ####
fileList <- list.files("./CloneFiles.v2/")
#fileList <- c("convert.Sample_RMS248_C14C7ACXX.clones.txt")
AllClonesData             <- rbindlist( lapply(fileList, function(x){
  print(x)
  exomeData <- read.csv( paste("./CloneFiles.v2/", x, sep=""), sep="\t", header = TRUE )
  if(nrow(exomeData)>0){
    exomeData$SampleName <- x
  } else {
    emptyDF$SampleName <- c(x)
    exomeData <- emptyDF
  }
  return(exomeData)
}) )

### Filter Clones by clone types ####
cloneObjIG                <- filterSpecificCloneTypes(cloneData = AllClonesData, cloneType = "IGH")
cloneObjIG.Expansion.GE3  <- cloneObjIG %>%  dplyr::filter(grepl("IGH",v) & count >= 3)
cloneObjTCR               <- filterSpecificCloneTypes(cloneData = AllClonesData, cloneType = "TRB")
cloneObjTCR.Expansion.GE3  <- cloneObjTCR %>% dplyr::filter(grepl("TRB",v) & count >= 3)

### Read normalised counts for each sample
readCounts <- readRDS("T:/Sivasish_Sindiri/R_workspace/MiXCR/RNASeq.readcounts.rds")
readCountsSum <- apply(readCounts, 2, sum)
readCountsSum <- as.data.frame(readCountsSum)
readCountsSum.df <- readCountsSum %>% tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")

### write the inputs for entropy data ####
# fileList <- list.files("./CloneFilesNitin/")
# AllClonesEntropyData             <- sapply(fileList, makeEntropyInput, inputDir="./CloneFilesNitin/", outputDir="./CloneFilesEntropyNitin/" )

fileList <- list.files("./CloneFiles.v2/")
cloneType = "IGH"
AllClonesEntropyData             <- sapply(fileList, makeEntropyInput,  cloneType=cloneType, inputDir="./CloneFiles.v2/", 
                                           outputDir=paste0("./CloneFilesEntropy.", cloneType, ".v2/") )

### write the inputs for ImmunoseqV2 entropy data ####
fileList <- list.files("./immunoseqv2/")
ImmunoseqV2EntropyData             <- sapply(fileList, immunoseqv2 )

# ### List files and read data into a single data matrix for Entropy results####
# fileList <- list.files("./finalResults_H_CL_JS_Nitin/")
# AllEntropyData             <- rbindlist( lapply(fileList, function(x){
#   print(x)
#   exomeData <- read.csv( paste("./finalResults_H_CL_JS_Nitin/", x, sep=""), sep="\t", header = TRUE )
#   if(nrow(exomeData) == 0){
#     exomeData <- emptyDFEntropyResults
#   }
#   return(exomeData[1,])
# }) )
# write.table(AllEntropyData, "AllEntropyData_H_CL_JS.Nitin.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# 
# ### List files and read data into a single data matrix for Entropy results####
# fileList <- list.files("./finalResults_H_VL_JS_Nitin.v3/")
# AllEntropyData             <- rbindlist( lapply(fileList, function(x){
#   print(x)
#   exomeData <- read.csv( paste("./finalResults_H_VL_JS_Nitin.v3/", x, sep=""), sep="\t", header = TRUE )
#   print(dim(exomeData))
#   if(nrow(exomeData) == 0){
#     exomeData <- emptyDFEntropyResults
#   }
#   return(exomeData[1,])
# }) )
# write.table(AllEntropyData, "AllEntropyData_H_CL_JS.Nitin.ImmunoseqV4.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### List files and read data into a single data matrix for Entropy results####
folderName = "./finalResults_H_VL_JS_landscape.TRB.v3/"
fileList <- list.files(folderName)
AllEntropyData             <- rbindlist( lapply(fileList, function(x){
  print(x)
  exomeData <- read.csv( paste(folderName, x, sep=""), sep="\t", header = TRUE )
  print(dim(exomeData))
  if(nrow(exomeData) == 0){
    exomeData <- emptyDFEntropyResults
  }
  return(exomeData[1,])
}) )
write.table(AllEntropyData, "AllEntropyData_H_CL_JS.landscape.TRB.v3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

AllEntropyData$FileName <- gsub(".Entropy.tx","",AllEntropyData$FileName)
AllEntropyData %<>% dplyr::rename(Sample.ID=FileName)
AllEntropyData.annot <- dplyr::full_join(AllEntropyData, metaData[,c("Sample.ID", "SAMPLE_ID.Alias","LIBRARY_TYPE","DIAGNOSIS.Alias")], by="Sample.ID")
  



############################################ Section 1b Select DF ################################################################
### To Do  ####
### Add switch case to select Clone df based on user input ###

### For now Select DF manually ####
#cloneType = "IGHClones"  ; countObj <- cloneObjIG %>% as.data.frame()
cloneType = "TRBClones"  ; countObj <- cloneObjTCR %>% as.data.frame()

### Attach metadata and generate countObj ####
countObj <- countObj %>% dplyr::rename(Sample.Data.ID=SampleName); 
countObj$Sample.Data.ID <- gsub("convert.|.clones.txt","", countObj$Sample.Data.ID)
#countObj$SAMPLE_ID <- gsub("-","_", countObj$SAMPLE_ID)
countObj.Annot <- dplyr::left_join(countObj, metaData, by="Sample.Data.ID") %>% 
  dplyr::select_(.dots=c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart", "Sample.Biowulf.ID.GeneExp", "Sample.ID.Alias",
                         "LIBRARY_TYPE","DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Color.Substatus" )) ; head(countObj.Annot)
### Plot the clone expansion
countObj.Annot.NoCL <- countObj.Annot %>% filter(!grepl('CellLine',LIBRARY_TYPE)) %>% filter(!grepl('^NS', DIAGNOSIS.Substatus.Tumor.Normal.Tissue) )

countObj.Annot.complete <- countObj.Annot.NoCL[complete.cases(countObj.Annot.NoCL),]
## sanity check
dim(countObj)
dim(countObj.Annot)
dim(countObj.Annot.NoCL)
dim(countObj.Annot.complete)


######## Make Step Plots to show expansion ############
countObj.Annot.NoCL.totalReads <- dplyr::left_join(countObj.Annot.complete, readCountsSum.df, by="Sample.Biowulf.ID.GeneExp")
countObj.Annot.NoCL.totalReads.complete <- countObj.Annot.NoCL.totalReads[complete.cases(countObj.Annot.NoCL.totalReads),]
dim(countObj.Annot.NoCL.totalReads.complete)

countObj.Annot.NoCL.totalReads$ReadsPerMillion <- ( countObj.Annot.NoCL.totalReads$count/countObj.Annot.NoCL.totalReads$readCountsSum)*1000000

countObj.Annot.NoCL.totalReads <- countObj.Annot.NoCL.totalReads %>% dplyr::select(Sample.Biowulf.ID.GeneExp, 
                                                                                   DIAGNOSIS.Substatus.Tumor.Normal.Tissue, 
                                                                                   count, readCountsSum, ReadsPerMillion,
                                                                                   Color.Substatus)

#toPlotDF <- countObj.Annot.NoCL.totalReads %>% dplyr::mutate(ReadsPerMillion = if_else(ReadsPerMillion >= 2, 2, ReadsPerMillion))

toPlotDF <- countObj.Annot.NoCL.totalReads %>% dplyr::filter(count > 0) %>% 
  arrange(Sample.Biowulf.ID.GeneExp, count) %>%
  group_by(Sample.Biowulf.ID.GeneExp) %>% 
  mutate(rank = dense_rank( -count )) %>% 
  distinct()
  # mutate(good_ranks = order(order(order_values, decreasing=TRUE)))
View(toPlotDF)

ggplot(toPlotDF[,c(1,3,5,6,7)]) +
  geom_step(aes(y = rank, x = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp,colour=Color.Substatus), 
            size=0.5 ) +
  facet_wrap(~toPlotDF$DIAGNOSIS.Substatus.Tumor.Normal.Tissue) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.text=element_text(size=12),
        axis.text = element_text(size=11),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(trans = "reverse", breaks = seq(1,max(toPlotDF$rank),by=4) ) +
  coord_trans(x = "log2" ) +
  scale_x_continuous(minor_breaks = c(),
                     breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4),
                     labels = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4) ) 

  
toPlotDF.NB.MYCN.NA <- countObj.Annot.NoCL.totalReads %>% filter(grepl('NB.MYCN.NA', DIAGNOSIS.Substatus.Tumor.Normal.Tissue))
NB.A <- ggplot(toPlotDF.NB.MYCN.NA[,c(1,3,5,6)]) +
  geom_step(aes(y = count, x = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp,colour=Color.Substatus), size=0.6 ) +
  facet_zoom(xy = ReadsPerMillion < 50 & count < 2500 ) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle(label = "NB.MYCN.NA") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size=15),
        axis.title = element_text(size = 18))

toPlotDF.NB.MYCN.A <- countObj.Annot.NoCL.totalReads %>% filter(grepl('NB.MYCN.A', DIAGNOSIS.Substatus.Tumor.Normal.Tissue))
NB.NA <- ggplot(toPlotDF.NB.MYCN.A[,c(1,3,5,6)]) +
  geom_step(aes(y = count, x = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp,colour=Color.Substatus), size=0.6 ) +
  facet_zoom(xy = ReadsPerMillion < 50 & count < 1000  ) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle(label = "NB.MYCN.A") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size=15),
        axis.title = element_text(size = 18))

NB.Unknown <- countObj.Annot.NoCL.totalReads %>% filter(grepl('NB.Unknown', DIAGNOSIS.Substatus.Tumor.Normal.Tissue))
NB.U <- ggplot(NB.Unknown[,c(1,3,5,6)]) +
  geom_step(aes(y = count, x = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp,colour=Color.Substatus), size=0.6 ) +
  theme_bw() +
  theme(legend.position="none") +
  ggtitle(label = "NB.Unknown") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size=15),
        axis.title = element_text(size = 18))

pdf("IGH.color.step.noFilter.NB.v4.pdf", height = 10, width = 25)
cowplot::plot_grid(NB.A, NB.NA, NB.U, nrow = 1, align = 'h', axis = 'h')
dev.off()
############################################ Section 2 Aggregating and summarizing ###############################################

### Summarize By Samples ####
countObj.Annot.NoCL <- countObj.Annot %>% filter(!grepl('CellLine|Normal',LIBRARY_TYPE)) %>% filter(!grepl('^NS', DIAGNOSIS.Substatus.Tumor.Normal.Tissue) )
countObj.Annot.NoCL <- countObj.Annot.NoCL[complete.cases(countObj.Annot.NoCL),]
countObj.Annot.gb.Samples <- countObj.Annot.NoCL  %>% dplyr::group_by(Sample.Biowulf.ID.GeneExp) %>% 
  dplyr::summarise(
    TotalClones=n(),
    TotalCloneSum=sum(count),
    Diagnosis= paste(unique(DIAGNOSIS.Substatus.Tumor.Normal.Tissue), collapse = ',' )
  )  %>% 
  dplyr::rename_(.dots=setNames(list("TotalClones"),c(cloneType))); dim(countObj.Annot.gb.Samples); head(countObj.Annot.gb.Samples); tbl_df(countObj.Annot.gb.Samples)
countObj.Annot.gb.Samples.Annotate  <- left_join(countObj.Annot.gb.Samples, metaData[,c("Sample.Biowulf.ID.GeneExp","Sample.ID.Alias","LIBRARY_TYPE","DIAGNOSIS.Substatus.Tumor.Normal.Tissue")], 
                                                 by="Sample.Biowulf.ID.GeneExp")
###Replace NA by 0
countObj.Annot.gb.Samples.Annotate[which(is.na(countObj.Annot.gb.Samples.Annotate$TotalCloneSum)), c(cloneType,"TotalCloneSum")] <- 0
countObj.Annot.gb.Samples.Annotate[which(countObj.Annot.gb.Samples.Annotate$TotalCloneSum == 0), "IGHClones"] <- 0

### Saving files
saveRDS(countObj.Annot.gb.Samples, paste("./Results/countObj.Annot.gb.Samples.Annotate.v2",".", cloneType,".RDS" , sep="") )
write.table(countObj.Annot.gb.Samples, paste("./Results/countObj.Annot.gb.Samples.v2.",".", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)
write.table(countObj.Annot.gb.Samples.Annotate, paste("./Results/countObj.Annot.gb.Samples",".Annotate.v2.", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)

### Summarize By Diagnosis ####
countObj.Annot.gb.Diagnosis  <- countObj.Annot  %>% dplyr::group_by(cdr3aa, v, d, j, DIAGNOSIS.Alias) %>% 
  summarise(
    TotalSamples=n(), 
    Samples = paste(Sample.ID, collapse = ','), 
    CloneCount = paste(count, collapse = ','),
    MedianCloneCount = median(count)
  )
# saveRDS(countObj.Annot.gb.Diagnosis, "./Results/countObj.Annot.gb.Diagnosis.RDS")
# write.table(countObj.Annot.gb.Diagnosis, "./Results/countObj.Annot.gb.Diagnosis.txt", sep="\t", quote = F, row.names = F)

############################################ Section 3 Plots #####################################################################

###### plot for Correlation between ImmuneEnrichment Scores vs Total clones ####

# Read the enrichment score data ####
# ssGSEA <- read.table("C:/Users/sindiris/R Scribble/RNASeq.RSEM/GSEA/results/RPKM_Data_Filt.NoCLNS.zscore.meta.Broad.PROJ.txt",
#                      sep = "\t", header = T, row.names = 1, check.names = FALSE)
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("C:/Users/sindiris/R Scribble/RNASeq.RSEM/GSEA/New analysis/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.gct")
ssGSEA.zscore <- apply(ssGSEAScores, 1, zscore_All) ; ssGSEA.t <- ssGSEA.zscore %>% data.frame() %>% tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")

# Binding ImmuneScore with TCR ####
countObj.Annot.gb.Samples.Annotate.NoNS <- countObj.Annot.gb.Samples.Annotate %>% filter( ! LIBRARY_TYPE %in% c("Normal", "CellLine")) %>% dplyr::mutate(TotalCloneSum = log10(TotalCloneSum+1)) 
length(unique(countObj.Annot.gb.Samples.Annotate.NoNS$Diagnosis))
immuneScore.Clones <- left_join(countObj.Annot.gb.Samples.Annotate.NoNS[,c("Sample.Biowulf.ID.GeneExp", "TotalCloneSum", "Diagnosis")], ssGSEA.t, by="Sample.Biowulf.ID.GeneExp") %>% 
                      data.frame
#immuneScore.Clones %<>% tibble::column_to_rownames("Sample.Biowulf.ID.GeneExp")
#immuneScore.Clones$TotalCloneSum <- log10(immuneScore.Clones$TotalCloneSum + 1)

# plot & Save   ####

varNames <- colnames(immuneScore.Clones[,c(2,4:26)])  
plotLists <- lapply(varNames, correlationPlots, constName="TotalCloneSum",  df= data.frame(immuneScore.Clones))
ImmuneScorePlots <- lapply(plotLists, function(l) l[[1]] )

SBName =paste("ImmuneScore.vs.TotalCloneSum.v2",cloneType,"new.pdf",sep=".")
ggsave(SBName, marrangeGrob(ImmuneScorePlots, ncol=1, nrow=1), width = 15, height = 10)
dev.off()

## Save files
##write.table(immuneScore.Clones, paste("./Results/immuneScore",".", cloneType,".txt" , sep=""), sep="\t", quote = F, row.names = F)

###### plot for TCR COunt Bean plot ####

## Prepare data for one variable plot ####

selectCol="TotalCloneSum" ; StatsFinalCol="Diagnosis" ; SampleNames <- "Sample.ID.Alias"
tcrcloneCountPre          <- countObj.Annot.gb.Samples.Annotate.NoNS %>% 
                             dplyr::select_(.dots=c(paste0("selectCol"), paste0("StatsFinalCol"), paste0("SampleNames")))

tcrcloneCountPre.Diag   <- tcrcloneCountPre %>%  dplyr::rename_(.dots = setNames(list(SampleNames,StatsFinalCol),c("Samples","Diagnosis"))) 
ScoresPre               <- tcrcloneCountPre.Diag[,!(colnames(tcrcloneCountPre.Diag) %in% c("Samples")), drop=FALSE]
orderOfFactor           <- as.character( unique(ScoresPre$Diagnosis) )
orderOfSignature        <- colnames(ScoresPre)[-ncol(ScoresPre)]
colList                 <- c(1:(ncol(ScoresPre)-1)) ; Scores <- ScoresPre

## Plot and Save ####
plotLists <- OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =FALSE, yLab = paste0("Log clones (", cloneType, ")"),
                                 summaryHlines = T, sizeOfDots = 0.6)
tcrcloneCountPlots <- lapply(plotLists, function(l) l[[1]])
tcrcloneCountData  <- lapply(plotLists, function(l) l[[2]]) %>% data.frame(check.names = FALSE) %>% bind_rows()

SBName = filename=paste(date,selectCol,cloneType, "2.new.pdf",sep=".")
ggsave(SBName, marrangeGrob(tcrcloneCountPlots, ncol=1, nrow=1), width = 15, height = 10)
dev.off()

###### plot for percent vs cloneCopy ####

## Remove samples with no cdr3aa ####

countObj.Annot.NoNA <- countObj.Annot %>% dplyr::filter( count != 0) %>% dplyr::filter( ! LIBRARY_TYPE %in% c("Normal", "CellLine")) ; dim(countObj.Annot.NoNA)
countObj.Annot.PercentTCR <- countObj.Annot.NoNA  %>% dplyr::group_by(Sample.ID) %>% 
  dplyr::mutate( TotalCount = sum(count), percentinSample = (count/TotalCount)) %>% 
  dplyr::select(count, cdr3aa, DIAGNOSIS.Alias, TotalCount, percentinSample)
dim(countObj.Annot.PercentTCR) ; #View(countObj.Annot.PercentTCR)

## Plot and Save ####

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log2(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}

pctPlot <- ggplot( data = countObj.Annot.PercentTCR, aes( count, percentinSample) ) + 
  geom_point(aes(colour = factor(DIAGNOSIS.Alias)), size = 2.5) + 
  coord_trans(x="log10") +
  scale_colour_manual("Diagnosis", values=setNames( StatsFinal$Color, StatsFinal$Diagnosis )  )+
  scale_y_continuous( trans = log_trans(10), 
                      name = paste0(cloneType," clone percentage in each sample"),
                      breaks = c(0.01,0.1,0.25,0.50,0.75,1),
                      labels = scales::percent
  ) +
  scale_x_continuous( #trans = log_trans(10), 
    name =  paste0(cloneType," clone copies"),
    breaks = c(1,10,25,50,100,150,200)
  ) +
  theme_bw() +
  theme( panel.grid.major = element_line(colour = "grey50", size = 0.25), 
         panel.grid.minor = element_blank())  #element_line(colour = "grey50", size = 0.25) ) + 

SBName =paste(date,"Clone abundance vs expansion",cloneType,"new.pdf",sep=".")
ggsave(SBName, marrangeGrob(list(pctPlot), ncol=1, nrow=1), width = 15, height = 10)
dev.off()

######## Comparing TCR vs BCR counts ####

countObj.TCR <- countObj.Annot.gb.Samples.Annotate %>% filter( !Diagnosis %in% "NS")
countObj.IGH <- countObj.Annot.gb.Samples.Annotate %>% filter( !Diagnosis %in% "NS")

combined.IGH.TCR <- left_join(countObj.TCR, countObj.IGH, by=c("SAMPLE_ID", "DIAGNOSIS.Alias", "Diagnosis", "SAMPLE_ID.Alias"));
head(combined.IGH.TCR)

combined.IGH.TCR$TRBClones <- combined.IGH.TCR$TRBClones + 1
combined.IGH.TCR$IGHClones <- combined.IGH.TCR$IGHClones + 1

combined.IGH.TCR$TotalCloneSum.x <- combined.IGH.TCR$TotalCloneSum.x + 1
combined.IGH.TCR$TotalCloneSum.y <- combined.IGH.TCR$TotalCloneSum.y + 1

## Plot and Save ####

pctPlot <- ggplot( data = combined.IGH.TCR, aes( log10(TotalCloneSum.x), log10(TotalCloneSum.x) )) + 
  geom_point(aes(colour = factor(DIAGNOSIS.Alias)), size = 2.5) + 
  geom_smooth(method=lm,  fill="grey") +
  #coord_trans(x="log10", y="log10") +
  scale_colour_manual("Diagnosis", values=setNames( StatsFinal$Color, StatsFinal$Diagnosis )  )+
  theme_bw() +
  theme( panel.grid.major = element_line(colour = "grey50", size = 0.25), 
         panel.grid.minor = element_blank()) +
  xlab("TRB Total Clones ")+
  ylab("IGH Total Clones ")
  


