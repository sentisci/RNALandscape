setwd("T:/Sivasish_Sindiri/R_workspace/MiXCR/")

# setwd("T:/Sivasish_Sindiri/R Scribble/RNALandscape")

## Khanlab meta data
metaData <- read.csv("MetadataMapper.v3.txt", sep="\t")
# metaData <- read.csv("tcr_rnaseq_file.txt", sep="\t")

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
#fileList <- list.files("./tcr_rnaseq/")
AllClonesData             <- rbindlist( lapply(fileList, function(x){
  print(x)
  #exomeData <- read.csv( paste("./tcr_rnaseq/", x, sep=""), sep="\t", header = TRUE )
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

######################################################## ENTROPY Analysis ############################################
### write the inputs for entropy data #
# fileList <- list.files("./CloneFilesNitin/")
# AllClonesEntropyData             <- sapply(fileList, makeEntropyInput, inputDir="./CloneFilesNitin/", outputDir="./CloneFilesEntropyNitin/" )

fileList <- list.files("./CloneFiles.v2/")
cloneType = "TRB"
AllClonesEntropyData             <- sapply(fileList, makeEntropyInput,  cloneType=cloneType, inputDir="./CloneFiles.v2/", 
                                           outputDir=paste0("./CloneFilesEntropy.", cloneType, ".v2/") )

### write the inputs for ImmunoseqV2 entropy data 
fileList <- list.files("./immunoseqv2/")
ImmunoseqV2EntropyData             <- sapply(fileList, immunoseqv2 )

# ### List files and read data into a single data matrix for Entropy results #
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
# ### List files and read data into a single data matrix for Entropy results 
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

### List files and read data into a single data matrix for Entropy results
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
#write.table(AllEntropyData, "AllEntropyData_H_CL_JS.landscape.TRB.v3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

AllEntropyData$FileName <- gsub(".Entropy.tx","",AllEntropyData$FileName)
AllEntropyData %<>% dplyr::rename(Sample.Data.ID=FileName)
AllEntropyData.annot <- dplyr::full_join(AllEntropyData, metaData[,c("Sample.Data.ID", "Sample.ID.Alias",
                                                                     "LIBRARY_TYPE",
                                                                     "DIAGNOSIS.Substatus.Tumor.Normal.Tissue")], by="Sample.Data.ID")
  


############################################ Section 1b Select DF ################################################################
### To Do  ####
### Add switch case to select Clone df based on user input ###

#################################################################################### For RNASEq ################################################################
### For now Select DF manually ####
#cloneType = "IGHClones"  ; countObj <- cloneObjIG %>% as.data.frame()
cloneType = "TRBClones"  ; countObj <- cloneObjTCR %>% as.data.frame()

### Attach metadata and generate countObj ####
countObj <- countObj %>% dplyr::rename(Sample.Data.ID=SampleName); 
countObj$Sample.Data.ID <- gsub("convert.|.clones.txt","", countObj$Sample.Data.ID)
#countObj$SAMPLE_ID <- gsub("-","_", countObj$SAMPLE_ID)
countObj.Annot <- dplyr::left_join(countObj, metaData, by="Sample.Data.ID") %>% 
  dplyr::select(one_of("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart", "Sample.Biowulf.ID.GeneExp", "Sample.ID.Alias",
                         "LIBRARY_TYPE","DIAGNOSIS.Substatus.Tumor.Normal.Tissue", "Color.Jun" )) ; head(countObj.Annot)
### Plot the clone expansion
countObj.Annot.NoCL <- countObj.Annot %>% filter(!grepl('CellLine',LIBRARY_TYPE)) %>% filter(!grepl('^NS', DIAGNOSIS.Substatus.Tumor.Normal.Tissue) )
countObj.Annot.NoCL %<>% dplyr::rename(Diagnosis = DIAGNOSIS.Substatus.Tumor.Normal.Tissue)


countObj.Annot.complete <- countObj.Annot.NoCL[complete.cases(countObj.Annot.NoCL),]
## sanity check
dim(countObj)
dim(countObj.Annot)
dim(countObj.Annot.NoCL)
dim(countObj.Annot.complete)

countObj.Annot.NoCL.totalReads <- dplyr::left_join(countObj.Annot.complete, readCountsSum.df, by="Sample.Biowulf.ID.GeneExp")
countObj.Annot.NoCL.totalReads.complete <- countObj.Annot.NoCL.totalReads[complete.cases(countObj.Annot.NoCL.totalReads),]
dim(countObj.Annot.NoCL.totalReads.complete)

countObj.Annot.NoCL.totalReads$ReadsPerMillion <- ( countObj.Annot.NoCL.totalReads$count/countObj.Annot.NoCL.totalReads$readCountsSum)*1000000

#################################################################################### For TCRSeq ###################################
### For now Select DF manually 

#cloneType = "IGHClones"  ; countObj <- cloneObjIG %>% as.data.frame()
cloneType = "TRBClones"  ; countObj <- cloneObjTCR %>% as.data.frame()

### Attach metadata and generate countObj #
countObj <- countObj %>% dplyr::rename(Sample.Data.ID=SampleName); 
countObj$Sample.Data.ID <- gsub("convert.|.clones.txt","", countObj$Sample.Data.ID)
#countObj$SAMPLE_ID <- gsub("-","_", countObj$SAMPLE_ID)
countObj.Annot <- dplyr::left_join(countObj, metaData, by="Sample.Data.ID") %>% 
  dplyr::select_(.dots=c("count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "VEnd", "DStart", "DEnd", "JStart", "Sample.Biowulf.ID.GeneExp",
                         "Diagnosis", "Color.Substatus" )) ; head(countObj.Annot)
### Plot the clone expansion
countObj.Annot.NoCL <- countObj.Annot %>%  filter(!grepl('^NS|^NA', Diagnosis) )

countObj.Annot.complete <- countObj.Annot.NoCL[complete.cases(countObj.Annot.NoCL),]
## sanity check
dim(countObj)
dim(countObj.Annot)
dim(countObj.Annot.NoCL)
dim(countObj.Annot.complete)



######## Make Step Plots to show expansion ####
countObj.Annot.NoCL.totalReads <- countObj.Annot.NoCL.totalReads %>% dplyr::select(Sample.Biowulf.ID.GeneExp, 
                                                                                   Diagnosis, 
                                                                                   count, freq, readCountsSum, ReadsPerMillion,
                                                                                   Color.Jun)
                                                                                   #Color.Substatus)

#toPlotDF <- countObj.Annot.NoCL.totalReads %>% dplyr::mutate(ReadsPerMillion = if_else(ReadsPerMillion >= 2, 2, ReadsPerMillion))

toPlotDF <- countObj.Annot.NoCL.totalReads %>% dplyr::filter(count > 0) %>% 
  arrange(Sample.Biowulf.ID.GeneExp, count) %>%
  group_by(Sample.Biowulf.ID.GeneExp) %>% 
  mutate(rank = dense_rank( -count )) %>% 
  distinct()
  # mutate(good_ranks = order(order(order_values, decreasing=TRUE)))
View(toPlotDF)

val = c("NB.MYCN.NA", "ASPS", "HBL", "NB.Unknown", "RMS.FP", "RMS.FN", "NB.MYCN.A", "UDS", "OS", "EWS", "DSRCT", "SS", "CCSK", "ML", "WT", "YST", "Teratoma")
toPlotDF$Diagnosis <- factor(toPlotDF$Diagnosis,levels = val, ordered = TRUE)
toPlotDF %<>% dplyr::arrange(Diagnosis)
customColorsVector <- setNames( unique(as.character(toPlotDF$Color.Jun)), unique(as.character(toPlotDF$Diagnosis)) )
# %<>% dplyr::mutate(Color.Jun = factor(Color.Jun, levels = unique(Color.Jun), ordered = TRUE))
### For RNASeq in Reverse
pdf("RNASeq.color.jun.step.v6.pdf", height = 15, width = 25)
ggplot(toPlotDF[,c(1,3,5,6,7)]) +
  geom_step(aes(x = rank, y = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp, colour= as.character(toPlotDF$Diagnosis) ), 
            size=0.7 ) +
  scale_colour_manual(values=customColorsVector) +
  facet_wrap(~toPlotDF$Diagnosis) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.text=element_text(size=16, face = "bold"),
        axis.text = element_text(size=14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "skyblue", fill=NA, size=1)) +
  scale_x_continuous( breaks = seq(1,max(toPlotDF$rank),by=4) ) +
  coord_trans(y = "log2" ) +
  scale_y_continuous(minor_breaks = c(),
                     breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4),
                     labels = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4) 
  ) 
dev.off()

### For RNASeq in original
ggplot(toPlotDF[,c(1,3,5,6,7)]) +
  geom_step(aes(y = rank, x = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp,colour=Color.Jun), 
            size=0.7 ) +
  facet_wrap(~toPlotDF$Diagnosis) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.text=element_text(size=16, face = "bold"),
        axis.text = element_text(size=14, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(colour = "skyblue", fill=NA, size=1)) +
  scale_y_continuous(trans = "reverse", breaks = seq(1,max(toPlotDF$rank),by=4) ) +
  coord_trans(x = "log2" ) +
  scale_x_continuous(minor_breaks = c(),
                     breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4),
                     labels = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4)
  ) 


### For TCRSeq
pdf("TCRSeq.color.step.v4.pdf", height = 15, width = 25)
ggplot(toPlotDF[,c(1,3,5,6,7)]) +
  geom_step(aes(y = rank, x = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp,colour=Color.Jun), 
            size=0.7 ) +
  facet_wrap(~toPlotDF$Diagnosis) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.text=element_text(size=16, face = "bold"),
        axis.text = element_text(size=12, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.border = element_rect(colour = "skyblue", fill=NA, size=1)) +
  scale_y_continuous(trans = "reverse", breaks = seq(1,max(toPlotDF$rank),by=25) ) +
  coord_trans(x = "log10" ) +
  scale_x_continuous(minor_breaks = c()                   
                     #breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4),
                     #labels = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4) 
  )
dev.off()

pdf("IGH.color.step.v5.pdf", height = 15, width = 25)
ggplot(toPlotDF[,c(1,3,5,6,7)]) +
  geom_step(aes(y = log2(rank), x = ReadsPerMillion, group=Sample.Biowulf.ID.GeneExp,colour=Color.Jun), 
            size=0.7 ) +
  facet_wrap(~toPlotDF$Diagnosis) +
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.text=element_text(size=16, face = "bold"),
        axis.text = element_text(size=12, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.border = element_rect(colour = "skyblue", fill=NA, size=1)) +
  scale_y_continuous(name = "log rank",
                     trans = "reverse")
                     #,
                     #breaks = c(0,  2,  4,  6, 8),
                     #labels = c(1,  4,  16, 64, 256)) +
  coord_trans(x = "log2") +
  scale_x_continuous(minor_breaks = c(),
                     breaks = c(0.01, 0.04, 0.3, 2,10,50,200,600),
                     labels = c(0.01, 0.04, 0.3, 2,10,50,200,600) ) 
dev.off()

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

countObj.Annot.NoNA <- countObj.Annot.NoCL.totalReads %>% dplyr::filter( count != 0); dim(countObj.Annot.NoNA)
countObj.Annot.PercentTCR <- countObj.Annot.NoNA  %>% dplyr::group_by(Sample.Biowulf.ID.GeneExp) %>% 
  dplyr::mutate( TotalCloneSum = sum(count), percentinSample = freq) %>% 
  dplyr::select(count, cdr3aa, Diagnosis, TotalCloneSum, percentinSample,ReadsPerMillion)
dim(countObj.Annot.PercentTCR) ; #View(countObj.Annot.PercentTCR)

## Plot and Save ####

base_breaks <- function(n = 10){
  function(x) {
    axisTicks(log2(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
customColorsVector <- setNames( unique(as.character(countObj.Annot.NoNA$Color.Jun)), unique(as.character(countObj.Annot.NoNA$Diagnosis)))

pctPlot <- ggplot( data = countObj.Annot.PercentTCR, aes( ReadsPerMillion, percentinSample) ) + 
  geom_point(aes(colour = factor(Diagnosis)), size = 2.5) + 
  coord_trans(x="log10") +
  scale_colour_manual("Diagnosis", values=customColorsVector  )+
  scale_y_continuous( trans = log_trans(10), 
                      name = paste0("frequency of a TCRB clone"),
                      breaks = c(0.01,0.1,0.25,0.50,0.75,1),
                      labels = scales::percent
  ) +
  scale_x_continuous( #trans = log_trans(10), 
    name =  paste0(" Expression of each TCRB clone"),
    #breaks = c(0.1,0.5,1,5,10)
  ) +
  theme_bw() +
  theme( panel.grid.major = element_line(colour = "grey50", size = 0.25), 
         panel.grid.minor = element_blank())  #element_line(colour = "grey50", size = 0.25) ) + 

## Filtering based on Frequency of each clones and its expansion (percentage) in that sample
test <- countObj.Annot.PercentTCR %>% filter( count >= 10 & percentinSample >= 0.01 )
countOFSamples <- test %>% group_by(Diagnosis) %>% mutate(countOFSamples = length(unique(Sample.Biowulf.ID.GeneExp)),
                                                          CountOFCDR3 = n()) %>% dplyr::select(Diagnosis, countOFSamples, CountOFCDR3) %>% 
                                            distinct() %>% arrange(countOFSamples)
countOFSamples

## Filtering based on its expansion (percentage) in that sample and Expression quantile
Exp_95_pct_cuttoff <- quantile(countObj.Annot.PercentTCR$ReadsPerMillion, probs = 0.99)
test <- countObj.Annot.PercentTCR %>% filter( ReadsPerMillion >= Exp_95_pct_cuttoff & percentinSample >= 0.01 )
countOFSamples <- test %>% group_by(Diagnosis) %>% mutate(countOFSamples = length(unique(Sample.Biowulf.ID.GeneExp)),
                                                          CountOFCDR3 = n()) %>% dplyr::select(Diagnosis, countOFSamples, CountOFCDR3) %>% 
                                        distinct() %>% arrange(countOFSamples)

print(paste("Exp_95_pct_cuttoff ",  Exp_95_pct_cuttoff[[1]]))
countOFSamples

SBName =paste0("Clone.abundance.vs.expansion",cloneType,".Exoression.v3.pdf")
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


#### For Venn plots Public dataBases ####
## Group by TCR
### By TCR
#countObj.Annot.gb <- countObj.Annot.complete  %>% dplyr::group_by(cdr3aa, v, d, j) %>% 
countObj.Annot.gb <- countObj.Annot.complete  %>% dplyr::group_by(cdr3aa) %>% 
  dplyr::summarise(
    TotalSamples=n(), 
    Samples = paste(Sample.Biowulf.ID.GeneExp, collapse = ','), 
    CloneCount = paste(count, collapse = ','),
    MedianCloneCount = median(count),
    Diagnosis= paste(unique(DIAGNOSIS.Substatus.Tumor.Normal.Tissue), collapse = ',' )
  )

## Using VDJdb
vdjdb <- read.csv("./public/vdjdb-2018-01-17/vdjdb_all_CDR3aa.txt", sep="\t")
vdjdbHealthy <- vdjdb %>% filter(grepl("healthy",meta.subject.cohort))
vdjdb.CDR3Beta <- unique(as.character(vdjdbHealthy$cdr3.beta))
length(vdjdb.CDR3Beta)
vdjdb.CDR3Beta.Complete <- unlist(stringr::str_extract_all(vdjdb.CDR3Beta, "^C.*F")); 
length(vdjdb.CDR3Beta.Complete)

## Using Waren
warenetal <- read.csv("./warren_et_al.Sushma/warren_etal.blood.cdr3aa.txt", sep="\t", header = F)
warenetal.CDR3Beta <- unique(as.character(warenetal$V1))
length(warenetal.CDR3Beta)
warenetal.CDR3Beta.Complete <- unlist(stringr::str_extract_all(warenetal.CDR3Beta, "^C.*F")); 
length(warenetal.CDR3Beta.Complete)


## All TCGA ## Using Bo et al/TCGA 
TCGA.Tumor.Normal <- read.csv("./public/TCGA/Bo.et.al.CDR3.fa", sep="\t", header = F); dim(TCGA.Tumor.Normal)
TCGATumorNormalTab <- data.frame(Names=gsub(">","",TCGA.Tumor.Normal[seq(1, 1366836,2),]), Sequence=TCGA.Tumor.Normal[seq(2, 1366836,2),])
head(TCGATumorNormalTab); dim(TCGATumorNormalTab)

#Normal
TCGANormal <- TCGATumorNormalTab %>% dplyr::filter(grepl('^N', Names)) ; dim(TCGANormal)
#Tumor
TCGATumor  <- TCGATumorNormalTab %>% dplyr::filter(grepl('^T', Names)) ; dim(TCGATumor)

## Getting All the TCGA Tumor CDR3aa
TCGATumor.CDR3Beta.Complete <- unique(unlist(stringr::str_extract_all(TCGATumor$Sequence, "^C[A|S]+.*F")))
length(TCGATumor.CDR3Beta.Complete)
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
#countObj.Annot.gb.Cancer <- countObj.Annot.gb %>% filter(!Diagnosis %in% c("NS"))
countObj.Annot.gb.Cancer <- countObj.Annot.gb
countObj.Annot.gb.Cancer$length <- sapply(as.character(countObj.Annot.gb$cdr3aa), nchar)
cancerShared             <- countObj.Annot.gb.Cancer %>% filter(TotalSamples > 1); length(cancerShared$cdr3aa)
ccancerShared.Complete    <- unlist(stringr::str_extract_all(cancerShared$cdr3aa, "^C.*F")); length(cancerShared.Complete)
cancerPrivate            <- countObj.Annot.gb.Cancer %>% filter(TotalSamples == 1); length(cancerPrivate$cdr3aa)
cancerPrivate.Complete    <- unlist(stringr::str_extract_all(cancerPrivate$cdr3aa, "^C.*F")); length(cancerPrivate.Complete)

countObj.Annot.gb.NS <- countObj.Annot %>%  filter(grepl('^NS', DIAGNOSIS.Substatus.Tumor.Normal.Tissue) )
# countObj.Annot.gb.NS <- countObj.Annot.gb %>% filter(Diagnosis %in% c("NS")); length(countObj.Annot.gb.NS$cdr3aa)
inHouseNormal.Complete    <- unlist(stringr::str_extract_all(countObj.Annot.gb.NS$cdr3aa, "^C.*F")); length(inHouseNormal.Complete)

vennCDR3aaList <- list("Cohort, Private" = unique(paste(cancerPrivate$cdr3aa)),
                       "Cohort, Shared" = unique(paste(cancerShared$cdr3aa)),
                       "TCGA,  Tumor" = TCGATumorDistinct,
                       "TCGA,  Normal"    =   TCGANormalDistinct,
                       "Warren et al (healthy)" = warenetal.CDR3Beta, 
                       "Chudakov et al (healthy)"=vdjdb.CDR3Beta,
                       "In House Normal" = unique(paste(countObj.Annot.gb.NS$cdr3aa)) 
)

# vennCDR3aaList <- list("Tumor, Private" = unique(cancerPrivate.Complete),
#                        "Tumor, Shared" = unique(cancerShared.Complete),
#                        "Warreb et al (healthy)" = unique(warenetal.CDR3Beta.Complete), 
#                        "Chudakov et al (healthy)"=unique(vdjdb.CDR3Beta.Complete),
#                        "In House Normal" = unique(inHouseNormal.Complete) 
# )

library(venn)
v.table <- venn::venn(vennCDR3aaList, ilab=TRUE, zcolor = "style", size = 15, cexil = 1, cexsn = 1)
TumorPrivate <- attr(v.table,"intersections")[["Tumor, Private"]]; length(TumorPrivate)
TumorShared  <- attr(v.table,"intersections")[["Tumor, Shared"]]; length(TumorShared)
InHouseNormals  <- attr(v.table,"intersections")[["In House Normal"]]; length(InHouseNormals)
warenetalNormals  <- attr(v.table,"intersections")[["Warreb et al (healthy)"]]; length(warenetalNormals)


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

############################# GLIPH analysis ###################################################
gliph.Khanlab <- read.table("./gliph/GliphInput.khanlab-convergence-groups.v2.txt", sep="\t", header = T)
gliph.Khanlab %<>% arrange(-Column1) %<>% separate(Column2, c("Column2a", "Column2b"), sep="-")
gliph.Warren <- read.table("./gliph/GliphInput.Warren-convergence-groups.v2.txt", sep="\t", header = T)
gliph.Warren %<>% arrange(-Column1) %<>% separate(Column2, c("Column2a", "Column2b"), sep="-")

vennCDR3aaList <- list('khanlab'=as.character(gliph.Khanlab$Column2b),
                       'Waren'= as.character(gliph.Warren$Column2b))
vennCDR3aaList <- list('khanlab'=as.character('a','b','c'),
                       'Waren'= as.character('a','b','c'))
v.table <- venn::venn(vennCDR3aaList, ilab=TRUE, zcolor = "style", size = 15, cexil = 1, cexsn = 1)

############################# Disect consensus with respect to diagnosis ######################
gliph.Khanlab <- read.table("./gliph/GliphInput.khanlab-convergence-groups.v2.txt", sep="\t", header = T)
#countObj.Annot.complete
customColorsVector <- setNames( unique(as.character(toPlotDF$Color.Jun)), unique(as.character(toPlotDF$Diagnosis)) )
consesus <- lapply(seq(1:24579), function(x){
  first_consensus <- unlist(strsplit(as.character(gliph.Khanlab[x,3]),split = " "))
  length(first_consensus)-gliph.Khanlab[x,1]
  first_consensus_indexes <- which( as.character(countObj.Annot.complete$cdr3aa) %in% first_consensus )
  first_countObj.Annot.complete <- countObj.Annot.complete[first_consensus_indexes,]
  first_countObj.Annot.complete$ConsensusSeq <- as.character(gliph.Khanlab[x,2])
  return(first_countObj.Annot.complete[,c("Sample.Biowulf.ID.GeneExp","Diagnosis", "ConsensusSeq")])
})
consensus_final <- do.call(rbind, consesus)
consensus_final %<>% group_by(ConsensusSeq) %>% mutate(Count = n(), 
                                                       Samples=paste0(Sample.Biowulf.ID.GeneExp, collapse = ", "),
                                                       Samples_Count=length(unique(Sample.Biowulf.ID.GeneExp)),
                                                       Histology=paste0(Diagnosis, collapse = ", "),
                                                       Histology_Count=length(unique(Diagnosis)) ) %>% 
                      dplyr::select(-one_of(c("Sample.Biowulf.ID.GeneExp","Diagnosis"))) %>%
                      dplyr::distinct() %>%
                      arrange(-(Count))

write.table(consensus_final, "consensus_final.txt", sep="\t", row.names = FALSE, col.names = TRUE)
consensus_final$ConsensusSeq <- factor(consensus_final$ConsensusSeq, levels = unique(consensus_final$ConsensusSeq), ordered = TRUE)
g <- ggplot(consensus_final, aes(ConsensusSeq)) + 
     geom_bar(aes(fill = factor(Diagnosis))) +
     scale_fill_manual(values = customColorsVector) +
     theme(title = element_text(size=13, face="bold")
         ,axis.title.x = element_text(size=13, face="bold")
         ,axis.title.y = element_text(size=13, face="bold")
         ,axis.text.x = element_text(size=10, face="bold", angle=90, vjust=1))
plot(g)

## Pie chart plot
for_pie_data <- consensus_final[,c("Histology_Count"), drop= FALSE]; head(for_pie_data)
for_pie_dataDF <- data.frame(table(for_pie_data))

## Waffle plot
consensus_Count_waffle <- consensus_final %>% group_by(Histology_Count) %>% mutate(Sum_hist_count = sum(Count)) %>%
              dplyr::select(one_of(c("Histology_Count", "Sum_hist_count"))) %>% 
              arrange(desc(Sum_hist_count)) %>% 
              distinct()
totalCDR3s <- sum(consensus_Count_waffle$Sum_hist_count)
consensus_Count_waffle <- consensus_Count_waffle %>% dplyr::mutate(percentage = (Sum_hist_count/totalCDR3s)*100)
consensus_Count_waffle$Histology_Count_name <- paste0(consensus_Count_waffle$Histology_Count," ( ", 
                                                     signif(consensus_Count_waffle$percentage,2), "% )")
consensus_Count_waffle

sample_group <- as.numeric( consensus_Count_waffle$Sum_hist_count)
names(sample_group) <- as.character(consensus_Count_waffle$Histology_Count_name)

pdf("Waffle_plot_CDR3_sequence_sharing_30_v2.pdf")
colors <- c("#4594D0", "#E84A9A", "#808000", "#BD6354",
            "#A5CE39", "#FFC91D", "#4258A7", "#61C29E",
            "#9A6324", "#8E5C97", "#469990", "#f58231")
waffle(sample_group/30, size=0.05, colors = colors,
       title="Number of CDR3 sequences shared among histologies", 
       xlab="1 square == 30 CDR3 sequences")
dev.off()
















