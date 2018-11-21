### TASK 1

## load the file
downstreamTargets <- read.csv("T:/Sivasish_Sindiri/Collaboration/TCR.Grants/NCI_GeneSet_RMS and EWS.txt", sep = "\t", header = TRUE,
                              stringsAsFactors = FALSE)
expressionTMM.RPKM <- readRDS("C:/Users/sindiris/R Scribble/RNASeq.RSEM/GeneRDSOutput/RPKM/RPKM_Data_Filt_Consolidated.GeneNames.2018-11-14.rds")
metadata           <- data.frame(lapply(rnaseqProject$validMetaDataDF,as.character), stringsAsFactors =FALSE)

metadata.RMS_FP.EWS <- metadata %>% dplyr::filter(DIAGNOSIS.Substatus.Tumor.Tissue %in% c("RMS.FP.Tumor","EWS.Tumor"))
metadata.RMS_FP.EWS.MYCN.A.Spleen.Tumor <- metadata %>% dplyr::filter(DIAGNOSIS.Substatus.Tumor.Normal.Tissue %in% c("RMS.FP","EWS", "NB.MYCN.A", "NS.spleen") & 
                                                           LIBRARY_TYPE %in% c("Tumor", "Normal") )
metadata.RMS_FP.EWS.MYCN.A.Spleen.CellLine <- metadata %>% dplyr::filter(DIAGNOSIS.Substatus.Tumor.Normal.Tissue %in% c("RMS.FP","EWS", "NB.MYCN.A", "NS.spleen") & 
                                                                  LIBRARY_TYPE %in% c("CellLine", "Normal") )

annotationRDS           = readRDS("C:/Users/sindiris/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS")
pcRDS                   = readRDS("C:/Users/sindiris/R Scribble/Annotation RDS/pc.other.HGNCTableFlat.rds")
tfRDS                   = readRDS("C:/Users/sindiris/R Scribble/Annotation RDS/TFs_no_epimachines.RDS")
csRDS                   = readRDS("C:/Users/sindiris/R Scribble/Annotation RDS/CellSurface.RDS")
cgaRDS                  = readRDS("C:/Users/sindiris/R Scribble/Annotation RDS/cancerGermlineAntigens.rds")

customGeneList <- list("tfRDS"=tfRDS, "csRDS"=csRDS,  "cgaRDS"=cgaRDS )

querryFoundList <- apply(downstreamTargets, 2, function(x){
    
    queryGenes = x
    foundList <- lapply(customGeneList, function(y){
                        
                        return(queryGenes[which(queryGenes %in% y$GeneName)])
              })
    return(foundList)
})

uniqueGenes <- unique(unlist(querryFoundList))
expressionTMM.RPKM.selected <- expressionTMM.RPKM %>% dplyr::filter(GeneName %in% uniqueGenes) %>% 
                                                      dplyr::select_(.dots=c("Chr","Start","End", "Strand","GeneID","GeneName", "Length",
                                                                                as.character(metadata.RMS_FP.EWS$SAMPLE_ID.Alias)))
stopifnot(nrow(expressionTMM.RPKM.selected)-6 != length(uniqueGenes))

write.table(expressionTMM.RPKM.selected, "T:/Sivasish_Sindiri/Collaboration/TCR.Grants/expressionTMM.RPKM.selected.txt", sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

## Ways of writing list to the file
## Method 1
lapply(unlist(querryFoundList, recursive = FALSE), write, "T:/Sivasish_Sindiri/Collaboration/TCR.Grants/RMS_FP.EWS.downstreamTargets.txt", append=TRUE)
## Method 2
foundGenesFromGeneSet <- do.call(rbind, lapply(unlist(querryFoundList, recursive = FALSE), as.list))
write.table(foundGenesFromGeneSet, "T:/Sivasish_Sindiri/Collaboration/TCR.Grants/RMS_FP.EWS.downstreamTargets.list.txt", sep = "\t", col.names = TRUE,
            row.names = TRUE, quote = FALSE)


### TASK 3

expressionTMM.RPKM_HLA           <- expressionTMM.RPKM
rownames(expressionTMM.RPKM_HLA) <- expressionTMM.RPKM$GeneName
## Toggle between cellline or tumor
expressionTMM.RPKM_HLA           <- expressionTMM.RPKM_HLA %>%  dplyr::select_(.dots=as.character(metadata.RMS_FP.EWS.MYCN.A.Spleen.Tumor$SAMPLE_ID.Alias))
#expressionTMM.RPKM_HLA           <- expressionTMM.RPKM_HLA %>%  dplyr::select_(.dots=as.character(metadata.RMS_FP.EWS.MYCN.A.Spleen.CellLine$SAMPLE_ID.Alias))
cytolyticScore                   <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM_HLA)
expressionTMM.RPKM_HLA_Cyto      <- t(rbind(expressionTMM.RPKM_HLA[c("HLA-A", "HLA-B", "HLA-C"),], cytolyticScore))
AliasNames_df  <- dplyr::left_join( data.frame("SAMPLE_ID.Alias"=rownames(expressionTMM.RPKM_HLA_Cyto)), 
                                    metadata.RMS_FP.EWS.MYCN.A.Spleen.Tumor[,c("SAMPLE_ID.Alias", "DIAGNOSIS.Substatus.Tumor.Normal.Tissue"),],
                                    #metadata.RMS_FP.EWS.MYCN.A.Spleen.CellLine[,c("SAMPLE_ID.Alias", "DIAGNOSIS.Substatus.Tumor.Normal.Tissue"),],
                                      by="SAMPLE_ID.Alias")

#### Zscore Function  
zscore_All = function( x = NA) {
medX <- median(x)
sdX <- sd(x, na.rm = FALSE)
y <- (x - medX) / sdX
return(y)
}

#### factorize columns
factorizeColumn = function( toFactor, asFactor){
  factorColumn <- factor(toFactor, levels=asFactor, ordered = TRUE )
  return(factorColumn)
}

#### One variable Plot
OneVariablePlotSort <- function(colList, Scores, orderOfFactor, orderOfSignature, standardize=FALSE, logit =FALSE,
                                plotType="StringBean",customColorDF=NA, yLab="Score", summaryHlines =FALSE, 
                                sizeOfDots = 1, legendDisplay=TRUE){
  
  #if (unique(is.na(customColors))) { customColors = setNames( StatsFinal$Color, StatsFinal$Diagnosis) }
  #function
  seqfunc   <- function(x, start){ return(seq(start, start+1, length.out = x))}
  xaxisSeq  <- function(x) {
    counts                  <- x %>% dplyr::group_by(Diagnosis) %>% dplyr::summarise(n= n())
    xvals                   <- c()
    for(count in counts$n){ 
      xvalsNew <- seqfunc(count, start = 1);
      xvals <- c(xvals, xvalsNew)
    }
    return(xvals)
  }
  drawStringBeanPlot <- function(x, tidyScoresPre){
    
    print(paste(x))
    tidyScores             <- tidyScoresPre %>% filter(Signatures == x) %>%  dplyr::group_by(Signatures,Diagnosis) %>% 
      dplyr::mutate(Med=median(Scores)) %>% 
      arrange(Signatures,Diagnosis,Scores) %>% 
      arrange(desc(Med)) %>% 
      ungroup() %>% 
      #mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor ),
      mutate( Diagnosis = factorizeColumn(Diagnosis, unique(Diagnosis) ),
              Color =  factorizeColumn(Color, unique(as.character(Color) ) ) ) %>% arrange(Diagnosis)
    tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)
    
    ##Make median Segment
    medianY <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianY=median(Scores))  %>% dplyr::arrange(Diagnosis,Signatures))$medianY
    medianX <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianX=median(SNONorm)) %>% dplyr::arrange(Diagnosis,Signatures))$medianX
    segmentDF <- data.frame( xstart = medianX-0.05, ystart=medianY, xend=medianX+0.15, yend=medianY)
    segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),Signatures=unique(tidyScores$Signatures)))
    
    summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
    scoreSummary <- summary(tidyScores$Scores)
    
    plot <- ggplot() +
      geom_point(data=tidyScores, aes(SNONorm, Scores, colour = factor(Diagnosis) ),show.legend = F, size=sizeOfDots) +
      scale_colour_manual(values=Color  ) +
      #geom_violin(data=tidyScores, aes(SNONorm, Scores, color = factor(Diagnosis)),show.legend = F)+
      facet_grid(Signatures~Diagnosis ,switch = "both") +
      #ylim( min(summaryStats$minV)-0.05,max(summaryStats$maxV)+0.05) + 
      #ylim(-3,2.5) +
      labs( title= x ) +
      ylab( yLab ) +
      theme_bw() +
      geom_hline(yintercept=0, size=0.1) +
      theme( title = element_text(size=13, face="bold")
             ,axis.title.x = element_blank()
             ,axis.title.y=element_text(size=10, face="bold")
             ,axis.text.x = element_blank()
             ,axis.text.y = element_text(size=10, face="bold")
             ,axis.ticks.x =element_blank()
             ,strip.text.y= element_blank()
             ,strip.text.x=element_text(size=10,face="bold", angle=90, vjust=1)
             ,strip.background=element_blank()
             ,panel.grid.major.x=element_blank()
             ,panel.grid.minor.x=element_blank()
             ,panel.border = element_rect(colour = "black", fill=NA, size=0.0000000002, linetype = 1)
             ,panel.spacing = unit(0, "cm")
             ,strip.switch.pad.grid = unit(0, "cm")
      ) +
      geom_segment(data = segmentDF, aes(x = xstart, xend = xend, y = ystart, yend = yend), size=1.5
                   #, linetype=2
                   , colour="red"
                   ,inherit.aes=FALSE
      ) + 
      scale_y_continuous( expand = c(0, -0.05),
                          limits = c(min(tidyScores$Scores)-0.25, max(tidyScores$Scores)+0.25)) +
      scale_x_continuous( expand = c(0.1, 0))
    
    
    if(summaryHlines) { plot <- plot + 
      geom_hline(yintercept = scoreSummary[[2]], linetype="dashed", colour="#8888ff", size=0.4) +
      geom_hline(yintercept = scoreSummary[[4]], linetype="dashed", colour="#8888ff", size=0.4) + 
      geom_hline(yintercept = scoreSummary[[5]], linetype="dashed", colour="#8888ff", size=0.4)
    }     
    colnames(segmentDF)[3] <- x ; segmentDF$Diagnosis <- factor(segmentDF$Diagnosis, levels = orderOfFactor, ordered = TRUE);
    segmentDF <- segmentDF %>% arrange(Diagnosis)
    return(list(plot, segmentDF[,c(1,3)]) )
  }
  drawDensityPlot <- function(x, tidyScoresPre=NA , orderOfFactor=NA, customColors=NA, yLab =yLab){
    
    print(paste(x))
    tidyScores             <- tidyScoresPre %>% filter(Signatures == x) %>%  dplyr::group_by(Signatures,Diagnosis) %>% 
      dplyr::mutate(Med=median(Scores)) %>% arrange(Signatures,Diagnosis,Scores) %>% 
      arrange(desc(Med)) %>% 
      ungroup() %>%
      mutate( Diagnosis = factorizeColumn(Diagnosis, unique(Diagnosis) ) ) %>% arrange(Diagnosis)
      #mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor )) %>% arrange(Diagnosis)
    tidyScores[,"SNONorm"] <- xaxisSeq(tidyScores)
    
    ##Make median Segment
    medianY <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianY=median(Scores)) %>% dplyr::arrange(Diagnosis,Signatures))$medianY
    medianX <- (tidyScores %>% dplyr::group_by(Diagnosis, Signatures) %>% dplyr::summarise(medianX=median(SNONorm))  %>% dplyr::arrange(Diagnosis,Signatures))$medianX
    segmentDF <- data.frame( xstart = medianX-0.05, ystart=medianY, xend=medianX+0.15, yend=medianY)
    segmentDF <- cbind(segmentDF, expand.grid(Diagnosis=unique(tidyScores$Diagnosis),Signatures=unique(tidyScores$Signatures)))
    
    summaryStats <- tidyScores %>% group_by(Diagnosis) %>% summarise(maxV = max(Scores), minV =min(Scores) )
    
    plot <- ggplot(data=tidyScores, aes(x = Scores, y = Diagnosis, height = ..density..)) +
      # to avoid overlaps of mountains , rel_min_height = 0.005, scale=0.9
      geom_density_ridges2(aes(fill = Diagnosis)) +
      scale_fill_manual(values=Color, guide=FALSE) +
      geom_vline(data=tidyScores, mapping=aes(xintercept=0), linetype = "dashed", colour = "maroon", size=1 ) +
      scale_y_discrete(expand = c(0.01, 0), limits = unique(rev(tidyScores$Diagnosis))) +
      scale_x_continuous(expand = c(0.01, 0)) +
      theme_ridges() + 
      theme(legend.title = element_text(size=15, face="bold") ) +
      labs( title= x ) +
      ylab("") +
      xlab(yLab)
    
    colnames(segmentDF)[3] <- x ; segmentDF$Diagnosis <- factor(segmentDF$Diagnosis, levels = orderOfFactor, ordered = TRUE);
    segmentDF <- segmentDF %>% arrange(Diagnosis)
    return(list(plot, segmentDF[,c(1,3)]) )
  }
  
  if(logit == TRUE) {
    Scores[,colList] <- apply(Scores[,colList,drop=FALSE] + 1 , 2, log2 )
  }
  
  if(standardize==TRUE){
    Scores[,colList] <- apply(Scores[,colList,drop=FALSE], 2, zscore_All )
  }
  
  tidyScoresPre <- Scores %>% tidyr::gather(Signatures, Scores, colList);
  mergeDF       <-  merge(tidyScoresPre, customColorDF, by.x="Diagnosis", by.y="Diagnosis", all.x=TRUE)
  tidyScoresPre  <- mergeDF[,c(1:4)] ; # tidyScoresPre$Diagnosis <- factor(tidyScoresPre$Diagnosis, levels = orderOfFactor, ordered = TRUE)
  
  if( plotType =="StringBean") {
    plotLists <- lapply(orderOfSignature, drawStringBeanPlot, tidyScoresPre)
  } else {
    customColorsVector <- setNames(as.character(customColorDF$Diagnosis), customColorDF$Color)
    plotLists <- lapply(orderOfSignature, drawDensityPlot, tidyScoresPre=tidyScoresPre, orderOfFactor=orderOfFactor, customColors=customColorsVector,
                        yLab =yLab)
  }
  return(plotLists)
}   

#### Plot the Immune Scores
orderOfFactor <- c("EWS", "RMS.FP", "NB.MYCN.A", "NS.spleen" )
Scores <- cbind(expressionTMM.RPKM_HLA_Cyto, AliasNames_df[,"DIAGNOSIS.Substatus.Tumor.Normal.Tissue", drop=FALSE]) %>% 
                  dplyr::rename(Diagnosis=DIAGNOSIS.Substatus.Tumor.Normal.Tissue) %>% 
                  dplyr::arrange(Diagnosis) %>% 
                  dplyr::mutate(Diagnosis = factor(Diagnosis, ordered = TRUE, levels = orderOfFactor))
orderOfSignature <- colnames(Scores)[-ncol(Scores)]
colList <- c(1:(ncol(Scores)-1))
#Generate custom colors
DiagFreq <- table(factor(Scores$Diagnosis,levels=orderOfFactor))
Color <- c("#990033", "#b36b00", "#0086b3", "#248f24")
customColorDF <- data.frame("Diagnosis"=names(DiagFreq), "Color"=Color)
#customColors <- unlist(sapply(seq(1:length(DiagFreq)), function(x){rep(customColors[x],DiagFreq[x])}))

plotLists <- OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =TRUE, logit =TRUE, 
                                 yLab = "Standardised gene expression", legendDisplay = FALSE, customColorDF = customColorDF,sizeOfDots = 1.5 )
ExpressionScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName ="T:/Sivasish_Sindiri/Collaboration/TCR.Grants/TMM.RPKM.GP.log2.zscore.RMS.EWS.MYCN.Spleen.Tumor.pdf"
ggsave(SBName, marrangeGrob(ExpressionScorePlots,ncol=4,nrow=1 ), width = 20, height = 8 )



### TASK 2
colNamesTumorsNormal.TCGA <- read.table("C:/Users/sindiris/R Scribble/RNASeq.RSEM/MetadataMapper.khanlab.TCGA.txt",
                                        sep = "\t", header = T, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
colNamesTumorsNormal.TCGA  <- colNamesTumorsNormal.TCGA  %>% arrange(desc(Source)) %>% dplyr::filter(!LIBRARY_TYPE %in% c("CellLine"))

ssGSEA <- read.table("C:/Users/sindiris/R Scribble/RNASeq.RSEM/GSEA/results/txt/RPKM_Data_Filt.log2.meta.Broad.PROJ.gct.txt",sep = "\t", 
                     header = T, row.names = 1, check.names = FALSE)

## Checking metadata vs data ###
metaData.colnames <- colNamesTumorsNormal.TCGA$SAMPLE_ID.Alias
dataMatrix.colnames <- colnames(ssGSEA)
TCGA.Colnames[which(!dataMatrix.colnames %in% metaData.colnames )]
## Check complete

selected.khanlab.TCGA.metadata  <- colNamesTumorsNormal.TCGA  %>% filter( DIAGNOSIS.Alias %in% c("LAML", "LUAD", "SKCM", "KIRC", "EWS", "RMS") )   %>% 
  dplyr::select( SAMPLE_ID.Alias, DIAGNOSIS.Alias )

ssGSEA.Selected <- ssGSEA %>% dplyr::select_(.dots = selected.khanlab.TCGA.metadata$SAMPLE_ID.Alias)

## Checking metadata vs data ##
stopifnot( ncol(ssGSEA.Selected) == length(as.character(selected.khanlab.TCGA.metadata$SAMPLE_ID.Alias)) )
#write.table(ssGSEA.Selected, "T:/Sivasish_Sindiri/Collaboration/TCR.Grants/ssGSEA.Selected.txt", sep = "\t", col.names = TRUE,
#            row.names = TRUE, quote = FALSE)
## Checking complete

orderOfFactor <- c("LAML", "LUAD", "SKCM", "KIRC", "EWS", "RMS") 
Scores <- cbind(t(ssGSEA.Selected), selected.khanlab.TCGA.metadata[,"DIAGNOSIS.Alias", drop=FALSE]) %>% 
  dplyr::rename(Diagnosis=DIAGNOSIS.Alias) %>% 
  dplyr::arrange(Diagnosis) #%>%
  #dplyr::mutate(Diagnosis = factor(Diagnosis, ordered = TRUE, levels = orderOfFactor))

orderOfSignature <- colnames(Scores)[-ncol(Scores)]
colList <- c(1:(ncol(Scores)-1))
#Generate custom colors
DiagFreq <- table(factor(Scores$Diagnosis,levels=orderOfFactor))
Color <- c("#a6a6a6", "#a6a6a6", "#a6a6a6", "#a6a6a6", "#990033", "#b36b00")
customColorDF <- data.frame("Diagnosis"=names(DiagFreq), "Color"=Color)
#customColors <- unlist(sapply(seq(1:length(DiagFreq)), function(x){rep(customColors[x],DiagFreq[x])}))

plotLists <- OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =TRUE, logit =FALSE, plotType = "density",
                                 yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF,sizeOfDots = 1 )
ExpressionScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName ="T:/Sivasish_Sindiri/Collaboration/TCR.Grants/TMM.RPKM.GP.log2.ssGSEA.zscore.khanlab.density.TCGA.pdf"
ggsave(SBName, marrangeGrob(ExpressionScorePlots,ncol=3,nrow=1 ), width = 20, height = 8 )




