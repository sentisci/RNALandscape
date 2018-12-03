expressionTMM.RPKM_HLA           <- expressionTMM.RPKM[,-c(1:11)]
rownames(expressionTMM.RPKM_HLA) <- expressionTMM.RPKM$GeneName.x
## Toggle between cellline or tumor
expressionTMM.RPKM_HLA           <- expressionTMM.RPKM_HLA %>%  dplyr::select_(.dots=as.character(rnaseqProject$metaDataDF$SAMPLE_ID.Alias))
#expressionTMM.RPKM_HLA           <- expressionTMM.RPKM_HLA %>%  dplyr::select_(.dots=as.character(metadata.RMS_FP.EWS.MYCN.A.Spleen.CellLine$SAMPLE_ID.Alias))
cytolyticScore                   <- corUtilsFuncs$cytolyticScore(expressionTMM.RPKM_HLA)
# expressionTMM.RPKM_HLA_Cyto      <- t(rbind(expressionTMM.RPKM_HLA[c("HLA-A", "HLA-B", "HLA-C"),], cytolyticScore))
# AliasNames_df  <- dplyr::left_join( data.frame("SAMPLE_ID.Alias"=rownames(expressionTMM.RPKM_HLA_Cyto)),
#                                     metadata.RMS_FP.EWS.MYCN.A.Spleen.Tumor[,c("SAMPLE_ID.Alias", "DIAGNOSIS.Substatus.Tumor.Normal.Tissue"),],
#                                     #metadata.RMS_FP.EWS.MYCN.A.Spleen.CellLine[,c("SAMPLE_ID.Alias", "DIAGNOSIS.Substatus.Tumor.Normal.Tissue"),],
#                                     by="SAMPLE_ID.Alias")

## Read ImmuneSignature data
geneSignatrue <- read.table("../Sondel.data/GSEA/results/only PC genes/RPKM_Data_Filt_ssGSEA.2018-11-29.ImmuneSignatures.gct",sep="\t", header = TRUE)
rownames(geneSignatrue) <- geneSignatrue$Name
geneSignatrueForPlot <- geneSignatrue[,-c(1)]
colnames(geneSignatrueForPlot) <- dplyr::left_join(data.frame(SAMPLE_ID=colnames(geneSignatrue)[-1]), 
                                            rnaseqProject$metaDataDF, by="SAMPLE_ID")[,"SAMPLE_ID.Alias"] %>% as.character()
geneSignatrueForPlotFinal <- rbind(geneSignatrueForPlot,cytolyticScore)
geneSignatrueForPlotFinalZscored <- t(apply(geneSignatrueForPlotFinal,1,zscore_All))
write.table(geneSignatrueForPlotFinalZscored, "../Sondel.data/GSEA/results/geneSignatrueForPlotFinalZscored.txt", sep = "\t", row.names = T, col.names = T, quote = FALSE)

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
      mutate( Diagnosis = factorizeColumn(Diagnosis, orderOfFactor ),
      #mutate( Diagnosis = factorizeColumn(Diagnosis, unique(Diagnosis) ),
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
  tidyScoresPre  <- mergeDF[,c(1:4)] ;  tidyScoresPre$Diagnosis <- factor(tidyScoresPre$Diagnosis, levels = orderOfFactor, ordered = TRUE)
  
  if( plotType =="StringBean") {
    plotLists <- lapply(orderOfSignature, drawStringBeanPlot, tidyScoresPre)
  } else {
    customColorsVector <- setNames(as.character(customColorDF$Diagnosis), customColorDF$Color)
    plotLists <- lapply(orderOfSignature, drawDensityPlot, tidyScoresPre=tidyScoresPre, orderOfFactor=orderOfFactor, customColors=customColorsVector,
                        yLab =yLab)
  }
  return(plotLists)
}   
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

#### Plot the Immune Scores
orderOfFactor <- c("ML.A", "ML.B" )
Scores <- cbind(t(geneSignatrueForPlotFinal), rnaseqProject$metaDataDF[,"DIAGNOSIS.Alias.substatus.T", drop=FALSE]) %>% 
  dplyr::rename(Diagnosis=DIAGNOSIS.Alias.substatus.T) %>% 
  dplyr::arrange(Diagnosis) %>% 
  dplyr::mutate(Diagnosis = factor(Diagnosis, ordered = TRUE, levels = orderOfFactor))
orderOfSignature <- colnames(Scores)[-ncol(Scores)]
colList <- c(1:(ncol(Scores)-1))
#Generate custom colors
DiagFreq <- table(factor(Scores$Diagnosis,levels=orderOfFactor))
Color <- c("#0086b3", "#248f24")
customColorDF <- data.frame("Diagnosis"=names(DiagFreq), "Color"=Color)
#customColors <- unlist(sapply(seq(1:length(DiagFreq)), function(x){rep(customColors[x],DiagFreq[x])}))

plotLists <- OneVariablePlotSort(colList, Scores=Scores, orderOfFactor, orderOfSignature, standardize =TRUE, logit =FALSE, 
                                 yLab = "Standardised enrichment score", legendDisplay = FALSE, customColorDF = customColorDF,sizeOfDots = 2.5 )
ExpressionScorePlots <- lapply(plotLists, function(l) l[[1]])
SBName ="../Sondel.data/Plots/Stringplot.v2.pdf"
ggsave(SBName, marrangeGrob(ExpressionScorePlots,ncol=2,nrow=1 ), width = 10, height = 8 )
