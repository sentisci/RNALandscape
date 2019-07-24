## Read the metadata 
metadata <- read.csv("../RNASeq.RSEM/MetadataMapper.v3.txt", sep="\t", header = T )

## STEP.1 Read the ssGSEA dataset
ssGSEAScores            <- corUtilsFuncs$parseBroadGTCOutFile("../RNASeq.RSEM/GSEA/New analysis/RPKM_Data_Filt_Consolidated.GeneNames.all.Khanlab.pc.log2.2019-06-14.PROJ.gct")
ssGSEAScores.t <- t(ssGSEAScores) %>% data.frame() %>% tibble::rownames_to_column(var="Sample.Biowulf.ID.GeneExp")

## STEP.2 Read the neoantigen dataset
neoantigen <- read.csv("../RNASeq.RSEM/NeoantigenFromSamples.Fusions.Variants.txt", sep="\t")

## STEP.3 Read the TCR count dataset
IGH.Count <- read.csv("../RNASeq.RSEM/TCR.Results/TCR.ResultscountObj.Annot.gb.Samples.IGH.txt", sep="\t")
#IGH.Count <- IGH.Count %>% dplyr::rename(Sample.Biowulf.ID=Sample.ID, IGH.TotalCloneSum=TotalCloneSum, IGH.UniqueClones=IGH )
IGH.Count <- IGH.Count %>% dplyr::rename(IGH.TotalCloneSum=TotalCloneSum, IGH.UniqueClones=IGH )

## STEP.4
TRB.Count <- read.csv("../RNASeq.RSEM/TCR.Results/TCR.ResultscountObj.Annot.gb.Samples.TRB.txt", sep="\t")
#TRB.Count <- TRB.Count %>% dplyr::rename(Sample.Biowulf.ID=Sample.ID, TRB.TotalCloneSum=TotalCloneSum, TRB.UniqueClones=TRB  )
TRB.Count <- TRB.Count %>% dplyr::rename(TRB.TotalCloneSum=TotalCloneSum, TRB.UniqueClones=TRB  )

## STEP.5 Read the Entropy
IGH.Entropy <- read.csv("../RNASeq.RSEM/TCR.Clones.Entropy/AllEntropyData_H_CL_JS.landscape.IGH.v3.txt", sep="\t")
IGH.Entropy <- IGH.Entropy %>% dplyr::rename(IGH.Entropy=Htot..Entropy., IGH.Richness=Num_CDR3..Richness., IGH.Clonality=Cltot..Clonality.)

## STEP.6
TRB.Entropy <- read.csv("../RNASeq.RSEM/TCR.Clones.Entropy/AllEntropyData_H_CL_JS.landscape.TRB.v3.txt", sep="\t")
TRB.Entropy <- TRB.Entropy %>% dplyr::rename(TRB.Entropy=Htot..Entropy., TRB.Richness=Num_CDR3..Richness., TRB.Clonality=Cltot..Clonality.)

## STEP.7 Read the fusion dataset
## Fusions
fusions <- read.csv("../RNASeq.Fusion.data/FinalFilteredfusionResultMatrix.txt", sep="\t")
fusions.Sample <- fusions %>% dplyr::group_by(Sample.ID) %>% dplyr::mutate(FusionsPerSample=n()) %>% 
                              dplyr::select(Sample.ID,FusionsPerSample ) %>% distinct(); View(fusions.Sample)
## STEP.8 Fusions for Neoantigens
fusionsForNeoantigens <- read.csv("../RNASeq.Fusion.data/FinalFusionResultMatrixForNeoantigens.txt", sep="\t")
fusionsForNeoantigens <- fusionsForNeoantigens %>% dplyr::group_by(Sample.ID) %>% dplyr::mutate(FusionsPerSample.Neoantiens=n()) %>% 
                                    dplyr::select(Sample.ID,FusionsPerSample.Neoantiens ) %>% distinct(); View(fusionsForNeoantigens)
 
## STEP.9 Read the mutation dataset
preFinalMutation <- read.csv("../RNASeq.Mutation.data/outputTXTOutput/4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10pc.Indels.LTE.1pc.v3.txt",
                             sep="\t")
FinalMutation <- read.csv("../RNASeq.Mutation.data/Fig1_manual_filter_tier1_minusfaultygermline_addselectTier2Tier3_indelback_05_21_19_manual_indel_splice_Selected.txt", sep="\t")

mutationForReporting <- dplyr::left_join( FinalMutation[,c("variantKey", "Patient.ID")], preFinalMutation[,c("variantKey", "Patient.ID", "Sample.ID")], 
                                       by=c("variantKey","Patient.ID") ); View(mutationForReporting)
ReportedPerSample <- mutationForReporting %>% dplyr::group_by(Sample.ID) %>% dplyr::mutate(variantsPerSampleFigure = n()) %>% 
                                  dplyr::select(Sample.ID, variantsPerSampleFigure) %>% 
                                  dplyr::distinct(); View(ReportedPerSample)
## STEP.10
preFinalReportedPerSample <- preFinalMutation %>% dplyr::group_by(Sample.ID) %>% dplyr::mutate(variantsPerSampleAll = n()) %>% 
                                                  dplyr::select(Sample.ID, variantsPerSampleAll) %>% 
                                                  dplyr::distinct(); View(preFinalReportedPerSample)

## STEP.11 Read the mutation for neoantigen
MutationForNeoantigens <- read.csv("../RNASeq.Mutation.data/outputTXTOutput/4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.10e4_propInTumor.LTE.10pc.Indels.LTE.1pc.Neoantigen.v3.txt",
                             sep="\t")
MutationForNeoantigensPerSample <- MutationForNeoantigens %>% dplyr::group_by(Sample.ID)  %>% dplyr::mutate(variantsPerSampleNeoantigens = n()) %>% 
                                                              dplyr::select(Sample.ID, variantsPerSampleNeoantigens) %>% 
                                                              dplyr::distinct(); View(MutationForNeoantigensPerSample)
#Merge all data together

## Merge ssGSEA
Step.1.Meta <- dplyr::left_join( metadata[,c("Sample.Data.ID","Sample.Biowulf.ID","Sample.Biowulf.ID.GeneExp","Sample.ID", "LIBRARY_TYPE", "DIAGNOSIS.Alias")], 
                                  ssGSEAScores.t, by="Sample.Biowulf.ID.GeneExp" )
indx <- apply(Step.1.Meta, 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.Meta[indx,1];View(Step.1.Meta[indx,c(1:6, ncol(Step.1.Meta))])

## Merge Neoantigen
Step.1.2.Meta <- dplyr::left_join( Step.1.Meta, neoantigen[,c("Sample.Biowulf.ID","VariantNeoAntigenCount","FusionNeoAntigenCount",
                                                                         "TotalNeoantigenCount")], by="Sample.Biowulf.ID" )
indx <- apply(Step.1.2.Meta, 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.2.Meta[indx,1];View(Step.1.2.Meta[indx,c(1:6, ncol(Step.1.2.Meta))])

## Merge IGH Count
Step.1.2.3.Meta <- dplyr::left_join( Step.1.2.Meta, IGH.Count[,c("Sample.ID", "IGH.TotalCloneSum", "IGH.UniqueClones")], by="Sample.ID")
indx <- apply(Step.1.2.3.Meta , 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.2.3.Meta [indx,1];View(Step.1.2.3.Meta [indx,c(1:6, ncol(Step.1.2.3.Meta ))])

## Merge TRB Count
Step.1.2.3.4.Meta <- dplyr::left_join( Step.1.2.3.Meta, TRB.Count[,c("Sample.ID", "TRB.TotalCloneSum", "TRB.UniqueClones")], by="Sample.ID")
indx <- apply(Step.1.2.3.4.Meta, 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.2.3.4.Meta[indx,1];View(Step.1.2.3.4.Meta[indx,c(1:6, ncol(ssScore.neoantigen.IGH.TRB.Meta))])

## Merge IGH Entropy
Step.1.2.3.4.5.Meta <- dplyr::left_join( Step.1.2.3.4.Meta, IGH.Entropy[,c("Sample.Data.ID", "IGH.Entropy", "IGH.Richness", "IGH.Clonality")], by="Sample.Data.ID")
indx <- apply(Step.1.2.3.4.5.Meta, 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.2.3.4.5.Meta[indx,1];View(Step.1.2.3.4.5.Meta[indx, c(1:6, ncol(Step.1.2.3.4.5.Meta))])

## Merge TRB Entropy
Step.1.2.3.4.5.6.Meta <- dplyr::left_join( Step.1.2.3.4.5.Meta, TRB.Entropy[,c("Sample.Data.ID", "TRB.Entropy", "TRB.Richness", "TRB.Clonality")], by="Sample.Data.ID")
indx <- apply(Step.1.2.3.4.5.6.Meta, 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.2.3.4.5.6.Meta[indx,1];View(Step.1.2.3.4.5.6.Meta[indx, c(1:6, ncol(Step.1.2.3.4.5.6.Meta))])

## Merge Fusions for Neoantigens
Step.1.2.3.4.5.6.7.Meta <- dplyr::left_join(Step.1.2.3.4.5.6.Meta, fusionsForNeoantigens, by="Sample.ID")
indx <- apply(Step.1.2.3.4.5.6.7.Meta, 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.2.3.4.5.6.7.Meta[indx,1];View(Step.1.2.3.4.5.6.7.Meta[indx, c(1:6, ncol(Step.1.2.3.4.5.6.7.Meta))])

## Merge Fusions
Step.1.2.3.4.5.6.8.Meta <- dplyr::left_join(Step.1.2.3.4.5.6.7.Meta, fusions.Sample, by="Sample.ID")
indx <- apply(Step.1.2.3.4.5.6.8.Meta, 1, function(x) any(is.na(x) | is.infinite(x)))
Step.1.2.3.4.5.6.8.Meta[indx,1];View(Step.1.2.3.4.5.6.8.Meta[indx, c(1:6, ncol(Step.1.2.3.4.5.6.8.Meta))])

## Merge mutations Final mutations
Step.1.2.3.4.5.6.8.9.Meta <- dplyr::left_join(Step.1.2.3.4.5.6.8.Meta, ReportedPerSample, by="Sample.ID")
Step.1.2.3.4.5.6.8.9.Meta[indx,1];View(Step.1.2.3.4.5.6.8.9.Meta[indx, c(1:6, ncol(Step.1.2.3.4.5.6.8.9.Meta))])

## Merge mutations prefinal Final mutations
Step.1.2.3.4.5.6.8.9.10.Meta<- dplyr::left_join(Step.1.2.3.4.5.6.8.9.Meta, preFinalReportedPerSample, by="Sample.ID")
Step.1.2.3.4.5.6.8.9.10.Meta[indx,1];View(Step.1.2.3.4.5.6.8.9.10.Meta[indx, c(1:6, ncol(Step.1.2.3.4.5.6.8.9.10.Meta))])

## Merge mutations for neoantigens
Step.1.2.3.4.5.6.8.9.10.11.Meta<- dplyr::left_join(Step.1.2.3.4.5.6.8.9.10.Meta, MutationForNeoantigensPerSample, by="Sample.ID")
Step.1.2.3.4.5.6.8.9.10.Meta[indx,1];View(Step.1.2.3.4.5.6.8.9.10.Meta[indx, c(1:6, ncol(Step.1.2.3.4.5.6.8.9.10.Meta))])


## print the matrix
write.table(Step.1.2.3.4.5.6.8.9.10.11.Meta, 
            "../RNASeq.RSEM/ssScore.neoantigen.IGH_C.TRB_C.IGH_E.TRB_E.Neo.Fusion.Final_Mut.preFinal_Mut.Meta.v6.txt", sep="\t")

## Make plots
mergeDataSet <- read.table("../RNASeq.RSEM/ssScore.neoantigen.IGH_C.TRB_C.IGH_E.TRB_E.Neo.Fusion.Final_Mut.preFinal_Mut.Meta.v7.txt", sep="\t",
                           header = TRUE)

## Count neoantigens per snv and indel
fileList <- list.files("../RNASeq.Mutation.data/mutation_neoantigen_files/")
countNeoantigesnAlt = function(x) {
  fileData <- read.table(paste0("../RNASeq.Mutation.data/mutation_neoantigen_files/",x), sep = "\t", header = TRUE)
  snv.Neoantigen = count(grepl(">",fileData$HGVSc))
  indel.Neoantigen = dim(fileData)[1] - snv.Neoantigen
  name = basename(gsub(".filtered.condensed.ranked.tsv","",x))
  return(list("Sample.Biowulf.ID"=name, "indel.Neoantigen"= indel.Neoantigen, "snv.Neoantigen"=snv.Neoantigen ))
}
countNeoantigesnList <- lapply(fileList, countNeoantigesnAlt)
countNeoantigesnDF <- data.table::rbindlist( countNeoantigesnList )

## Merge the above data with the merged dataset
mergeDataSet.CountNeo <- dplyr::full_join(mergeDataSet, countNeoantigesnDF, by= "Sample.Biowulf.ID")
mergeDataSet.CountNeo[is.na(mergeDataSet.CountNeo)] <- 0

## Keep data only for Tumor
mergeDataSet.T <- mergeDataSet.CountNeo %>% dplyr::filter(LIBRARY_TYPE == "Tumor")
# mergeDataSet.Neoantigens <- mergeDataSet.T %>% dplyr::select(Sample.Biowulf.ID, DIAGNOSIS.Alias, FusionNeoAntigenCount,
#                                                              VariantNeoAntigenCount, fs.Neoantigen, snv.Neoantigen)
mergeDataSet.Neoantigens <- mergeDataSet.T %>% dplyr::select(Sample.Biowulf.ID, DIAGNOSIS.Alias, FusionNeoAntigenCount,
                                                             indel.Neoantigen, snv.Neoantigen)
dataSetNeo <- tidyr::gather(mergeDataSet.Neoantigens, "Alteration", "Count", indel.Neoantigen, snv.Neoantigen, FusionNeoAntigenCount )
View(dataSetNeo)

dataSetNeo1 <- dataSetNeo
dataSetNeo1$Count <- log2(dataSetNeo$Count+1)
dataMean <- dataSetNeo1 %>% group_by(DIAGNOSIS.Alias, Alteration) %>%
            mutate(CountMean= sum(Count)) %>% 
            dplyr::select(DIAGNOSIS.Alias,Alteration,CountMean) %>% distinct()
dataMean$Alteration <- factor(dataMean$Alteration, ordered = TRUE, levels =c("FusionNeoAntigenCount","indel.Neoantigen","snv.Neoantigen"))
View(dataMean)

ggplot(data=dataMean,aes(x=DIAGNOSIS.Alias,y=CountMean, fill=Alteration ))+
  geom_bar(stat="identity", width = 1, colour = "black") +
  coord_polar(theta = "x")+
  theme_bw() +
  scale_fill_brewer(palette="Set2")+
  #scale_fill_hue(l=40) +
  xlab("")+ylab("Log2 Total Counts")+ggtitle("Total sum of tumor specific neoantigen counts")+
  theme(legend.position="bottom",
        text = element_text(size=15 ),
        plot.title = element_text(hjust = 0.5))


## Testing
FinalMutationEdit1 <- FinalMutation %>% dplyr::mutate(Alteration = gsub("frameshift deletion|frameshift insertion|nonframeshift insertion|nonframeshift deletion","indel",Exonic.function)) %>% 
  dplyr::mutate(func = gsub("nonsynonymous SNV|stopgain|stoploss","SNV",Alteration)) %>% group_by(DIAGNOSIS.Alias) %>% 
  dplyr::summarize(CountMean = mean(n()))
View(FinalMutationEdit1)

test <- FinalMutationEdit1 %>% filter(DIAGNOSIS.Alias == "OS")
test %>% group_by(Patient.ID, func) %>% dplyr::summarise(AltSumPerPat = n()) %>% group_by(func) %>% dplyr::summarise(AltSumPerPatMean= mean(AltSumPerPat))


