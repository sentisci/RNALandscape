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
mergeDataSet <- read.table("../RNASeq.RSEM/ssScore.neoantigen.IGH_C.TRB_C.IGH_E.TRB_E.Neo.Fusion.Final_Mut.preFinal_Mut.Meta.v6.txt", sep="\t",
                           header = TRUE)
mergeDataSet.T <- mergeDataSet %>% dplyr::filter(LIBRARY_TYPE == "Tumor")
mergeDataSet.Neoantigens <- mergeDataSet.T %>% dplyr::select(Sample.ID, DIAGNOSIS.Alias, VariantNeoAntigenCount, 
                                                             FusionNeoAntigenCount )
dataSetNeo <- tidyr::gather(mergeDataSet.Neoantigens, "Diagnosis", "Count", VariantNeoAntigenCount, 
                            FusionNeoAntigenCount )
View(dataSetNeo)


ggplot(data=dataSetNeo,aes(x=DIAGNOSIS.Alias,y=log2(Count+1), fill=Diagnosis ))+
  geom_bar(stat="identity") +
  coord_polar()+
  scale_fill_brewer(palette="Set1")+
  scale_fill_hue(l=40) + 
  xlab("")+ylab("")

dataSetNeo1 <- dataSetNeo
dataSetNeo1$Count <- log2(dataSetNeo$Count+1)
dataMean <- dataSetNeo1 %>% group_by(DIAGNOSIS.Alias) %>% mutate(CountMean= mean(Count)) %>% 
                    dplyr::select(DIAGNOSIS.Alias,Diagnosis,CountMean) %>% distinct()

ggplot(data=dataMean,aes(x=DIAGNOSIS.Alias,y=CountMean, fill=Diagnosis ))+
  geom_bar(stat="identity") +
  coord_polar()+
  scale_fill_brewer(palette="Set2")+
  #scale_fill_hue(l=40) +
  xlab("")+ylab("log2 Counts")+ggtitle("Average tumor specific neoantigen counts.")+
  theme(legend.position="bottom") +
  theme(panel.grid.major = element_line(colour = "Blue")) + 
  theme_bw()




