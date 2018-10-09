## Source all classes and packages

source("./utilityPackages.R")
source("./statisticalPackages.R")
source("./class.R")

## Expression Analysis for Landscape paper

rnaseqProject <- ProjectSetUp$new(
  
  date                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[1]],
  time                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[2]],
  projectName             = "RNASeq.RSEM",
  annotationRDS           = "C:/Users/sindiris/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS",
  outputPrefix            = "landscape",
  filterGenes             = TRUE,
  filterGeneMethod        = "bySum",
  factorName              = "DIAGNOSIS.Substatus.Tumor.Normal.Tissue",
  metaDataFileName        = "MetadataMapper.txt",
  outputdirRDSDir         = "GeneRDSOutput",
  outputdirTXTDir         = "GeneTXTOutput",
  gseaDir                 = "GSEA",
  plotsDir                = "Figures",
  plotsDataDir            = "FigureData",
  diffGeneExpAnaDir       = "DiffExpResults"
  
)

corUtilsFuncs <- CoreUtilities$new(  ProjectSetUpObject = rnaseqProject )

rm(mergeObjectsNoDup)
mergeObjectsNoDup <- corUtilsFuncs$getMergedMatrix(dir               = "TPM_Genes", 
                                                   fileFormat        = "txt", 
                                                   colNameSelect     = "expected_count", 
                                                   isRowNames        = TRUE, 
                                                   rowNamesColInFile = 1,
                                                   fileSuffix        = ".genes.results",
                                                   primaryID         = "gene_id")

setDT(mergeObjectsNoDup, keep.rownames = TRUE)
mergeObjectsConso <- corUtilsFuncs$consolidateDF(mergeObjectsNoDup, funcName = "sum", featureName = "rn")
mergeObjectsConso <- mergeObjectsConso %>% data.frame() %>% tibble::column_to_rownames(var = "rn") %>% as.matrix()
















