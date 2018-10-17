library(R6)
library(assertthat)
library(microbenchmark)

# ### Project SetUp Class ####
ProjectSetUp <- R6Class(
  classname = "ProjectSetUP",
  portable = TRUE,
  private  = list(
    
    checkDirExists = function() {
      
      projectDirPath = paste(self$workDir, self$projectName, sep="/")
      ## Check & create if the project directory exists ?? 
      if(!dir.exists(file.path(self$workDir))) {
        warning(paste(self$workDir, "No Such directory exsist !!"))
        paste0("Creating working dir for this project ", projectDirPath)
        dir.create(self$workDir)
      }
      
      ## Create all the required dirs for the project.
      dir.create(projectDirPath)
      fileDirs <- paste(projectDirPath, c(self$outputdirRDSDir, self$outputdirTXTDir, self$gseaDir, 
                                          self$plotsDir, self$plotsDataDir, self$DiffGeneExpAnaDir,
                                          "GeneCountsInput", "TranscriptCountsInput", "ExonCountsInput" ), sep="/")
      lapply(fileDirs, function(x){ if(!dir.exists(x)) dir.create(x) })
    },
    ## Generate filters to exclude given list
    generateFiltersToExclude = function(factorsToExclude=NA){
      
      lapply(factorsToExclude, function(x){
        subList=x
        sublistNames=names(x)
        paste0("!(",
               paste0(sapply(sublistNames, function(y){
                 
                 #print(subList[y])
                 filter <- paste0(y," %in% ","\"",unlist(unname(subList[y])), "\"" )
                 return(filter)
               }), collapse = " & " ),
               ")"
        )
      })
    },
    readMetaData = function() {
      self$metaDataDF <- read.csv(paste0(self$workDir,"/",self$projectName, "/",self$metaDataFileName), sep="\t", header = T);
      print(paste0("Dimension of metadata is ", paste(dim(self$metaDataDF))))
      if(!is.na(self$factorsToExclude)){
        filters <- private$generateFiltersToExclude(factorsToExclude=self$factorsToExclude)
        test <- sapply(filters, function(x){
          self$metaDataDF <- dplyr::filter_(self$metaDataDF, .dots=list(x))
        })
      }
      print(paste0("Dimension of metadata after applying parameter 'factorsToExclude' is ", paste(dim(self$metaDataDF)[1])))
    },
    readAnnotation = function() {
      self$annotationDF <- readRDS(self$annotationRDS) %>% as.data.frame()
    },
    readProteinCoding = function(){
      self$pcDF <- as.data.frame( readRDS(self$pcRDS) )
      colnames(self$pcDF)[c(1,2)] <- c("GeneName","GeneID")
    }
    ,
    readTranscriptionFactor = function(){
      self$tfDF <- readRDS(self$tfRDS)
    },
    readCellSurface = function(){
      self$csDF <- readRDS(self$csRDS)  %>% dplyr::filter(NewCount >= 5)
    },
    readCancerGermlineAntigens = function(){
      self$cgaDF <- readRDS(self$cgaRDS) 
    }
  ),
  public = list(
    date                    = NULL,
    time                    = NULL,
    projectName             = NULL,
    annotationRDS           = NULL,
    outputPrefix            = NULL,
    filterGenes             = NULL,
    filterGeneMethod        = NULL,
    factorName              = NULL,
    metaDataFileName        = NULL,
    workDir                 = NULL,
    outputdirRDSDir         = NULL,
    outputdirTXTDir         = NULL,
    gseaDir                 = NULL,
    plotsDir                = NULL,
    plotsDataDir            = NULL,
    DiffGeneExpAnaDir       = NULL,
    DiffGeneExpRDS          = NULL,
    metaDataDF              = NULL,
    annotationDF            = NULL,
    pcDF                    = NULL,
    pcRDS                   = NULL,
    tfDF                    = NULL,
    tfRDS                   = NULL,
    csDF                    = NULL, 
    csRDS                   = NULL,
    cgaDF                   = NULL,
    cgaRDS                  = NULL,
    factorsToExclude        = NULL,
    initialize              = function(date = NA, time =NA, projectName = NA, annotationRDS = NA, outputPrefix = NA,
                                       filterGenes = NA, filterGeneMethod = NA, factorName = NA, metaDataFileName = NA, 
                                       workDir = NA, outputdirRDSDir = NA, outputdirTXTDir = NA,gseaDir = NA, plotsDir = NA, 
                                       plotsDataDir = NA, DiffGeneExpAnaDir = NA, DiffGeneExpRDS = NA, pcRDS = NA,tfRDS=NA,
                                       csRDS =NA, cgaRDS=NA,factorsToExclude=NA) {
      
      self$date <- date
      self$time <- time
      self$projectName <- projectName
      self$annotationRDS <- annotationRDS
      self$outputPrefix <- outputPrefix
      self$filterGenes <- filterGenes
      self$filterGeneMethod <- filterGeneMethod
      self$factorName <- factorName
      self$metaDataFileName <- metaDataFileName
      self$workDir <- "C:/Users/sindiris/R Scribble/"
      self$outputdirRDSDir <- outputdirRDSDir
      self$outputdirTXTDir <- outputdirTXTDir
      self$gseaDir <- gseaDir
      self$plotsDir <- plotsDir
      self$plotsDataDir <- plotsDataDir
      self$DiffGeneExpAnaDir <- DiffGeneExpAnaDir
      self$DiffGeneExpRDS <- DiffGeneExpRDS
      self$pcRDS <- pcRDS
      self$tfRDS <- tfRDS
      self$csRDS <- csRDS
      self$cgaRDS <- cgaRDS
      self$factorsToExclude <- factorsToExclude
      private$checkDirExists()
      private$readMetaData()
      private$readAnnotation()
      if (!is.na(pcRDS)){ private$readProteinCoding() }
      if (!is.na(tfRDS)){ private$readTranscriptionFactor() }
      if (!is.na(csRDS)){ private$readCellSurface() }
      if (!is.na(csRDS)){ private$readCancerGermlineAntigens() }
      print(paste0("Changing to working directory ", self$workDir))
    }
  ) 
)

# ### Utilities Class ####

CoreUtilities <- R6Class(
  classname = "CoreUtilities",
  portable  = TRUE,
  ## Read CSV or TXT files
  private   = list(
    readTXTFiles  = function(x, fileSuffix=NA, colNameSelect=NA, primaryID=NA ){
      if(!is.na(fileSuffix)) 
      {
        sampleName <- gsub(fileSuffix, "", basename(x)) 
      } else { 
        sampleName <- basename(x) 
      }
      ## print(paste0("I am in readTXTFiles ", sampleName, "  ", colNameSelect))
      rdsObj <- fread(x, sep="\t", header = TRUE) %>% rename_(.dots=setNames(colNameSelect, sampleName))
      return(rdsObj[, c(primaryID,sampleName), with=FALSE])
      
    },
    ## Merge CSV or TXT files
    mergeTXTFiles = function( x, fileSuffix=NA, colNameSelect=NA, primaryID=NA ){
      selectedFileList    <- x[which(basename(x) %in% self$allFileList)]
      print(paste0("Selecting ", length(selectedFileList), " files out of ", length(x), "from the given folder"))
      dataMatrixLists     <- lapply(selectedFileList, private$readTXTFiles, fileSuffix=fileSuffix,colNameSelect=colNameSelect, primaryID=primaryID)
      dataMatrix          <- purrr::reduce(dataMatrixLists, full_join, by=primaryID)
      return(dataMatrix)
    }
  ),
  public    = list(
    workDir                = NULL,
    projectName            = NULL,
    allFileList            = NULL,
    initialize             = function(ProjectSetUpObject = NA ){
      
      assert_that("ProjectSetUP" %in% class(ProjectSetUpObject), 
                  msg="Please setup Project using ProjectSetUp class !!\nProjectSetUpObject cannot be NA !!")
      self$workDir     <- ProjectSetUpObject$workDir
      self$projectName <- ProjectSetUpObject$projectName
    },
    ## Get merged matrix
    getMergedMatrix = function(dir = NA, fileFormat = NA, colNameSelect = NA, colIndexSelect = NA, isRowNames = FALSE, rowNamesColInFile = NA,
                               fileSuffix=NA, primaryID=NA, metadata=NA, metadataFileRefCol=NA){
      
      if( !is.na(colNameSelect) & !is.na(colIndexSelect)) stopifnot("Use either colNameSelect or colIndexSelect but not both")
      if( is.na(colNameSelect) & is.na(colIndexSelect)) stopifnot("Use either colNameSelect or colIndexSelect but not both")
      assert_that( !is.na(fileFormat)  , msg ="Please provide file format to parse. Following are supported, \"rds\", \"csv or text\", \"feather\" ")
      #
      if( is.na(colNameSelect) & is.na(colIndexSelect)) {
        print("First column of every file will be selected.Or use either \"colNameSelect\" or \"colIndexSelect\"")
        colIndexSelect = 1
      }
      if ( is.na(colNameSelect) ) { 
        colNameSelect  =  colIndexSelect
      } 
      if(isRowNames) assert_that(!is.na(rowNamesColInFile), msg= "Please provide row names index , \"rowNamesColInFile\" can't be NA ")
      
      self$allFileList = paste0(as.character(metadata[,metadataFileRefCol]),fileSuffix)
      
      dirs                <- list.dirs(paste0(self$workDir,"/",self$projectName,"/",dir))[-1]
      folderNames         <- basename(dirs)
      
      switch(fileFormat,
             
             #        "rds"     = {
             #          mergedDataList        <- lapply( lapply(dirs, list.files, full.names=T),
             #                                           private$mergeRDSFiles, fileFormat = fileFormat, colInterest = colIndexSelect, fileSuffix = fileSuffix,
             #                                           colNameSelect = colNameSelect )
             #          rowNames              <- rownames(mergedDataList[[1]])
             #          mergedData            <- dplyr::bind_cols( mergedDataList)
             #          rownames(mergedData)  <- rowNames
             #          return(mergedData)
             #        },
             "txt"    = {
               mergedDataList        <- lapply( lapply(dirs, list.files, full.names=T),
                                                private$mergeTXTFiles, fileSuffix = fileSuffix, colNameSelect = colNameSelect, primaryID=primaryID )
               mergedData            <-  purrr::reduce( mergedDataList,  full_join, by=primaryID )
               mergedData            <-  tibble::column_to_rownames(mergedData, var=primaryID)
               return(mergedData)
             }
      )
    },
    ## Consolidate the matrix
    consolidateDF = function(df = NA, featureName=NA, funcName=NA, colsToExclude = NA){
      assert_that( funcName %in% c("sum", "max", "mean"), msg= "For now duplicate featureName can be consolidated only by sum, max or mean" )
      ### Get character or factor columns
      listColsToExclude = names(which(sapply(df, class) == "character" | sapply(df, class) == "factor"))
      assert_that(length(listColsToExclude) != 0, msg="Character and Factor columns found. Please make sure the object has all numeric columns !!")
      if(!is.na(colsToExclude)) assert_that(colsToExclude %in% colnames(df), msg = "colsToExclude not found in the data frame !!")
      
      assert_that(!is.na(featureName), msg="featureName cannot be NA")
      
      if(funcName == "sum")  sumFunc <- function(x)  sum(x)
      if(funcName == "max")  sumFunc <- function(x)  max(x)
      if(funcName == "mean") sumFunc <- function(x)  mean(x)
      
      if(is.na(colsToExclude)) df[, lapply(.SD, sumFunc), by=featureName]
      else df[, lapply(.SD, sumFunc), by=featureName, .SDcols = -colsToExclude]
    },
    ## matchAndChangeColNames columns
    matchAndChangeColNames = function(){
      metaDT <-  setDT(rnaseqProject$metaDataDF)
      setDT(metaDT)
    },
    ## Factorise columns
    factorizeColumn = function( toFactor, asFactor){
      factorColumn <- factor(toFactor, levels=asFactor, ordered = TRUE )
      return(factorColumn)
    },
    ## Keep Protein Coding
    keepProteinCoding = function(geneMatrix=NULL, annotation=NULL, featureType="GeneID") {
      ###Keep only Protein Coding Genes
      PC <- read.table("C:/Users/sindiris/R Scribble/Annotation RDS/HGNC-protein-coding-List.txt", header = T, sep="\t")
      annotationPC <- annotation %>% filter(GeneName %in% PC$Genes)
      foundGenes <- which(row.names(geneMatrix) %in% annotationPC[,featureType] )
      geneMatrixNew <- geneMatrix %>% data.frame() %>% .[foundGenes,]
      annotationPC <- annotationPC %>% mutate(GeneID = self$factorizeColumn(annotationPC$GeneID, row.names(geneMatrixNew)))
      return(list(geneMatrixNew, annotationPC) )
    },
    ## convert FPKM to TPM
    fpkmToTpm = function(fpkm){
      exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
    },
    ## subsetMetaData
    subsetMetaData = function(colnamesDF=NA){
      df <- dplyr::left_join(colnamesDF,rnaseqProject$metaDataDF, by="SAMPLE_ID") %>% filter(complete.cases(.))
      if(!is.null(rnaseqProject$metaDataDF)) {
        rnaseqProject$metaDataDF <- df
      }
    },
    ## Annotate a gene expression df with "GeneID" as primary key
    featureNameAnot = function(querryDF=NA, identifier=NA){
      annotDF <- dplyr::left_join(rnaseqProject$annotationDF, queryDF, by=identifier)
      print(paste("Annotating expression DF"))
      print(dim(annotDF))
      return(annotDF)
    }
  )
)

# ### GeneExpression normalization using multiple tools ####
# "keepProteinCoding(), fpkmToTpm()" function in another class

GeneExpNormalization <- R6Class(
  classname = "GeneExpNormalization",
  inherit   = CoreUtilities, 
  portable  = TRUE,
  private   = list(
    GeneDF               = NULL,
    GeneDFNorm           = NULL,
    miniAnnotationDF     = NULL,
    countObj             = NULL,
    featureType          = NULL,
    packageRNAseq        = NULL,
    annotationDF         = NULL,
    design               = NULL,
    proteinCodingOnly    = NULL,
    setUp                = function(){
      
      if( private$featureType == "Exon")      { private$featureType = "ExonID"       }
      if( private$featureType == "Transcript"){ private$featureType = "TranscriptID" }
      if( private$featureType == "Gene")      { private$featureType = "GeneID"       }
      
      private$annotationDF %<>% mutate(!!private$featureType := super$factorizeColumn(private$annotationDF[[private$featureType]], 
                                                                                      row.names(private$countObj))) %>% 
        dplyr::arrange_(.dots=c(private$featureType))
      
      if(private$proteinCodingOnly) {
        keepProteinCodingOut     <- super$keepProteinCoding(geneMatrix = private$countObj, 
                                                            annotation=private$annotationDF, featureType=private$featureType)
        private$countObj         <- keepProteinCodingOut[[1]]
        private$annotationDF     <- keepProteinCodingOut[[2]]
      }
      
      private$miniAnnotationDF <- private$annotationDF[,c(private$featureType, "Length")]
      colnames(private$miniAnnotationDF) <- c("GeneID", "Length")
      
      switch(private$packageRNAseq,
             
             "edgeR"= {
               
               ##  Make EdgeR Object
               private$GeneDF                  <- DGEList(counts=private$countObj, genes=private$miniAnnotationDF, 
                                                          group = as.factor(as.character(private$design)) )
               ## Estimate Normalising Factors
               private$GeneDFNorm              <- calcNormFactors(private$GeneDF)
               
             },
             "deseq2"= {
               
               ##  Make deseq2 Object
               condition                       <- private$design
               private$GeneDF                  <- DESeqDataSetFromMatrix(countData = private$countObj, colData = DataFrame(condition), design = ~condition)
               mcols(private$GeneDF)$basepairs <- private$miniAnnotationDF$Length-174
               ## Estimate SizeFactors
               private$GeneDFNorm              <- estimateSizeFactors(private$GeneDF)
             }
      )
      
    }
  ),
  public    = list(
    
    initialize           = function(countObj = NA, featureType = NA, packageRNAseq = NA, annotationDF = NA, design = NA,
                                    proteinCodingOnly = TRUE ) 
    {
      
      assert_that( class(countObj) == "matrix", msg = "Please provide raw count in matrix format")
      assert_that( nrow(countObj) == nrow(annotationDF), msg = "Please provide gene annotation as a \"data.frame\" ")
      assert_that( is.logical(proteinCodingOnly) , msg = "proteinCodingOnly parameter should be boolean")
      
      private$countObj <- countObj
      private$featureType <- featureType
      private$packageRNAseq <- packageRNAseq
      private$annotationDF  <- annotationDF
      private$design        <- design
      private$proteinCodingOnly <- proteinCodingOnly
      private$setUp()
    },
    
    edgeRMethod      = function(x) {
      
      assert_that(private$packageRNAseq == "edgeR", msg = paste0("GeneExpNormalization Object was created for ",private$packageRNAseq,
                                                                 ". Please create a new GeneExpNormalization Object for edgeR"))
      
      assert_that(x %in% c("CPM", "TMM-RPKM", "TPM", "NormFactorDF"), msg = "This function can only generate \"CPM\", \"TMM-RPKM\", \"TPM\" values ")
      
      if(x == "NormFactorDF") return(private$GeneDFNorm)
      if(x == "CPM" )         return( as.data.frame(cpm(private$GeneDFNorm,  normalized.lib.sizes = TRUE,log = FALSE))   )
      if(x == "TMM-RPKM" )   return( as.data.frame(rpkm(private$GeneDFNorm, normalized.lib.sizes = TRUE, log = FALSE))  )
      if(x == "TPM" )         return(apply(rpkm(private$GeneDFNorm, normalized.lib.sizes = TRUE), 2 , super$fpkmToTpm)         )
      
    },
    deseq2           = function(x) {
      
      assert_that(private$packageRNAseq == "deseq2", msg = paste0("GeneExpNormalization Object was created for ",private$packageRNAseq,
                                                                  ". Please create a new GeneExpNormalization Object for deseq2"))
      
      assert_that(x %in% c("FPM", "TMM-RPKM", "TPM", "RLOG", "VST"), 
                  msg = "This function can only generate \"CPM\", \"TMM-RPKM\", \"TPM\", \"RLOG\", \"VST\" values ")
      
      
      if(x == "NormFactorDF") return(private$GeneDFNorm)
      if(x == "FPM" )         return( as.data.frame(fpm(object = private$GeneDFNorm,  robust = TRUE))    )
      if(x == "TMM-RPKM" )   return( as.data.frame(fpkm(object = private$GeneDFNorm,  robust = TRUE))   )
      if(x == "TPM" )         return( apply(as.data.frame(fpkm(object = private$GeneDFNorm,  robust = TRUE)), 2 , fpkmToTpm))
      if(x == "RLOG" )        return( rlog(private$GeneDFNorm, blind = TRUE)              )
      if(x == "VST" )         return( vst(private$GeneDFNorm, blind = TRUE))
    }  
  )
)

# ### Differntial Gene Expression Analysis ####

DifferentialGeneExp <- R6Class(
  classname = "DifferentialGeneExp",
  portable  = TRUE,
  private   = list
  (
  
  countObj            = NULL,
  metadataDF          = NULL,
  group1Flatten       = NULL,
  group2Flatten       = NULL,
  GeneDF_DiffExp      = NULL,
  pairGroup1Name      = NULL,
  group1Samples       = NULL,
  pairGroup2Name      = NULL,
  group2Samples       = NULL,
  setUP               = function(x) {},
  changeGroupName     = function(vectorOfNames, toChangeName, inVectorName){
    return(gsub( paste0(paste0("^",vectorOfNames) %>% paste0(.,"$"), collapse = "|"),toChangeName,inVectorName))
  },
  flattenGroup        = function(x){
    unlist(
      lapply(x, function(y){  
        if(y$each){
          return( as.list(sapply(unlist(unname(y[1])), function(y){ return(y) }) ) )
        } else {
          return(y[1])
        }
      }), recursive = FALSE)
  },
  makePairs           = function(x, y){
    pairs = unlist(lapply(x, function(X) {
      lapply(y, function(Y) {
        c(X, Y)
      })
    }), recursive=FALSE)
    return(pairs)
  },
  DiffGeneExp         = function( group1Count = NA, group2Count = NA, pairGroup1Name = NA,pairGroup2Name = NA, DGEobj = NA  ){
    
    ## Use the appropriate method depending upon the size of the replicates
    if( group1Count > 1 & group2Count > 1 )   {

      print("Groups have replicates")
      ## Generate model & design for differential gene expression (edgeR only for now)
      
      #Using exactTest()
      modelGroup   <- factor( c( rep(pairGroup1Name, group1Count),rep(pairGroup2Name, group2Count)) )
      modelDesign  <- model.matrix( ~modelGroup )
      print(paste(pairGroup1Name, group1Count,pairGroup2Name,group2Count ))
      ## Estimate Dispersion using estimateDisp() 
      GeneDF_Dispersion  <- estimateDisp(DGEobj, design = modelDesign )
      print("Predicting differntially expressed genes using exactTest()")
      private$GeneDF_DiffExp <- exactTest(GeneDF_Dispersion, pair = c(unique(modelGroup)))$table
      
      #print(head(private$GeneDF_DiffExp))
      #Using Quasilikelihood ratio test()
      #modelGroup   <- factor( c(rep(pairGroup2Name, group2Count), rep(pairGroup1Name, group1Count)) )
      #modelDesign  <- model.matrix( ~modelGroup )
      #print("Predicting differntially expressed genes using Quasi Likelihood ratio test glmQLFit() & glmQLFTest()")
      #fit                        <- glmQLFit(GeneDF_Dispersion, design = modelDesign )
      #private$GeneDF_DiffExp     <- glmQLFTest(fit, coef=2)$table
      
      #print("Predicting differntially expressed genes using LimmaVoom")
      # GeneDF_Dispersion_v   <- voom(DGEobj, modelDesign)
      # fit_v                 <- lmFit(GeneDF_Dispersion_v, modelDesign)
      # fit_v                 <- eBayes(fit_v)
      # modelGroup_v   <- factor( c(rep(pairGroup2Name, group2Count), rep(pairGroup1Name, group1Count)) )
      # modelDesign_v  <- model.matrix( ~modelGroup_v )
      # GeneDF_Dispersion_vwts  <- voomWithQualityWeights(DGEobj, design=modelDesign_v, normalization="none", plot=FALSE)
      # fit_vwts                <- eBayes(lmFit(GeneDF_Dispersion_vwts,design=modelDesign_v))
      # GeneDF_DiffExp_V        <- topTable(fit_vwts,coef=2,number=length(fit_vwts$genes[,1]),sort.by="none")
      # print(head(GeneDF_DiffExp_V))
      
      
      
    }
    else  {
      ## Generate model & design for differential gene expression (edgeR only for now)
      modelGroup   <- factor(c(rep(pairGroup1Name, group1Count), rep( pairGroup2Name, group2Count   )))
      modelDesign  <- model.matrix( ~modelGroup )
      print(modelGroup)
      print(modelDesign)
      ## Estimate Dispersion using estimateDisp() 
      GeneDF_Dispersion  <- estimateDisp(DGEobj, design = modelDesign )
      print("Predicting differntially expressed genes using exact test")
      private$GeneDF_DiffExp     <- exactTest(GeneDF_Dispersion, pair = c(unique(modelGroup)))$table
      #fit                        <- glmQLFit(GeneDF_Dispersion, design = modelDesign )
      #private$GeneDF_DiffExp     <- glmLRT(fit, coef=2)$table 
    }
    
    private$GeneDF_DiffExp["FDR"]   <- p.adjust(private$GeneDF_DiffExp$PValue, method="BH")
    private$GeneDF_DiffExp          <- tibble::add_column(private$GeneDF_DiffExp, GeneID = rownames(private$GeneDF_DiffExp), .before = 1)
    
    return(private$GeneDF_DiffExp)
  },
  MeanGroupExpression = function(group1Count = NA, group2Count = NA, geneExpMatrix = NA){
    
    ## GeneExpression Mean
    if(group1Count > 1) {
      private$GeneDF_DiffExp[private$pairGroup1Name]        <- apply(geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE], 1, mean)
      #group1MeanExp                                        <- apply(geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE], 1, mean)
      #print(paste(head(group1MeanExp)))
    } else {
      private$GeneDF_DiffExp[private$pairGroup1Name]        <- geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE]
    }
    
    if(group2Count > 1) {
      private$GeneDF_DiffExp[private$pairGroup2Name]        <- apply(geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE], 1, mean)
      #group2MeanExp                                        <- apply(geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE], 1, mean)
      #print(paste(head(group2MeanExp)))
      
    } else {
      private$GeneDF_DiffExp[private$pairGroup2Name]        <- geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE]
    }
    
    private$GeneDF_DiffExp["FDR"]   <- p.adjust(private$GeneDF_DiffExp$PValue, method="BH")
    private$GeneDF_DiffExp["myLogFC"]   <- log2(private$GeneDF_DiffExp[private$pairGroup2Name]+1) - log2(private$GeneDF_DiffExp[private$pairGroup1Name]+1)
    return(private$GeneDF_DiffExp)
  },
  executeDiffGeneExp  = function( pairedList = NA){
    
    ## Get the group names
    pairGroup1List <- private$group1Flatten[  pairedList[1] ]
    private$pairGroup1Name <- names(pairGroup1List)
    pairGroup1 <- unname(unlist(pairGroup1List))
    
    pairGroup2List <- private$group2Flatten[  pairedList[2] ]
    private$pairGroup2Name <- names(pairGroup2List)
    pairGroup2 <- unname(unlist(pairGroup2List))
    
    ## Make Custom Filters
    group1.group2.Filt   <- interp(~y %in% x, .values=list(y = as.name(self$groupColumnName), x = c(pairGroup1, pairGroup2) ))
    
    factorizeGroupColumn <- interp( ~factor(groupColumnName, levels = val, ordered = TRUE), 
                                    groupColumnName=as.name(self$groupColumnName), val=c(pairGroup1, pairGroup2) )
    
    ## Filter the metadata table for group1 and group2
    subSetSamples <- private$metadataDF %>% 
      filter_(group1.group2.Filt) %>% 
      mutate_(.dots = setNames( list(factorizeGroupColumn) , self$groupColumnName )) %>%
      dplyr::arrange_(.dots = self$groupColumnName)
    
    private$group1Samples <- subSetSamples %>% filter_(interp(~y %in% x, .values=list(y = as.name(self$groupColumnName), x = pairGroup1))) %>%
      select_(.dots=list(self$samplesColumnName))
    private$group2Samples <- subSetSamples %>% filter_(interp(~y %in% x, .values=list(y = as.name(self$groupColumnName), x = pairGroup2))) %>%
      select_(.dots=list(self$samplesColumnName))
    
    ## Subset Count matrix
    pairCountMatrix <- private$countObj %>% data.frame() %>% dplyr::select_(.dots=as.character(subSetSamples[,self$samplesColumnName]))
    
    ## Getting count of each group
    groupCount  <- table(subSetSamples[,self$groupColumnName])
    group1Count <- sum(groupCount[pairGroup1])
    group2Count <- sum(groupCount[pairGroup2])
    
    ## Generate model & design for differential gene expression (edgeR only for now)
    modelGroup   <- factor(c(rep(private$pairGroup1Name, group1Count), rep( private$pairGroup2Name, group2Count   )))
    modelDesign  <- model.matrix( ~modelGroup )
    
    ## Generate Gene expression in the desired Unit.
    expressionObj <- GeneExpNormalization$new(
      countObj          = as.matrix(pairCountMatrix), 
      featureType       = self$featureType, 
      packageRNAseq     = self$packageRNAseq, 
      annotationDF      = rnaseqProject$annotationDF, 
      design            = modelDesign[,2], 
      proteinCodingOnly = FALSE
    )
    geneExpression = expressionObj$edgeRMethod(self$expressionUnit)
    
    ## perform differential gene expression.
    private$GeneDF_DiffExp <- private$DiffGeneExp(group1Count = group1Count, group2Count = group2Count, DGEobj = expressionObj$edgeRMethod("NormFactorDF"), 
                                                  pairGroup1Name = private$pairGroup1Name, pairGroup2Name = private$pairGroup2Name)
    
    ## calculat mean expression
    private$GeneDF_DiffExp                <- private$MeanGroupExpression(group1Count = group1Count, group2Count = group2Count, geneExpMatrix = geneExpression)
    private$GeneDF_DiffExp["AvglogFPKM"]  <- apply(geneExpression , 1, mean)
    
    ## Annotate Genes 
    private$GeneDF_DiffExp                <- private$featureNameAnot(querryDF=private$GeneDF_DiffExp, identifier="GeneID")
    
    ## Filter Genes
    private$filterGenes(filterName="all")
    private$filterGenes(filterName="proteinCoding")
    private$filterGenes(filterName="cellsurface")
    private$filterGenes(filterName="transcriptionFactor")
    private$filterGenes(filterName="cancergermlineantigen")
    
    return(private$GeneDF_DiffExp)
  },
  ## Annotate a gene expression df with "GeneID" as primary key
  featureNameAnot = function(querryDF=NA, identifier=NA){
    annotDF <- dplyr::left_join(rnaseqProject$annotationDF, querryDF, by=identifier)
    print(paste("Annotating expression DF"))
    print(dim(annotDF))
    return(annotDF)
  },
  filterGenes = function(filterName="NA"){
    switch(filterName,
           
           "proteinCoding"= {
             
             GeneDF_DiffExp_PC <- na.omit(dplyr::left_join(rnaseqProject$pcDF, private$GeneDF_DiffExp, by="GeneID"))
             #print("Filter matched ", dim(GeneDF_DiffExp_PC)[1], " out of  ", dim(rnaseqProject$pcDF)[1], " given protein coding genes")
             print("filtering for pritein coding genes")
             #pcDF <- private$GeneDF_DiffExp %>% filter(GeneName %in% rnaseqProject$pcDF)
             saveRDS(GeneDF_DiffExp_PC, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                               "/", rnaseqProject$DiffGeneExpRDS,"/",private$pairGroup1Name,"_",
                                               private$pairGroup2Name,"_",filterName,"_",".rds"))
             write.table(GeneDF_DiffExp_PC, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                   "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                   private$pairGroup2Name,"_",filterName,".txt"), sep="\t", row.names = FALSE,
                         quote=FALSE)
           },
           "cellsurface"= {
             
             GeneDF_DiffExp_csDF <- private$GeneDF_DiffExp %>% filter(GeneName %in% as.character(rnaseqProject$csDF[,"GeneName"]))
             GeneDF_DiffExp_csDF <- dplyr::left_join(GeneDF_DiffExp_csDF, rnaseqProject$pcDF, by = "GeneID")
             #print("Filter matched ", dim(GeneDF_DiffExp_csDF)[1], " out of  ", dim(rnaseqProject$csDF)[1], " given CS genes")
             print("filtering for cellsurface genes")
             saveRDS(GeneDF_DiffExp_csDF, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                 "/", rnaseqProject$DiffGeneExpRDS,"/",private$pairGroup1Name,"_",
                                                 private$pairGroup2Name,"_",filterName,".rds"))
             write.table(GeneDF_DiffExp_csDF, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                     "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                     private$pairGroup2Name,"_",filterName,".txt"), sep="\t", row.names = FALSE,
                         quote=FALSE)
             
           },
           "transcriptionFactor"= {
             
             GeneDF_DiffExp_tfDF <- private$GeneDF_DiffExp %>% filter(GeneName %in% as.character(rnaseqProject$tfDF[,"GeneName"]))
             GeneDF_DiffExp_tfDF <- dplyr::left_join(GeneDF_DiffExp_tfDF, rnaseqProject$pcDF, by = "GeneID")
             #print("Filter matched ", dim(GeneDF_DiffExp_tfDF)[1], " out of  ", dim(rnaseqProject$tfDF)[1], " given TF genes")
             print("filtering for transcriptionFactor genes")
             saveRDS(GeneDF_DiffExp_tfDF, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                 "/", rnaseqProject$DiffGeneExpRDS,"/",private$pairGroup1Name,"_",
                                                 private$pairGroup2Name,"_",filterName,".rds"))
             write.table(GeneDF_DiffExp_tfDF, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                     "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                     private$pairGroup2Name,"_",filterName,".txt"), sep="\t", row.names = FALSE,
                         quote=FALSE)
           },
           "cancergermlineantigen"= {
             
             GeneDF_DiffExp_cgaDF <- private$GeneDF_DiffExp %>% filter(GeneName %in% as.character(rnaseqProject$cgaDF[,"GeneName"]))
             GeneDF_DiffExp_cgaDF <- dplyr::left_join(GeneDF_DiffExp_cgaDF, rnaseqProject$pcDF, by = "GeneID")
             #print("Filter matched ", dim(GeneDF_DiffExp_cgaDF)[1], " out of  ", dim(rnaseqProject$tfDF)[1], " given TF genes")
             print("filtering for cancergermlineantigen genes")
             saveRDS(GeneDF_DiffExp_cgaDF, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                 "/", rnaseqProject$DiffGeneExpRDS,"/",private$pairGroup1Name,"_",
                                                 private$pairGroup2Name,"_",filterName,".rds"))
             write.table(GeneDF_DiffExp_cgaDF, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                     "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                     private$pairGroup2Name,"_",filterName,".txt"), sep="\t", row.names = FALSE,
                         quote=FALSE)
           },
           "all" = {
             saveRDS(private$GeneDF_DiffExp, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                    "/", rnaseqProject$DiffGeneExpRDS,"/",private$pairGroup1Name,"_",
                                                    private$pairGroup2Name,"_",filterName,".rds"))
             write.table(private$GeneDF_DiffExp, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                        "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                        private$pairGroup2Name,"_",filterName,".txt"), sep="\t", row.names = FALSE,
                         quote=FALSE)
           }
    )
  }
  ),
  public    = list
  (
  
  group1            = NULL,
  group2            = NULL,
  packageRNAseq     = NULL,
  factorsExclude    = NULL,
  groupColumnName   = NULL,
  samplesColumnName = NULL,
  pairedList        = NULL,
  diffGeneExpList   = NULL,
  expressionUnit    = NULL,
  featureType       = NULL,
  writeFiles        = NULL,
  initialize        = function(countObj = NA, metadataDF = NA, packageRNAseq = NA, expressionUnit = NA, featureType = NA,
                               group1   = NA, group2   = NA, groupColumnName = NA, samplesColumnName = NA,
                               writeFiles = FALSE){
    
    private$countObj           <- countObj
    private$metadataDF         <- metadataDF
    self$expressionUnit        <- expressionUnit
    self$packageRNAseq         <- packageRNAseq
    self$featureType           <- featureType
    self$group1                <- group1
    self$group2                <- group2
    self$groupColumnName       <- groupColumnName
    self$samplesColumnName     <- samplesColumnName
    self$writeFiles            <- writeFiles
    
    ## Get the group names
    private$group1Flatten = private$flattenGroup(self$group1)
    private$group2Flatten = private$flattenGroup(self$group2)
    
    ## Pair elements from two groups
    self$pairedList    = private$makePairs(names(private$group1Flatten), names(private$group2Flatten))
    
    ## printing
    print("printing Group1")
    print(private$group1Flatten)
    print("printing Group2")
    print(private$group2Flatten)
    
  },
  performDiffGeneExp   = function(){
    self$diffGeneExpList <- lapply(self$pairedList, private$executeDiffGeneExp)
    return(self$diffGeneExpList)
  }
  )
)

