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
    readMetaData = function() {
      self$metaDataDF <- read.csv(paste0(self$workDir,"/",self$projectName, "/",self$metaDataFileName), sep="\t", header = T);
    },
    readAnnotation = function() {
      self$annotationDF <- readRDS(self$annotationRDS) %>% as.data.frame()
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
    metaDataDF              = NULL,
    annotationDF            = NULL,
    initialize              = function(date = NA, time =NA, projectName = NA, annotationRDS = NA, outputPrefix = NA,
                                       filterGenes = NA, filterGeneMethod = NA, factorName = NA, metaDataFileName = NA, 
                                       workDir = NA, outputdirRDSDir = NA, outputdirTXTDir = NA,
                                       gseaDir = NA, plotsDir = NA, plotsDataDir = NA, 
                                       DiffGeneExpAnaDir = NA) {
      
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
      private$checkDirExists()
      private$readMetaData()
      private$readAnnotation()
      
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
      print(paste0("I am in readTXTFiles ", sampleName, "  ", colNameSelect))
      rdsObj <- fread(x, sep="\t", header = TRUE) %>% rename_(.dots=setNames(colNameSelect, sampleName))
      return(rdsObj[, c(primaryID,sampleName), with=FALSE])
      
    },
    ## Merge CSV or TXT files
    mergeTXTFiles = function( x, fileSuffix=NA, colNameSelect=NA, primaryID=NA ){
      dataMatrixLists     <- lapply(x, private$readTXTFiles, fileSuffix=fileSuffix,colNameSelect=colNameSelect, primaryID=primaryID)
      dataMatrix          <- purrr::reduce(dataMatrixLists, full_join, by=primaryID)
      return(dataMatrix)
    }
  ),
  public    = list(
    workDir                = NULL,
    projectName            = NULL,
    initialize             = function(ProjectSetUpObject = NA ){
      
      assert_that("ProjectSetUP" %in% class(ProjectSetUpObject), 
                  msg="Please setup Project using ProjectSetUp class !!\nProjectSetUpObject cannot be NA !!")
      self$workDir     <- ProjectSetUpObject$workDir
      self$projectName <- ProjectSetUpObject$projectName
    },
    format = function(...) {
      c(
        paste0("Person:"),
        paste0("  Name: "),
        paste0("  Age:  ")
      )
    },
    ## Get merged matrix
    getMergedMatrix = function(dir = NA, fileFormat = NA, colNameSelect = NA, colIndexSelect = NA, isRowNames = FALSE, rowNamesColInFile = NA,
                               fileSuffix=NA, primaryID=NA){
      
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
  DiffGeneExp         = function( group1Count = NA, group2Count = NA, DGEobj = NA, modelDesign = NA){
    
    ## Use the appropriate method depending upon the size of the replicates
    if( group1Count > 1 & group2Count > 1 )   {
      GeneDF_Dispersion  <- estimateDisp(DGEobj, design = modelDesign )
      fit                        <- glmQLFit(GeneDF_Dispersion, design = modelDesign )
      private$GeneDF_DiffExp     <- glmQLFTest(fit, coef=2)$table
    }
    else  {
      print("Estimating dispersion")
      GeneDF_Dispersion  <- estimateGLMCommonDisp(DGEobj, method="deviance", robust="TRUE",subset=NULL )
      print("Predicting Differential Gene Expression")
      private$GeneDF_DiffExp     <- exactTest(GeneDF_Dispersion, pair = c(unique(modelGroup)))$table
    }
    
    private$GeneDF_DiffExp["FDR"]   <- p.adjust(private$GeneDF_DiffExp$PValue, method="BH")
    private$GeneDF_DiffExp          <- tibble::add_column(private$GeneDF_DiffExp, GeneID = rownames(private$GeneDF_DiffExp), .before = 1)
    
    return(private$GeneDF_DiffExp)
  },
  MeanGroupExpression = function(group1Count = NA, group2Count = NA, geneExpMatrix = NA){
    
    ## GeneExpression Mean
    if(group1Count > 1) {
      private$GeneDF_DiffExp[private$pairGroup1Name]        <- apply(geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE], 1, mean)
    } else {
      private$GeneDF_DiffExp[private$pairGroup1Name]        <- geneExpMatrix[,as.character(private$group1Samples[,1]), drop=FALSE]
    }
    
    if(group2Count > 1) {
      private$GeneDF_DiffExp[private$pairGroup2Name]        <- apply(geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE], 1, mean)
    } else {
      private$GeneDF_DiffExp[private$pairGroup2Name]        <- geneExpMatrix[,as.character(private$group2Samples[,1]), drop=FALSE]
    }
    private$GeneDF_DiffExp["FDR"]   <- p.adjust(private$GeneDF_DiffExp$PValue, method="BH")
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
    modelGroup   <- c(rep(private$pairGroup1Name, group1Count ), rep( private$pairGroup2Name, group2Count ))
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
                                    modelDesign = modelDesign)
    
    ## calculat mean expression
    private$GeneDF_DiffExp                <- private$MeanGroupExpression(group1Count = group1Count, group2Count = group2Count, geneExpMatrix = geneExpression)
    private$GeneDF_DiffExp["AvglogFPKM"]  <- apply(geneExpression , 1, mean)
    
    ## Annotate Genes 
    private$GeneDF_DiffExp                <- private$featureNameAnot(querryDF=private$GeneDF_DiffExp, identifier="GeneID")
    
    ## Filter Genes
    private$filterGenes(filterName="proteinCoding")
    private$filterGenes(filterName="cellsurface")
    private$filterGenes(filterName="transcriptionFactor")
    
    return(private$GeneDF_DiffExp)
  },
  ## Annotate a gene expression df with "GeneID" as primary key
  featureNameAnot = function(querryDF=NA, identifier=NA){
    annotDF <- dplyr::left_join(rnaseqProject$annotationDF, querryDF, by=identifier)
    print(paste("Annotating expression DF"))
    print(dim(annotDF))
    return(annotDF)
  },
  filterGenes = function(filterName=NA){
    switch(private$packageRNAseq,
           
           "proteinCoding"= {
             
             write.table(private$GeneDF_DiffExp, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                        "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                        private$pairGroup2Name,".txt"), sep="\t", row.names = FALSE, quote=FALSE)
             
             
           },
           "cellsurface"= {
             
             write.table(private$GeneDF_DiffExp, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                        "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                        private$pairGroup2Name,".txt"), sep="\t", row.names = FALSE, quote=FALSE)
             
             
           },
           "transcriptionFactor"= {
             
             write.table(private$GeneDF_DiffExp, paste0(rnaseqProject$workDir,"/",rnaseqProject$projectName,
                                                        "/", rnaseqProject$DiffGeneExpAnaDir,"/",private$pairGroup1Name,"_",
                                                        private$pairGroup2Name,".txt"), sep="\t", row.names = FALSE, quote=FALSE)
             
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

