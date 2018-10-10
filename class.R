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
                                          self$plotsDir, self$plotsDataDir, self$diffGeneExpAnaDir,
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
    diffGeneExpAnaDir       = NULL,
    metaDataDF              = NULL,
    annotationDF            = NULL,
    initialize              = function(date = NA, time =NA, projectName = NA, annotationRDS = NA, outputPrefix = NA,
                                       filterGenes = NA, filterGeneMethod = NA, factorName = NA, metaDataFileName = NA, 
                                       workDir = NA, outputdirRDSDir = NA, outputdirTXTDir = NA,
                                       gseaDir = NA, plotsDir = NA, plotsDataDir = NA, 
                                       diffGeneExpAnaDir = NA) {
      
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
      self$diffGeneExpAnaDir <- diffGeneExpAnaDir
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
    readTXTFiles  = function(x, fileSuffix=NA, colNameSelect=NA, primaryID=NA )
    {
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
    mergeTXTFiles = function( x, fileSuffix=NA, colNameSelect=NA, primaryID=NA )
    {
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
    getMergedMatrix = function(dir = NA, fileFormat = NA, colNameSelect = NA, colIndexSelect = NA, isRowNames = FALSE, rowNamesColInFile = NA,
                               fileSuffix=NA, primaryID=NA)
    {
      
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
                 mergedDataList        <- lapply( lapply(dirs, list.files, full.names=T)[2:3],
                                                  private$mergeTXTFiles, fileSuffix = fileSuffix, colNameSelect = colNameSelect, primaryID=primaryID )
                 mergedData            <-  purrr::reduce( mergedDataList,  full_join, by=primaryID )
                 mergedData            <-  tibble::column_to_rownames(mergedData, var=primaryID)
                 return(mergedData)
              }
       )
    },
    consolidateDF = function(df = NA, featureName=NA, funcName=NA, colsToExclude = NA)
    {
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
    matchAndChangeColNames = function()
    {
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
    fpkmToTpm = function(fpkm){
      exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
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
      
      if( private$featureType == "Exon")  { private$featureType = "ExonID"  }
      if( private$featureType == "Transcript"){ private$featureType = "TranscriptID" }
      if( private$featureType == "Gene"){ private$featureType = "GeneID" }
      
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
      
      switch (private$packageRNAseq,
              
              "edgeR"= {
                
                ##  Make EdgeR Object
                private$GeneDF                  <- DGEList(counts=private$countObj, genes=private$miniAnnotationDF, group = private$design)
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
  
  DGEobj        = NULL,
  expressionDF  = NULL,
  metadataDF    = NULL,
  packageRNAseq = NULL,
  groupColumnName      = NULL,
  samplesColumnName    = NULL,
  modelDesign          = NULL,
  GeneDF_Dispersion    = NULL,
  GeneDF_DiffExp       = NULL,
  expressionUnit       = NULL,
  modelGroup           = NULL,  
  group1Name           = NULL, 
  group1Count          = NULL,
  group2Name           = NULL,
  group2Count          = NULL,
  subSetSamples        = NULL,
  setUP         = function(x) {
    
    ## Get the group names
    private$group1Name = names(self$group1)
    private$group2Name = names(self$group2)
    
    ## Make Custom Filters
    group1.group2.Filt   <- interp(~y %in% x, .values=list(y = as.name(private$groupColumnName), x = c(private$group1Name, private$group2Name)))
    factorizeGroupColumn <- interp( ~factor(groupColumnName, levels = val, ordered = TRUE), 
                                    groupColumnName=as.name(private$groupColumnName), val=c(private$group1Name, private$group2Name) )
    
    ## Renaming the groups by groupName in annotationObj
    private$metadataDF[,private$groupColumnName] <- private$changeGroupName(self$group1, private$group1Name, private$metadataDF[,private$groupColumnName])
    private$metadataDF[,private$groupColumnName] <- private$changeGroupName(self$group2, private$group2Name, private$metadataDF[,private$groupColumnName])
    
    ## Renaming the groups by groupName in expressionObj
    private$DGEobj$samples$group <- private$changeGroupName(self$group1, private$group1Name, private$DGEobj$samples$group)
    private$DGEobj$samples$group <- private$changeGroupName(self$group2, private$group2Name, private$DGEobj$samples$group)
    
    
    ## Filter the metadata table for group1 and group2
    private$subSetSamples <- private$metadataDF %>% 
      filter_(group1.group2.Filt) %>% 
      mutate_(.dots = setNames( list(factorizeGroupColumn) , private$groupColumnName )) %>% 
      dplyr::arrange_(.dots = private$groupColumnName)
    
    ## Getting count of each group
    groupCount <- table(private$subSetSamples[,private$groupColumnName])
    private$group1Count <- groupCount[[private$group1Name]]
    private$group2Count <- groupCount[[private$group2Name]]
    
    ## Generate model & design for differential gene expression (edgeR only for now)
    private$modelGroup = c(rep(private$group1Name, private$group1Count ), rep( private$group2Name, private$group2Count ))
    private$modelDesign                    <- model.matrix( ~private$modelGroup )
    
  },
  changeGroupName = function(vectorOfNames, toChangeName, inVectorName){
    return(gsub( paste0(paste0("^",vectorOfNames) %>% paste0(.,"$"), collapse = "|"),toChangeName,inVectorName))
  },
  flattenGroup = function(x){
    unlist( 
      sapply(group1, function(x){  
        if(x$each){
          return(sapply(unlist(unname(x[1])), function(x){ return(list(x)) }))
        } else {
          return(x[1])
        }
      }),
      recursive = FALSE)
  }
  ),
  public    = list
  (
  
  group1         = NULL,
  group2         = NULL,
  filterGenes    = NULL,
  factorsExclude = NULL,
  initialize     = function(DGEobj = NA, expressionDF = NA, metadataDF = NA, packageRNAseq = NA,
                            group1   = NA, group2   = NA, groupColumnName = NA, samplesColumnName = NA, 
                            factorsExclude = NA ){
    
    private$DGEobj             <- DGEobj
    private$expressionDF       <- expressionDF
    private$metadataDF         <- metadataDF
    private$packageRNAseq      <- packageName
    self$group1                <- group1
    self$group2                <- group2
    private$groupColumnName    <- groupColumnName
    private$samplesColumnName  <- samplesColumnName
    self$factorsExclude        <- factorsExclude
    private$setUP()
  },
  edgeRMethod   = function(){
    
    ## Use the appropriate method depending upon the size of the replicates
    if( private$group1Count > 1 & private$group2Count > 1 )   { 
      
      private$GeneDF_Dispersion  <- estimateDisp(private$DGEobj, design = private$modelDesign )
      fit                        <- glmQLFit(private$GeneDF_Dispersion, design = private$modelDesign )
      private$GeneDF_DiffExp     <- glmQLFTest(fit, coef=2)$table
    }
    else  {
      print("Estimating dispersion")
      private$GeneDF_Dispersion  <- estimateGLMCommonDisp(private$DGEobj, method="deviance", robust="TRUE",subset=NULL )
      print("Predicting Differential Gene Expression")
      private$GeneDF_DiffExp     <- exactTest(private$GeneDF_Dispersion, pair = c(unique(private$modelGroup)))$table
    }
    private$GeneDF_DiffExp["FDR"]      <- p.adjust(private$GeneDF_DiffExp$PValue, method="BH")
    private$GeneDF_DiffExp <- tibble::add_column(private$GeneDF_DiffExp, GeneID = rownames(private$GeneDF_DiffExp), .before = 1)
    self$addMeanGroupExpression()
    return(private$GeneDF_DiffExp)
  },
  addMeanGroupExpression = function(){
    
    ## GeneExpression Mean
    expressionDFFilt    <- private$expressionDF %>% dplyr::select_( .dots =  as.character(private$subSetSamples[,private$samplesColumnName]) )
    expressionDFFilt    <- log2(expressionDFFilt+1)
    private$GeneDF_DiffExp[paste0("mean_",private$group1Name)]   <- apply(expressionDFFilt[,c(1:(private$group1Count)), drop=FALSE], 1, mean)
    private$GeneDF_DiffExp[paste0("mean_",private$group2Name)]   <- apply(expressionDFFilt[,(1+private$group1Count):
                                                                                             (private$group1Count+private$group2Count), drop=FALSE], 1, mean)
  }
  )
)

