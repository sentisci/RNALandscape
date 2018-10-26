### Transcriptome gene expression analysis (RNASeq Analysis). 

## Highlights
1. Start to end analysis of whole transciptome datasets
2. Easy use of multiple differential gene expression analysis tools like EdgeR, DESeq2 and LimmaVoom at once.
3. Perform multiple comparisions at once
   for example:
   To perform Differntial gene expression analysis for following groups
   
   ## Control groups ##
    Brain          <- c("NS.cerebellum","NS.cerebrum")
    germline       <- c("NS.testis","NS.ovary")
    vitalNormals   <- c("NS.heart","NS.kidney","NS.liver","NS.lung")


   ## Condition group
   Tumors         <-  c("ASPS","DSRCT", "EWS" ,"HBL", "ML", "NB" ,"OS", "RMS", "SS", "Teratoma" ,"UDS" ,"YST","WT", "CCSK")
   NB             <-  c("NB.MYCN.NA","NB.MYCN.A", "NB.Unknown")
   RMS            <-  c("RMS.FP" , "RMS.FN")
   
   
