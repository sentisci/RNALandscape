# Transcriptome gene expression analysis (RNASeq Analysis). 

## Highlights
1. Start to end analysis of whole transciptome datasets
2. Easy use of multiple differential gene expression analysis tools like EdgeR, DESeq2 and LimmaVoom at once.
3. Perform multiple comparisions at once
   for example:
   To perform Differntial gene expression analysis for following groups
   
   ### Control groups ##
    Brain          <- c("NS.cerebellum","NS.cerebrum")  
    germline       <- c("NS.testis","NS.ovary")  
    vitalNormals   <- c("NS.heart","NS.kidney","NS.liver","NS.lung") 
    muscle         <- c("skeletal muscles")
     
   ### Condition group
   Tumors         <-  c("ASPS","DSRCT", "EWS" ,"HBL", "ML", "NB" ,"OS", "RMS", "SS", "Teratoma" ,"UDS" ,"YST","WT", "CCSK")  
   NB             <-  c("NB.MYCN.NA","NB.MYCN.A", "NB.Unknown")  
   RMS            <-  c("RMS.FP" , "RMS.FN" )  
         
   ***
   ## Groups that can be compared 
   
   ### N:N comparisons
   
   Normal = list(list("Brain"=Brain,each=FALSE), list("vitalNormals"=vitalNormals,each=FALSE)  )  
   Tumors = list(list("Tumors"=Tumors,each=FALSE) )   
  
   #### With `each`= FALSE ( for Tumors )  
   | Group 1 | group 2 |
   | :---         |     :---:      |
   | Brain     | Tumor    |
   #### With `each`= TRUE  ( for Tumors )   
   | Group 1 | group 2 |
   | :---         |     :---:      |
   | Brain     | ASPS    |
   | Brain     | OS       |
   | Brain     | CCSK       |

   vitalNormals = list(list("vitalNormals"=vitalNormals,each=FALSE)  )  
   Tumors = list(list("Tumors"=Tumors,each=FALSE) )   
   
   #### With `each`= FALSE ( for Both vitalNormals & Tumors )  
   | Group 1 | group 2 |
   | :---         |     :---:      |
   | vitalNormals | Tumors |
   #### With `each`= TRUE  ( for Both vitalNormals & Tumors )  
   | Group 1 | group 2 |
   | :---         |     :---:      |
   | NS.heart     | ASPS    |
   | NS.heart     | OS       |
   | NS.heart     | CCSK       |
   | NS.lung     | ASPS    |
   | NS.lung     | OS       |
   | NS.lung     | CCSK       |
      
   ***
   
   ### 1:1 comparisons  
   Normal = list(list("Brain"=Brain,each=FALSE), list("muscle"=muscle,each=FALSE)  
   Tumors = list(list("NB"=NB,each=FALSE), list("RMS"=RMS,each=FALSE)  
   #### With `each`= FALSE & `OneToOne` = TRUE  ( for Both Normal & Tumors )  
   | Group 1 | group 2 |
   | :---         |     :---:      |
   | Normal     | Tumors    |
   | Brain     | NB       |
   | muscle     | RMS       |
   
