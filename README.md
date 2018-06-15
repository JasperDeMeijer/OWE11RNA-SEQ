# OWE11RNA-SEQ

Contributers:
- Heleen Severin
- Tim van de Kerkhof
- Shane Pullens
- Jasper de Meijer

RNA-SEQ analysis using R.

Script is developed in a R notebook and wil perform the following:
- Load and read the data into a dataframe
- Normalize the data based on the strain
- Determines the foldchange using edgeR
- Filters the most significant genes based on p-value and foldchange
- Retrieves the pathways of the most significant genes using KEGG
- Updates the annotation to the most recent annotation known by NCBI 
- Produces a .GMX which can be used for a GSEA
- Exports 4 files with foldchanges and annotation
- Exports 1 file with most significant pathways of WCFS1 (Proof concept)

To run this analysis: run RNAseq_analyse.Rmd in Rstudio
When running this script a filechoose will pop up, pick the file containing the RNA-SEQ count data (in this case: Data/RNA-Seq-counts.txt)
