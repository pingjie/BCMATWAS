Requirement:
1. R (V4.0+)
2. Linux/Unix for running Plink and S-PrediXcan

R packages: data.table, glmnet, sqldf

Steps:
1. Run JTI.R / PrediXcan.R / UTMOST.R to build different models. Example: Rscript JTI.R ENSG00000007237.18 Komen EUR
2. Merge Models by running: ./mergeModels.sh
3. Run S-PrediXcan
4. GCTA: ./RunGCTA.sh
5. COJO: ./COJO_assoc_analysis.sh
6. 
