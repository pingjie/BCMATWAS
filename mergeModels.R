library(data.table)
library(qvalue)

args <- commandArgs(trailingOnly = T)
tissueN <- as.numeric(args[1])

method <- args[1]
race <- args[2]
projID <- args[3]

prefix <- paste0(method, ".", projID, ".", race)

if (method == "UTMOST") {
    setwd(paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/models/", method, "/", race))
} else {
    setwd(paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/models/", method, "/", race, "/", projID))
}

### Extra
extraTab <- fread(paste0(prefix, ".Extras.tmp.txt"), header = F)
colnames(extraTab) <- c("gene", "genename", "pred.perf.R", "pred.perf.R2", "n.snps.in.model", "pred.perf.pval")
qobj <- qvalue(p = extraTab$pred.perf.pval, pi0 = 1)
extraTab$pred.perf.qval <- qobj$qvalues
fwrite(extraTab, file = paste0(prefix, ".Extras.all.txt"), sep = "\t")

### Weights
weightTab <- fread(paste0(prefix, ".Weights.tmp.txt"), header = F)
if (method == "UTMOST") {
    colnames(weightTab) <- c("gene", "rsid", "chr_bp", "counted_allele", "ref_allele", "weight", "R", "R2", "p", "Lambda", "r_rt", "r2_rt", "p_rtiter", "iteration_n")
    weightTab <- weightTab[, c("gene", "rsid", "chr_bp", "ref_allele", "counted_allele", "weight", "R", "R2", "p", "Lambda", "r_rt", "r2_rt", "p_rtiter", "iteration_n")]
} else if (method == "PrediXcan") {
    colnames(weightTab) <- c("gene", "genename", "rsid", "chr_bp", "ref_allele", "counted_allele", "weight", "R", "R2", "P", "Lambda")
} else if (method == "JTI") {
    colnames(weightTab) <- c("gene", "rsid", "chr_bp", "ref_allele", "counted_allele", "weight", "R", "R2", "P", "Lambda")
}

fwrite(weightTab, file = paste0(prefix, ".Weights.all.txt"), sep = "\t")

####

weightTab_R0.1 <- subset(weightTab, R2 > 0.01)
extraTab_R0.1 <- subset(extraTab, gene %in% weightTab_R0.1$gene)

fwrite(extraTab_R0.1, file = paste0(prefix, ".Extras.txt"), sep = "\t")
fwrite(weightTab_R0.1, file = paste0(prefix, ".Weights.txt"), sep = "\t")

weightTab4db <- weightTab_R0.1[, c("rsid", "gene", "weight", "ref_allele", "counted_allele")]
colnames(weightTab4db) <- c("rsid", "gene", "weight", "ref_allele", "eff_allele")

extraTab4db <- extraTab_R0.1

### Build DB
library(sqldf)

db <- dbConnect(SQLite(), dbname = paste0(prefix, ".db"))

dbWriteTable(conn = db, name = "extra", value = extraTab4db, overwrite = TRUE)
dbWriteTable(conn = db, name = "weights", value = weightTab4db, overwrite = TRUE)

### Cov
covTab <- fread(paste0(prefix, ".Covs.tmp.txt"), header = F)
colnames(covTab) <- c("GENE", "RSID1", "RSID2", "VALUE")

fwrite(covTab, file = paste0(prefix, ".Covs.all.txt"), sep = "\t")

covTab_R0.1 <- subset(covTab, GENE %in% extraTab_R0.1$gene)
fwrite(covTab_R0.1, file = paste0(prefix, ".Covs.txt"), sep = "\t")
