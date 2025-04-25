library(data.table)
library(sqldf)

gencode <- as.data.frame(fread("/nobackup/sbcs/pingj2/TWAS/db/gencode.v26.GRCh38.genes.txt"))

setwd("/nobackup/sbcs/pingj2/EURASN_TWAS/models")

for (race in c("EUR", "ASN")) {
    for (trptEvt in c("Exp", "APA", "spT")) {
        print(paste0(race, ".", trptEvt))
        maxR_models <- fread(paste0(race, ".", trptEvt, ".maxR.models.tsv"))
        model_prefixes <- unique(maxR_models$maxName)

        extras <- list()
        weights <- list()
        for (prefix in model_prefixes) {
            method <- strsplit(prefix, "\\.")[[1]][1]
            proj <- strsplit(prefix, "\\.")[[1]][2]
            extras[[prefix]] <- fread(paste0(method, "/", race, "/", proj, "/", method, ".", proj, ".", race, ".Extras.txt"))
            weights[[prefix]] <- fread(paste0(method, "/", race, "/", proj, "/", method, ".", proj, ".", race, ".Weights.txt"))
            if (method == "PrediXcan") {
                extras[[prefix]]$gene <- sapply(extras[[prefix]]$gene, function(x) strsplit(x, "\\.")[[1]][1])
                weights[[prefix]]$gene <- sapply(weights[[prefix]]$gene, function(x) strsplit(x, "\\.")[[1]][1])
            }
            if (method == "spTWAS") {
                extras[[prefix]]$gene <- sapply(extras[[prefix]]$gene, function(x) {
                    tmp <- strsplit(x, ":")[[1]]
                    return(paste0(tmp[1], ":", tmp[2], ":", tmp[3], ":", tmp[5]))
                })

                weights[[prefix]]$gene <- sapply(weights[[prefix]]$gene, function(x) {
                    tmp <- strsplit(x, ":")[[1]]
                    return(paste0(tmp[1], ":", tmp[2], ":", tmp[3], ":", tmp[5]))
                })
            }
            weights[[prefix]]$varID <- paste0("chr", weights[[prefix]]$rsid, "_b38")
            extras[[prefix]]$alpha <- 0.5
        }

        maxR_extra <- subset(extras[[model_prefixes[1]]], gene %in% maxR_models$gene[which(maxR_models$maxName == model_prefixes[1])])
        maxR_weight <- subset(weights[[model_prefixes[1]]][, c("gene", "rsid", "varID", "ref_allele", "counted_allele", "weight")], gene %in% maxR_models$gene[which(maxR_models$maxName == model_prefixes[1])])

        for (n in 2:length(model_prefixes)) {
            maxR_extra <- rbind(maxR_extra, subset(extras[[model_prefixes[n]]], gene %in% maxR_models$gene[which(maxR_models$maxName == model_prefixes[n])]))
            maxR_weight <- rbind(maxR_weight, subset(weights[[model_prefixes[n]]][, c("gene", "rsid", "varID", "ref_allele", "counted_allele", "weight")], gene %in% maxR_models$gene[which(maxR_models$maxName == model_prefixes[n])]))
        }

        if (trptEvt == "Exp") {
            maxR_extra$geneid <- maxR_extra$gene
            maxR_extra <- merge(maxR_extra, gencode[, c("geneid", "genetype")], by = "geneid", all.x = T)
        } else if (trptEvt == "spT") {
            maxR_extra$geneid <- sapply(maxR_extra$gene, function(x) strsplit(strsplit(x, ":")[[1]][4], "\\.")[[1]][1])
            maxR_extra <- merge(maxR_extra, gencode[, c("geneid", "genetype")], by = "geneid", all.x = T)
        } else if (trptEvt == "APA") {
            maxR_extra <- merge(maxR_extra, gencode[, c("genename", "genetype")], by = "genename", all.x = T)
        }

        maxR_extra <- maxR_extra[, c("gene", "genename", "genetype", "alpha", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval")]
        colnames(maxR_extra) <- c("gene", "genename", "gene_type", "alpha", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval")
        fwrite(maxR_extra, file = paste0(race, ".", trptEvt, ".Extras.txt"), sep = "\t")

        maxR_weight <- maxR_weight[, c("gene", "rsid", "varID", "ref_allele", "counted_allele", "weight")]
        colnames(maxR_weight) <- c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")
        fwrite(maxR_weight, file = paste0(race, ".", trptEvt, ".Weights.txt"), sep = "\t")

        db <- dbConnect(SQLite(), dbname = paste0(race, ".", trptEvt, ".db"))
        dbWriteTable(conn = db, name = "extra", value = as.data.frame(maxR_extra), overwrite = TRUE)
        dbWriteTable(conn = db, name = "weights", value = as.data.frame(maxR_weight), overwrite = TRUE)
    }
}

race <- "AFR"

for (trptEvt in c("Exp", "APA", "spT")) {
    print(paste0(race, ".", trptEvt))
    method <- ""
    if (trptEvt == "Exp") {
        method <- "PrediXcan"
    } else if (trptEvt == "APA") {
        method <- "APA"
    } else if (trptEvt == "spT") {
        method <- "spTWAS"
    }

    extra <- fread(paste0(method, "/", race, "/Komen/", method, ".Komen.", race, ".Extras.txt"))
    weight <- fread(paste0(method, "/", race, "/Komen/", method, ".Komen.", race, ".Weights.txt"))

    if (method == "PrediXcan") {
        extra$gene <- sapply(extra$gene, function(x) strsplit(x, "\\.")[[1]][1])
        weight$gene <- sapply(weight$gene, function(x) strsplit(x, "\\.")[[1]][1])
    }
    if (method == "spTWAS") {
        extra$gene <- sapply(extra$gene, function(x) {
            tmp <- strsplit(x, ":")[[1]]
            return(paste0(tmp[1], ":", tmp[2], ":", tmp[3], ":", tmp[5]))
        })
        weight$gene <- sapply(weight$gene, function(x) {
            tmp <- strsplit(x, ":")[[1]]
            return(paste0(tmp[1], ":", tmp[2], ":", tmp[3], ":", tmp[5]))
        })
    }
    weight$varID <- paste0("chr", weight$rsid, "_b38")
    extra$alpha <- 0.5

    if (trptEvt == "Exp") {
        extra$geneid <- extra$gene
        extra <- merge(extra, gencode[, c("geneid", "genetype")], by = "geneid", all.x = T)
    } else if (trptEvt == "spT") {
        extra$geneid <- sapply(extra$gene, function(x) strsplit(strsplit(x, ":")[[1]][4], "\\.")[[1]][1])
        extra <- merge(extra, gencode[, c("geneid", "genetype")], by = "geneid", all.x = T)
    } else if (trptEvt == "APA") {
        extra <- merge(extra, gencode[, c("genename", "genetype")], by = "genename", all.x = T)
    }
    extra <- extra[, c("gene", "genename", "genetype", "alpha", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval")]
    colnames(extra) <- c("gene", "genename", "gene_type", "alpha", "n.snps.in.model", "pred.perf.R2", "pred.perf.pval", "pred.perf.qval")
    fwrite(extra, file = paste0(race, ".", trptEvt, ".Extras.txt"), sep = "\t")

    weight <- weight[, c("gene", "rsid", "varID", "ref_allele", "counted_allele", "weight")]
    colnames(weight) <- c("gene", "rsid", "varID", "ref_allele", "eff_allele", "weight")
    fwrite(weight, file = paste0(race, ".", trptEvt, ".Weights.txt"), sep = "\t")

    db <- dbConnect(SQLite(), dbname = paste0(race, ".", trptEvt, ".db"))
    dbWriteTable(conn = db, name = "extra", value = as.data.frame(extra), overwrite = TRUE)
    dbWriteTable(conn = db, name = "weights", value = as.data.frame(weight), overwrite = TRUE)
}
