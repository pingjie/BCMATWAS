library(data.table)
library(dplyr)

knownloci <- fread("/nobackup/sbcs/pingj2/db/BC_known_loci_hg38.cyto.txt")
colnames(knownloci)[c(1:2)] <- c("CHR", "POS")

races <- c("EUR", "AFR", "ASN")
subtypes <- c("Overall", "ERpos", "ERneg")

for (race in races) {
    print(race)
    for (subtype in subtypes) {
        print(subtype)
        sumstat <- fread(paste0("/nobackup/sbcs/pingj2/sumstat/breast/", race, "/", race, "_", subtype, "_sumstat.txt"))
        sumstat$CHR <- as.numeric(sapply(sumstat$SNP, function(x) strsplit(x, "_")[[1]][1]))
        sumstat$POS <- as.numeric(sapply(sumstat$SNP, function(x) strsplit(x, "_")[[1]][2]))

        snp_in_region_all <- sumstat %>% filter(CHR == knownloci$CHR[1], between(POS, knownloci$POS[1] - 500000, knownloci$POS[1] + 500000))

        for (s in 2:nrow(knownloci)) {
            #print(knownloci$SNP[s])
            snp_in_region <- sumstat %>% filter(CHR == knownloci$CHR[s], between(POS, knownloci$POS[s] - 500000, knownloci$POS[s] + 500000))
            snp_in_region_all <- rbind(snp_in_region_all, snp_in_region)
        }

        snp_in_region_all <- snp_in_region_all %>% distinct(SNP, .keep_all = T)
        write.table(file = paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/", subtype, "/", race, "/", 'all_known_loci.500k.snps'),
                    snp_in_region_all$SNP, quote = F, col.names = F, row.names = F)
        
        snpMa <- snp_in_region_all[, c("SNP", "TEST", "OTHER", "TEST_FREQ", "BETA", "SE", "P", "N")]
        colnames(snpMa) <- c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")

        write.table(file = paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/", subtype, "/", race, "/", 'all_known_loci.500k.ma'), 
                    snpMa, sep = "\t", quote = F, col.names = T, row.names = F)
    }
}

### Generate a conditional sumstat

library(data.table)

races <- c("EUR", "AFR", "ASN")
subtypes <- c("Overall", "ERpos", "ERneg")

for (race in races) {
    print(race)
    for (subtype in subtypes) {
        print(subtype)
        sumstat <- fread(paste0("/nobackup/sbcs/pingj2/sumstat/breast/", race, "/", race, "_", subtype, "_sumstat.txt"))
        
        cojo_sumstat <- fread(paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/", race, ".", subtype, ".COJO.cma.txt"))
        cojo_sumstat <- cojo_sumstat[, c("SNP", "bC", "bC_se", "pC")]

        sumstat_conditional <- subset(sumstat, SNP %in% cojo_sumstat$SNP)
        sumstat_conditional <- merge(sumstat_conditional, cojo_sumstat, by = "SNP")
        sumstat_conditional <- sumstat_conditional[, c("SNP", "bC", "bC_se", "pC", "TEST", "OTHER", "TEST_FREQ", "N")]
        colnames(sumstat_conditional) <- colnames(sumstat)

        sumstat_other <- subset(sumstat, !(SNP %in% cojo_sumstat$SNP))
        sumstat_new <- rbind(sumstat_other, sumstat_conditional)

        fwrite(sumstat_new, file = paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/", race, ".", subtype, ".COJO.sumstat.txt"), 
               sep = "\t", quote = F, col.names = T, row.names = F)
    }
}
