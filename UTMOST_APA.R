# module load PLINK/1.9b_5.2 GCC/10.2.0 OpenMPI/4.0.5 R/4.0.5

library(glmnet)
library(foreach)

set.seed(1024)
source('/home/pingj2/soft/glasso/glasso.r')

mod <- 1  #mod=1 training(ABCD) tuning&testing(E); mod=2 training(ABC) tuning(D) testing(E)

args <- commandArgs(trailingOnly = TRUE)
targetGene <- args[1] # NM_001301849.2|SGSM3|chr22|+
race <- args[2] ### EUR / ASN

geneid <- strsplit(targetGene, "\\|")[[1]][1]
genesymbol <- strsplit(targetGene, "\\|")[[1]][2]

#parameters setting
ntune <- 5  #N of grids for each tuning parameter
fold <- 5  #N of fold

# Load gencode annotation file
gencode_path <- "/nobackup/sbcs/pingj2/db/hg38_refseq_extracted_3UTR.bed"
gencode <- read.table(gencode_path, header = F, stringsAsFactors = F)

cat('INFO Starting UTMOST for APA', genesymbol, targetGene, 'in', race, 'at MOD', mod, 'with', ntune, 'tuning parameters and', fold, 'folds', '\n')

### functions ###
loadRData <- function (fileName) {
  load(fileName)
  get(ls()[ls() != "fileName"])
}

### some preset
plink_path <- 'plink'

genotype_path <- paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/geno/", race)
SNPs_in_sumstat <- paste0("/nobackup/sbcs/pingj2/sumstat/breast/", race, "/", race, ".variants")

###
workDir <- paste0("/nobackup/sbcs/pingj2/EURASN_TWAS/models/UTMOST_APA/", race, "/")
setwd(workDir)
tmp_folder <- paste0(workDir, "/tmp/")

# mkdir tmp folder. will be cleaned
cat('INFO mkdir tmp folders\n')
options(warn = -1)
dir.create(tmp_folder)
dir.create(paste0(tmp_folder, '/', geneid))
tmp_folder <- paste0(tmp_folder, '/', geneid)
options(warn = 0)

### Generate Exp table
cat('INFO loading expression data for gene', genesymbol, targetGene, '\n')

expdir <- "/nobackup/sbcs/pingj2/EURASN_TWAS/APA/"

projIDs <- c("Komen", "GTEx")
if (race == "ASN") {
    projIDs <- c("Komen", "ABCC")
} else if (race == "EUR") {
    projIDs <- c("Komen", "GTEx")
}

exp <- list()
for (n in 1:length(projIDs)) {
  projID <- projIDs[n]
  print(projID)
  exptab <- loadRData(paste0(expdir, "/", projID, "/", projID, ".", race, ".PDUI.residuals.rda"))
  if (targetGene %in% colnames(exptab)) {
    tmpExpTab <- as.data.frame(matrix(nrow = nrow(exptab), ncol = 2))
    colnames(tmpExpTab) <- c("sampleid", "exp")
    tmpExpTab[, 1] <- rownames(exptab)
    tmpExpTab[, 2] <- exptab[, targetGene]
    
    exp[[projID]] <- tmpExpTab
  }
}

expressedProjs <- names(exp)
P <- length(expressedProjs) #N of tissues/omics

ssize <- unlist(lapply(exp, nrow))  #N of samples for each tissue
T_num <- length(expressedProjs)  #N of tissues

### Prepare genotyping data
#get chr pos for the gene

dist <- 500000
chr <- as.numeric(sub('^...', '', gencode$V1[which(gencode$V4 == targetGene)]))
pos_from <- gencode$V2[which(gencode$V4 == targetGene)]
pos_to <- gencode$V3[which(gencode$V4 == targetGene)]
pos_from_gene <- max(1, pos_from - dist)
pos_to_gene <- pos_to + dist

cat('INFO generating dosage genotype data\n')
#extract genotypes from plink file to dosage file (1 mb)
cmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, 
              ' --chr ', chr, ' --from-bp ', pos_from_gene, ' --to-bp ', pos_to_gene, ' --recode A --out ', tmp_folder, '/', geneid)
system(cmd, ignore.stdout = T, ignore.stderr = T)

cmd <- paste0(plink_path,' --bfile ', genotype_path, ' --extract ', SNPs_in_sumstat, 
              ' --chr ', chr, ' --from-bp ', pos_from_gene, ' --to-bp ', pos_to_gene, ' --make-bed --out ', tmp_folder, '/', geneid)
system(cmd, ignore.stdout = T, ignore.stderr = T)

#load dosage file (1 mb)
dosage_500k <- try(read.table(paste0(tmp_folder, '/', geneid, '.raw'), header = T, stringsAsFactors = F))
if ('try-error' %in% class(dosage_500k)) {
  stop('no SNP available for this gene')
}

dosage_500k <- dosage_500k[, -c(1, 3:6)]
names(dosage_500k) <- c('SUBJID', gsub('.{2}$', '', sub('.', '', names(dosage_500k)[-1])))
dosage_500k[, -1] <- round(apply(dosage_500k[, -1], 2, function(x) ifelse(is.na(x), mean(x, na.rm = T), x)), 3) #post imputation imputation. NA replaced by mean

for(j in 2:ncol(dosage_500k)){ #center dosage to 0, not sure why
  dosage_500k[, j] <- dosage_500k[, j] - mean(dosage_500k[, j])
}

allsampleIDs <- c(exp[[1]]$sampleid, exp[[2]]$sampleid)
dosage_500k <- subset(dosage_500k, SUBJID %in% allsampleIDs)

N <- nrow(dosage_500k) #overall genotyped sample size

#generate covariance matrix
tmp <- as.matrix(dosage_500k[, -1])
XX <- t(tmp) %*% as.matrix(tmp)/N #var cor matrix
Xnorm <- diag(XX)  #var of dosage
remove(tmp)
remove(XX)
sub_id <- dosage_500k[, 1]  #subject id
M <- ncol(dosage_500k) - 1 #no of snps

#subject id map for each tissue and each fold
sub_id_map <- list()  
for (t in 1:T_num) {
  tmp <- rep(0, nrow(exp[[t]]))
  for (j in 1:length(tmp)) {
    tmp[j] <- which(sub_id == exp[[t]][j, 1])
  }
  sub_id_map[[t]] <- tmp
}

#cv setting
cv_config <- cv_helper(N, fold) 
cv_perm <- cv_config$perm
cv_idx <- cv_config$idx

#result boxes
single_res_test <- list()
single_lam <- matrix(0, fold, P)
single_theta_est <- list()

multi_res_test <- list()
multi_lam <- matrix(0, fold, 2)
multi_theta_est <- list()

multi_res_test2 <- list()
multi_lam2 <- array(0, dim = c(fold, P, 2))
multi_theta_est2 <- list()

res_tune <- list()
rec_lamv <- matrix(0, fold, ntune)

avg_tune_res <- list() #average tuning error
single_initial_est <- list()

#--------------get start--------------
cat('INFO Start model training for gene', genesymbol, '\n')
for (f in 1:fold) {
  cat('INFO Fold =', f, '\n')
  test_index <- cv_perm[cv_idx[f, 1]:cv_idx[f,2]]
  test_id <- sub_id[test_index] #sample id in the testing set
  tuning_index <- cv_perm[cv_idx[f %% fold + 1, 1]:cv_idx[f %% fold + 1, 2]]
  tuning_id <- sub_id[tuning_index] #sample in the 'early stop' set
  
  #dfs for each set
  X_test <- list()
  Y_test <- list()
  X_tune <- list()
  Y_tune <- list()
  X_train <- list()
  Y_train <- list()
  X_all <- list()
  Y_all <- list()
  
  for (t in 1:T_num) {
    X_train_tmp <- sub_id_map[[t]][!(sub_id_map[[t]] %in% c(test_index, tuning_index))]
    Y_train_tmp <- !(sub_id_map[[t]] %in% c(test_index, tuning_index))
    X_tuning_tmp <- sub_id_map[[t]][(sub_id_map[[t]] %in% tuning_index)]
    Y_tuning_tmp <- (sub_id_map[[t]] %in% tuning_index)
    X_test_tmp <- sub_id_map[[t]][(sub_id_map[[t]] %in% test_index)]
    Y_test_tmp <- (sub_id_map[[t]] %in% test_index)
    X_train[[t]] <- apply(as.matrix(dosage_500k[X_train_tmp, -1]), 2, as.numeric)
    Y_train[[t]] <- exp[[t]][Y_train_tmp, 2]
    X_tune[[t]] <- apply(as.matrix(dosage_500k[X_tuning_tmp, -1]), 2, as.numeric)
    Y_tune[[t]] <- exp[[t]][Y_tuning_tmp, 2]
    X_test[[t]] <- apply(as.matrix(dosage_500k[X_test_tmp, -1]), 2, as.numeric)
    Y_test[[t]] <- exp[[t]][Y_test_tmp, 2]
    X_all_tmp <- sub_id_map[[t]]
    X_all[[t]] <- apply(as.matrix(dosage_500k[X_all_tmp, -1]), 2, as.numeric)
    Y_all[[t]] <- exp[[t]][,2]
    
    if (mod == 1) {
      #combine
      X_train[[t]] <- rbind(X_train[[t]], X_test[[t]])
      Y_train[[t]] <- c(Y_train[[t]], Y_test[[t]])
    }
  }
  
  #---single tissue model training---
  #model training using all of the samples to get lambda range and initial est
  if (f == 1) {
    single_initial_est_all <- matrix(0, ncol(X_train[[1]]), T_num) #N of snp by N tissue matrix
    single_summary_all <- list()
    for(t in 1:T_num){
      tt <- cv.glmnet(X_all[[t]], Y_all[[t]], alpha = 0.5, nfolds = 5, nlambda = 50, dfmax = 100) #cv in 100% of the data, do not find the lambda here
      single_summary_all[[t]] <- tt
      single_initial_est_all[, t] <- tt$glmnet.fit$beta[, which.min(tt$cvm)] #single tissue EN, beta and lambda with best cv error
    }
  }
  
  #single tissue elastic net for each fold
  single_initial_est[[f]] <- matrix(0, ncol(X_train[[1]]), T_num) #snp by tissue matrix
  
  for (t in 1:T_num) {
    tt <- cv.glmnet(X_train[[t]], Y_train[[t]], alpha = 0.5, nfolds = 5, nlambda = 50, dfmax = 100) 
    single_initial_est[[f]][, t] <- tt$glmnet.fit$beta[, which.min(tt$cvm)]
  }
  
  #---joint model training---
  # use elastic net ests row norm as weights 
  lam_range <- minmax_lambda(single_summary_all) #lambda range
  sig_norm <- apply(single_initial_est[[f]], 1, function(x) sqrt(sum(x^2)) ) #norm for each snp
  sig_norm[sig_norm == 0] <- rep(min(sig_norm[sig_norm > 0]), sum(sig_norm == 0))/2
  sig_norm <- sig_norm/sum(sig_norm) #??? length=N of snps
  weights2 <- 1/sig_norm
  weights2 <- weights2/sum(weights2)
  
  tis_norm <- apply(single_initial_est[[f]], 2, function(x) sum(abs(x)) )
  tis_norm[tis_norm == 0] <- rep(min(tis_norm[tis_norm > 0]), sum(tis_norm == 0))/2
  tis_norm <- tis_norm/sum(tis_norm)
  weights1 <- 1/tis_norm
  weights1 <- weights1/sum(weights1)
  lam_V <- seq(lam_range[1], lam_range[2], length.out = ntune)
  
  initial_numeric <- as.numeric(single_initial_est[[f]])
  
  ## preparation
  XY <- grad_prep(X_train, Y_train)
  XX_train <- lapply(X_train, function(x) t(x) %*% x / nrow(x) )
  spsz <- unlist(lapply(X_train, nrow)) #N of samples for each training set
  res_tune[[f]] <- array(NA, dim = c(ntune, ntune, P)) #5*5 matrix *49 tissue, to find the best comb of lam1 and lam2
  
  rec_lamv[f,] <- lam_V  #should be the same across 5 folds
  
  #go throught all of the combination of lam1 and lam2
  for (lam1 in 1:ntune) {
    for (lam2 in 1:ntune) {
      #print(paste0('lam1= ',lam1,' lam2= ',lam2))
      single_est <- matrix(initial_numeric, M, P) #snps by tissue matrix. initial iteration with single tissue estimates, save compute time
      ans <- glasso(X = X_train, Y = Y_train, X1 = X_tune, Y1 = Y_tune, XX = XX_train, XY = XY, Xnorm = Xnorm, 
                    lambda1 = lam_V[lam1]/spsz, lambda2 = lam_V[lam2], theta = single_est)
      #ans$est: snps by tissue matrix
      #ans$tune_err: tuning error for each tissue
      #ans$avg_tune_err: average tuning error among 49 tissues
      
      #mod 1: apply to tuning set / testing set
      #mod 2: apply to testing set 
      if (lam1 == 1 & lam2 == 1) {
        multi_res_test[[f]] <- list()
      }
      if (mod == 1) {
        multi_res_test[[f]][[paste0(lam1, '_', lam2)]] <- multi_mse(ans$est, X_tune, Y_tune)
      } else {
        multi_res_test[[f]][[paste0(lam1, '_', lam2)]] <- multi_mse(ans$est, X_test, Y_test)
      }

      if (sum(ans$est != 0) > 0) {
        res_tune[[f]][lam1, lam2, ] <- ans$tune_err #tuning err for each lam pair
        remove(single_est)
        remove(ans)
      } else {
        res_tune[[f]][lam1, lam2, ] <- ans$tune_err
        remove(single_est)
        remove(ans)
      }
    }
  }
  #average tuning error across the tissues for each lambda pair
  avg_tune_res[[f]] <- apply(res_tune[[f]], c(1,2), mean) 
}

#-----------------------------
#find the best comb of lambda1 and lambda2 across 5 folds
avg_tune_cross_fold <- matrix(rep(1, ntune * ntune), nrow = ntune, ncol = ntune)
for (i_tune in 1:ntune) {
  for (j_tune in 1:ntune) {
    avg_tune_cross_fold[i_tune, j_tune] <- mean(c(avg_tune_res[[1]][i_tune, j_tune],
                                                  avg_tune_res[[2]][i_tune, j_tune],
                                                  avg_tune_res[[3]][i_tune, j_tune],
                                                  avg_tune_res[[4]][i_tune, j_tune],
                                                  avg_tune_res[[5]][i_tune, j_tune]), na.rm = T)
    avg_tune_cross_fold[i_tune, j_tune] <- ifelse(is.na(avg_tune_cross_fold[i_tune, j_tune]), 1, avg_tune_cross_fold[i_tune, j_tune])
  }
}

#find the best comb for each fold
best.lam <- which(avg_tune_cross_fold == min(avg_tune_cross_fold), arr.ind = TRUE)[1, ]

cat('INFO The best lambda combination for each fold is', paste0(as.numeric(best.lam)[1], '_', as.numeric(best.lam)[2]), '\n')

#predicted expression in tuning set under best lambda comb
multi_res <- list()
for (f in 1:ntune) {
  multi_res[[f]] <- multi_res_test[[f]][[paste0(as.numeric(best.lam)[1], '_', as.numeric(best.lam)[2])]]
}

#combine 5 fold results
cv_df <- list()
cv_r <- cv_p <- c()

for (proj_i in 1:T_num) {
  cv_df[[proj_i]] <- multi_res[[1]][proj_i][[1]]
  
  for (fold_i in 2:fold) {
    cv_df[[proj_i]] <- rbind(cv_df[[proj_i]], multi_res[[fold_i]][proj_i][[1]])
  }
  #print(nrow(cv_df[[proj_i]]))
  fit <- cor.test(cv_df[[proj_i]][, 1], cv_df[[proj_i]][, 2])
  cv_r[proj_i] <- as.numeric(fit$estimate)
  cv_p[proj_i] <- as.numeric(fit$p.value)
  
  cv_r[proj_i] <- ifelse(is.na(cv_r[proj_i]), 0, cv_r[proj_i])
  cv_p[proj_i] <- ifelse(is.na(cv_p[proj_i]), 1, cv_p[proj_i])
}

#generate an estimate with whole data 
cat('INFO Model training using entire data \n')
XY <- grad_prep(X_all, Y_all)
XX_all <- lapply(X_all, function(x) t(x) %*% x / nrow(x) )
ans <- glasso_no_early_stopping(X = X_all, Y = Y_all, XX = XX_all, XY = XY, Xnorm = Xnorm,
                                lambda1 = lam_V[best.lam[1]]/spsz, lambda2 = lam_V[best.lam[2]], theta = single_initial_est_all)

info <- read.table(paste0(tmp_folder, '/', geneid, '.bim'), header = F, stringsAsFactors = F)
downstream_est <- data.frame(info[, c(2, 5, 6)], ans$est)
colnames(downstream_est)[4:ncol(downstream_est)] <- expressedProjs

#-------------------------------------------------------------
#calculate r and p (retrained)
pred_exp_list <- list()
for (proj_i in 1:length(expressedProjs)) {
  projName <- expressedProjs[proj_i]
  pred_exp_list[[proj_i]] <- t(apply(X_all[[proj_i]], function(x) x * downstream_est[, proj_i + 3], MARGIN = 1))
  cor_result <- cor.test(rowSums(pred_exp_list[[proj_i]]), Y_all[[proj_i]]) #retrained
  
  weight_df <- data.frame(gene = targetGene,
                          rsid = downstream_est[, 1],
                          chr_bp = paste0(info$V1[which(info$V2 %in% downstream_est[, 1])], "_", info$V4[which(info$V2 %in% downstream_est[, 1])]),
                          ref_allele = downstream_est[, 2],
                          counted_allele = downstream_est[, 3],
                          weight = downstream_est[, proj_i + 3],
                          r_cv = cv_r[proj_i],
                          r2_cv = (cv_r[proj_i]) ^ 2,  #cross-validation r2
                          p_cv = cv_p[proj_i], #cross-validation p
                          lambda = paste0(round(lam_V[best.lam[1]], 2), ';', round(lam_V[best.lam[2]], 2)),
                          r_rt = as.numeric(cor_result$estimate), #retrain r
                          r2_rt = as.numeric(cor_result$estimate ^ 2), #retrain r2
                          p_rt = cor_result$p.value, #retrain p
                          iter = ans$iter) # tracking the iterations

  #rm weight=0
  weight_df <- weight_df[weight_df$weight != 0, ]

  #if (nrow(weight_df) > 0 & cv_r[proj_i] > 0.1) { # & cv_p[proj_i] < 0.05
  if (nrow(weight_df) > 0) {
    #print(proj_i)
    write.table(weight_df, paste0(workDir, "res/Weight_", projName, '_', geneid, '.txt'), quote = F,row.names = F, sep = '\t')

    extra_res <- t(c(targetGene, genesymbol, cv_r[proj_i], (cv_r[proj_i]) ^ 2, length(unique(weight_df$rsid)), cv_p[proj_i]))
    colnames(extra_res) <- c("gene", "genename", "pred.perf.R", "pred.perf.R2", "n.snps.in.model", "pred.perf.pval")
    write.table(extra_res, file = paste0(workDir, "res/Extra_", projName, "_", geneid, ".txt"), sep = "\t", row.name = F, col.names = T, quote = F)

    #estimate covariance matrix
    if ( length(unique(weight_df$rsid)) == 1 ) {
        dosage <- as.data.frame(t(t(dosage_500k[, weight_df$rsid])))
        colnames(dosage) <- unique(weight_df$rsid)
    } else {
        dosage <- dosage_500k[, weight_df$rsid]
    }
    
    cov <- cov(as.matrix(dosage))
    cov[upper.tri(cov)] <- NA
    cov <- cbind(expand.grid(dimnames(cov)), value = as.vector(cov))
    colnames(cov) <- c('RSID1', 'RSID2', 'VALUE')
    cov <- cov[!is.na(cov$VALUE), ]
    cov$GENE <- targetGene
    cov <- cov[, c('GENE', 'RSID1', 'RSID2', 'VALUE')]
    #save covariance matrix
    write.table(cov, paste0(workDir, 'res/Cov_', projName, '_', geneid, '.txt'), quote = F, row.names = F, sep = '\t')
  }
}

#cleaning
cat('INFO cleaning the tmp folder ... \n')
#will only clean the subfolder under the tmp folder
system(paste0('rm -r ', tmp_folder), wait = T)

cat('INFO UTMOST for', genesymbol, targetGene, 'done \n')
