#!/bin/sh

ml load GCC/8.2.0 OpenMPI/3.1.4 pandas/0.24.2 R/3.6.0

### ABCC ASN
race=ASN
projID=ABCC

gctaDir=/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/

for model in UTMOST_spTWAS UTMOST_APA APA spTWAS UTMOST PrediXcan
do
  echo ${model}

  resDir=${gctaDir}/assoRes_cond/${race}

  dbDir=/nobackup/sbcs/pingj2/EURASN_TWAS/models/${model}/${race}/${projID}

  twasDB=${dbDir}/${model}.${projID}.${race}.db
  twasCov=${dbDir}/${model}.${projID}.${race}.Covs.txt

  for subtype in Overall ERpos ERneg
  do
    echo ${subtype}
    /nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
      --model_db_path ${twasDB} \
      --covariance ${twasCov} \
      --gwas_folder ${gctaDir} \
      --gwas_file_pattern "${race}.${subtype}.COJO.sumstat.txt" \
      --snp_column SNP \
      --effect_allele_column TEST \
      --non_effect_allele_column OTHER \
      --beta_column BETA \
      --pvalue_column P \
      --keep_non_rsid \
      --output_file ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv 2>&1 | tee ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.log

    Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.tsv
  done
done

### Komen ASN
race=ASN
projID=Komen

gctaDir=/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/

for model in UTMOST_spTWAS UTMOST_APA APA spTWAS UTMOST PrediXcan
do
  echo ${model}

  resDir=${gctaDir}/assoRes_cond/${race}

  dbDir=/nobackup/sbcs/pingj2/EURASN_TWAS/models/${model}/${race}/${projID}

  twasDB=${dbDir}/${model}.${projID}.${race}.db
  twasCov=${dbDir}/${model}.${projID}.${race}.Covs.txt

  for subtype in Overall ERpos ERneg
  do
    echo ${subtype}
    /nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
      --model_db_path ${twasDB} \
      --covariance ${twasCov} \
      --gwas_folder ${gctaDir} \
      --gwas_file_pattern "${race}.${subtype}.COJO.sumstat.txt" \
      --snp_column SNP \
      --effect_allele_column TEST \
      --non_effect_allele_column OTHER \
      --beta_column BETA \
      --pvalue_column P \
      --keep_non_rsid \
      --output_file ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv 2>&1 | tee ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.log

    Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.tsv
  done
done

### Komen EUR
race=EUR
projID=Komen

gctaDir=/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/

for model in UTMOST_spTWAS UTMOST_APA APA spTWAS UTMOST PrediXcan
do
  echo ${model}

  resDir=${gctaDir}/assoRes_cond/${race}

  dbDir=/nobackup/sbcs/pingj2/EURASN_TWAS/models/${model}/${race}/${projID}

  twasDB=${dbDir}/${model}.${projID}.${race}.db
  twasCov=${dbDir}/${model}.${projID}.${race}.Covs.txt

  for subtype in Overall ERpos ERneg
  do
    echo ${subtype}
    /nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
      --model_db_path ${twasDB} \
      --covariance ${twasCov} \
      --gwas_folder ${gctaDir} \
      --gwas_file_pattern "${race}.${subtype}.COJO.sumstat.txt" \
      --snp_column SNP \
      --effect_allele_column TEST \
      --non_effect_allele_column OTHER \
      --beta_column BETA \
      --pvalue_column P \
      --keep_non_rsid \
      --output_file ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv 2>&1 | tee ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.log

    Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.tsv
  done
done

### Komen AFR
race=AFR
projID=Komen

gctaDir=/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/

for model in APA spTWAS PrediXcan
do
  echo ${model}

  resDir=${gctaDir}/assoRes_cond/${race}

  dbDir=/nobackup/sbcs/pingj2/EURASN_TWAS/models/${model}/${race}/${projID}

  twasDB=${dbDir}/${model}.${projID}.${race}.db
  twasCov=${dbDir}/${model}.${projID}.${race}.Covs.txt

  for subtype in Overall ERpos ERneg
  do
    echo ${subtype}
    /nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
      --model_db_path ${twasDB} \
      --covariance ${twasCov} \
      --gwas_folder ${gctaDir} \
      --gwas_file_pattern "${race}.${subtype}.COJO.sumstat.txt" \
      --snp_column SNP \
      --effect_allele_column TEST \
      --non_effect_allele_column OTHER \
      --beta_column BETA \
      --pvalue_column P \
      --keep_non_rsid \
      --output_file ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv 2>&1 | tee ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.log

    Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.tsv
  done
done

### GTEx EUR
race=EUR
projID=GTEx

gctaDir=/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/

for model in UTMOST PrediXcan JTI UTMOST_spTWAS UTMOST_APA APA spTWAS
do
  echo ${model}

  resDir=${gctaDir}/assoRes_cond/${race}

  dbDir=/nobackup/sbcs/pingj2/EURASN_TWAS/models/${model}/${race}/${projID}

  twasDB=${dbDir}/${model}.${projID}.${race}.db
  twasCov=${dbDir}/${model}.${projID}.${race}.Covs.txt

  for subtype in Overall ERpos ERneg
  do
    echo ${subtype}
    /nobackup/sbcs/pingj2/soft/MetaXcan/software/SPrediXcan.py \
      --model_db_path ${twasDB} \
      --covariance ${twasCov} \
      --gwas_folder ${gctaDir} \
      --gwas_file_pattern "${race}.${subtype}.COJO.sumstat.txt" \
      --snp_column SNP \
      --effect_allele_column TEST \
      --non_effect_allele_column OTHER \
      --beta_column BETA \
      --pvalue_column P \
      --keep_non_rsid \
      --output_file ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv 2>&1 | tee ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.log

    Rscript /nobackup/sbcs/pingj2/TWAS/scripts/addFDRBonf.R ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.csv ${resDir}/${model}.${projID}.${race}.${race}.${subtype}.assocRes.tsv
  done
done
