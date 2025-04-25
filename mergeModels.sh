### UTMOST EUR Komen

method=UTMOST
race=EUR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST EUR GTEx

method=UTMOST
race=EUR
projID=GTEx

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST EUR
tar zcf ${method}.${race}.res.tar.gz res &
rm *.tmp.txt

### UTMOST ASN Komen

method=UTMOST
race=ASN
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST ASN ABCC

method=UTMOST
race=ASN
projID=ABCC

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST ASN
tar zcf ${method}.${race}.res.tar.gz res &
rm *.tmp.txt

### PrediXcan EUR Komen

method=PrediXcan
race=EUR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### PrediXcan EUR GTEx

method=PrediXcan
race=EUR
projID=GTEx

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### PrediXcan ASN Komen

method=PrediXcan
race=ASN
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### PrediXcan ASN ABCC

method=PrediXcan
race=ASN
projID=ABCC

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### JTI EUR GTEx

method=JTI
race=EUR
projID=GTEx

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### PrediXcan AFR Komen

method=PrediXcan
race=AFR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### APA EUR Komen

method=APA
race=EUR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### APA AFR Komen

method=APA
race=AFR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### APA ASN Komen

method=APA
race=ASN
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### APA EUR GTEx

method=APA
race=EUR
projID=GTEx

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### APA ASN ABCC

method=APA
race=ASN
projID=ABCC

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### spTWAS EUR Komen

method=spTWAS
race=EUR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### spTWAS AFR Komen

method=spTWAS
race=AFR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### spTWAS ASN Komen

method=spTWAS
race=ASN
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### spTWAS ASN ABCC

method=spTWAS
race=ASN
projID=ABCC

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### spTWAS EUR GTEx

method=spTWAS
race=EUR
projID=GTEx

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}/${projID}

find res/ -maxdepth 1 -type f -name 'Extra_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Weight_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name 'Cov_*' -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

rm -r *.tmp.txt tmp &
tar zcf ${method}.${projID}.${race}.res.tar.gz res &

### UTMOST APA EUR Komen

method=UTMOST_APA
race=EUR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### APA UTMOST EUR GTEx

method=UTMOST_APA
race=EUR
projID=GTEx

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${projID}/${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${projID}/${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${projID}/${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST EUR
tar zcf ${method}.${race}.res.tar.gz res &
rm *.tmp.txt

### APA UTMOST ASN Komen

method=UTMOST_APA
race=ASN
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${projID}/${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${projID}/${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${projID}/${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### APA UTMOST ASN ABCC

method=UTMOST_APA
race=ASN
projID=ABCC

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${projID}/${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${projID}/${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${projID}/${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST ASN
tar zcf ${method}.${race}.res.tar.gz res &
rm *.tmp.txt

### UTMOST spTWAS EUR Komen

method=UTMOST_spTWAS
race=EUR
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${projID}/${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${projID}/${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${projID}/${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### spTWAS UTMOST EUR GTEx

method=UTMOST_spTWAS
race=EUR
projID=GTEx

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${projID}/${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${projID}/${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${projID}/${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST EUR
tar zcf ${method}.${race}.res.tar.gz res &
rm *.tmp.txt

### spTWAS UTMOST ASN Komen

method=UTMOST_spTWAS
race=ASN
projID=Komen

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${projID}/${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${projID}/${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${projID}/${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### APA UTMOST ASN ABCC

method=UTMOST_spTWAS
race=ASN
projID=ABCC

cd /nobackup/sbcs/pingj2/EURASN_TWAS/models/${method}/${race}

find res/ -maxdepth 1 -type f -name "Extra_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Extras.tmp.txt &
find res/ -maxdepth 1 -type f -name "Weight_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Weights.tmp.txt &
find res/ -maxdepth 1 -type f -name "Cov_${projID}_*" -print0 | sort -z | xargs -0 cat -- >> ${projID}/${method}.${projID}.${race}.Covs.tmp.txt &

sed -i '/genename/d' ${projID}/${method}.${projID}.${race}.Extras.tmp.txt
sed -i '/gene/d' ${projID}/${method}.${projID}.${race}.Weights.tmp.txt
sed -i '/GENE/d' ${projID}/${method}.${projID}.${race}.Covs.tmp.txt

Rscript /nobackup/sbcs/pingj2/EURASN_TWAS/scripts/mergeModels.R ${method} ${race} ${projID}

### UTMOST ASN
tar zcf ${method}.${race}.res.tar.gz res &
rm *.tmp.txt
