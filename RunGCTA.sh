###Run GCTA scripts
#!/bin/bash

RACE=$1 # EUR, AFR, ASN

for SUBTYPE in Overall ERpos ERneg
do
    gctaDir=/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA/${SUBTYPE}/${RACE}

    mkdir ${gctaDir}/res/
    mkdir ${gctaDir}/tmp/

    ### 
    ml load PLINK/2.00-alpha2

    plink2 --bfile /nobackup/sbcs/pingj2/1KG/${RACE}/1KG.${RACE}.hg38 \
        --extract ${gctaDir}/all_known_loci.500k.snps \
        --make-bed \
        --out ${gctaDir}/tmp/all_known_loci

    while read line
    do
        arrIN=(${line//\t/ })
        chr=${arrIN[0]}
        pos=${arrIN[1]}
        rsid=${arrIN[3]}

        grep ${chr}_${pos} ${gctaDir}/all_known_loci.500k.ma | awk '{print $1}' >> ${gctaDir}/all_known_loci.cond.snplist
    done < /nobackup/sbcs/pingj2/db/BC_known_loci_hg38.cyto.txt

    for chr in {1..22}
    do
        /nobackup/sbcs/pingj2/soft/gcta_1.93.2beta/gcta64 \
          --bfile ${gctaDir}/tmp/all_known_loci \
          --chr ${chr} \
          --cojo-cond ${gctaDir}/all_known_loci.cond.snplist \
          --cojo-file ${gctaDir}/all_known_loci.500k.ma \
          --cojo-collinear 0.8 \
          --thread-num 36 \
          --out ${gctaDir}/res/all_known_loci_chr${chr}
    done
done

### Combine GCTA results
RACE=$1 # EUR, AFR, ASN

for SUBTYPE in Overall ERpos ERneg
do
    gctaDir=/nobackup/sbcs/pingj2/EURASN_TWAS/GCTA
    resDir=${gctaDir}/${SUBTYPE}/${RACE}/res
    
    cat ${resDir}/all_known_loci_chr1.cma.cojo > ${gctaDir}/${RACE}.${SUBTYPE}.COJO.cma.txt
    for chr in {2..22}
    do
        tail -n +2 ${resDir}/all_known_loci_chr${chr}.cma.cojo >> ${gctaDir}/${RACE}.${SUBTYPE}.COJO.cma.txt
    done
done
