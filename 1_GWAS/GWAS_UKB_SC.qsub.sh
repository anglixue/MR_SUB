#PBS -S /bin/bash
#PBS -N SC_GWAS_UKB
#PBS -A UQ-IMB-CNSG
#PBS -o /shares/compbio/Group-Yang/a.xue/proj/addiction/SC/stdout_gwas_SC
#PBS -e /shares/compbio/Group-Yang/a.xue/proj/addiction/SC/stderr_gwas_SC
#PBS -l select=1:ncpus=16:intel=True:mem=96GB
#PBS -l walltime=168:00:00

cd /shares/compbio/Group-Yang/a.xue/proj/addiction/

if [ ! -d "SC" ]; then
    mkdir SC
fi

cd /shares/compbio/Group-Yang/a.xue/proj/addiction/SC/
mv *.* std* /shares/compbio/Group-Yang/a.xue/trash/

bolt="/shares/compbio/Group-Yang/a.xue/bin/BOLT/bolt"
phen="/shares/compbio/Group-Yang/a.xue/proj/addiction/Addiction_Substance_QCed.phen"
exsnp="/gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EUR_impHRC/ukbEUR_imp_all_v2_imp_QC_HRC_maf0001.EXCLUDEvarList"

bim="/gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EUR_impHRC/ukbEUR_imp_chr{1:22}_v2_imp_QC_HRC.bim"
bed="/gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EUR_impHRC/ukbEUR_imp_chr{1:22}_v2_imp_QC_HRC.bed"
fam="/gpfs/gpfs01/polaris/Q0286/UKBiobank/v2EUR_impHRC/ukbEUR_imp_chr1_v2_imp_QC_HRC.fam"

covar="/gpfs/gpfs01/polaris/Q0286/uqaxue/phen/covar_sex_age_10PCs.txt"
table="/gpfs/gpfs01/polaris/Q0286/uqaxue/bin/BOLT/tables/LDSCORE.1000G_EUR.tab.gz"
map="/gpfs/gpfs01/polaris/Q0286/uqaxue/bin/BOLT/tables/genetic_map_hg19.txt.gz"
hm3="/gpfs/gpfs01/polaris/Q0286/uqaxue/LDprune/hm3_pruned.snplist"
j_hm3="/gpfs/gpfs01/polaris/Q0286/uqaxue/phen/ukbEURu_imp_v2_HM3_QC_R2_09.snplist"
out="/shares/compbio/Group-Yang/a.xue/proj/addiction/SC/ukbEUR_MAF1e-4_SC.bolt"

$bolt --bed=$bed --bim=$bim --fam=$fam --exclude $exsnp --phenoFile=$phen --phenoCol=SC --covarFile $covar --covarCol=Sex --qCovarCol=Age --qCovarCol=PC{1:10} --lmm --modelSnps=$j_hm3 --maxModelSnps=720000 --LDscoresFile=$table --geneticMapFile=$map --statsFile=$out --numThreads 16 >> SC_bolt.log 2>&1
