#PBS -S /bin/bash
#PBS -N addiction_UKB_diseases_GSMR2
#PBS -A UQ-IMB-CNSG
#PBS -o /home/uqaxue/90days/uqaxue/proj/addiction/GSMR/main/stdout_addiction_UKB_diseases
#PBS -e /home/uqaxue/90days/uqaxue/proj/addiction/GSMR/main/stderr_addiction_UKB_diseases
#PBS -l select=1:ncpus=1:intel=True:mem=50GB
#PBS -l walltime=8:00:00

cd /home/uqaxue/90days/uqaxue/proj/addiction/GSMR/main/
# this script will run GCTA-integrated GSMR

gcta="/home/uqaxue/90days/uqaxue/software/gcta64_20Mar2020"
exposure_file="../Exposure.list"
outcome_file="../Outcome.list"
output_path="addiction_UKB_diseases"
geno_file="/home/uqaxue/90days/uqaxue/proj/addiction/GSMR/gsmr_ref_UKB.txt"

$gcta --gsmr2-beta-2 --mbfile ${geno_file} --gsmr-file ${exposure_file} ${outcome_file} --gsmr-direction 0 --heidi-thresh 0.01 --heidi-filter-thresh 0 --gsmr-snp-min 1 --clump-r2 0.01 --out ${output_path} --effect-plot 


