#PBS -S /bin/bash
#PBS -N gsmr_addiction_public_disease_socio
#PBS -A UQ-IMB-CNSG
#PBS -o /shares/compbio/Group-Yang/a.xue/proj/addiction/GSMR/new_global_HEIDI/stdout_gsmr_addiction_public_disease_socio
#PBS -e /shares/compbio/Group-Yang/a.xue/proj/addiction/GSMR/new_global_HEIDI/stderr_gsmr_addiction_public_disease_socio
#PBS -l select=1:ncpus=1:intel=True:mem=20GB
#PBS -l walltime=24:00:00

cd /shares/compbio/Group-Yang/a.xue/proj/addiction/GSMR/new_global_HEIDI/

# this script will run GCTA-integrated GSMR

gcta="/shares/compbio/Group-Yang/uqzzhu1/bin/gcta64_20Mar2020"
exposure_file="../Exposure.list"
outcome_file="../Outcome_public_disease_socio.list"
output_path="addiction_public_disease_socio"
geno_file="/shares/compbio/Group-Yang/a.xue/proj/addiction/GSMR/gsmr_ref_UKB.txt"


${gcta} --gsmr2-beta-2 --mbfile ${geno_file}  --gsmr-file ${exposure_file} ${outcome_file} --gsmr-direction 0 --heidi-thresh 0.01 --heidi-filter-thresh 0.005 --gsmr-snp-min 5 --clump-r2 0.01 --out ${output_path} --effect-plot 



