#PBS -S /bin/bash
#PBS -N meta
#PBS -A UQ-IMB-CNSG
#PBS -o ./stdout_meta
#PBS -e ./stderr_meta
#PBS -l select=1:ncpus=1:intel=True:mem=10GB
#PBS -l walltime=5:00:00

cd /shares/compbio/Group-Yang/a.xue/proj/addiction/AC_not_ascertained_log2_adj_longitudinal_removed_underreport_all/


touch meta.log
Rscript meta.R >> meta.log




