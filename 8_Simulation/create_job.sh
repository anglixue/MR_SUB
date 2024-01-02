##
for i in {0.01,0.02,0.03,0.04,0.05,0.1}
do

# sed  's/VAR_P/'${i}'/g' simu-mr.R > simu-mr_${i}.R
sed  's/VAR_P/'${i}'/g' simu-mr_directional.R > simu-mr_directional_${i}.R
# sed  's/VAR_P/'${i}'/g' simu-mr_directional_oligo.R > simu-mr_directional_oligo_${i}.R


# sed  's/VAR_P/'${i}'/g' simulation_balanced.job > simulation_balanced_${i}.job
sed  's/VAR_P/'${i}'/g' simulation_directional.job > simulation_directional_${i}.job
# sed  's/VAR_P/'${i}'/g' simulation_directional_oligo.job > simulation_directional_oligo_${i}.job


#qsub simulation_directional_${i}.job
qsub simulation_directional_${i}.job
done
#


# for i in {0.03,0.05,0.1,0.15,0.3,0.4,0.5}; do rm ./Rsq_zp_${i}/results/data/*.txt -f; done
# for i in {0.03,0.05,0.1,0.15,0.3,0.4,0.5}; do rm ./Oligo_Rsq_zp_${i}/results/data/*.txt -f; done

# for i in {0.03,0.05,0.1,0.15,0.3,0.4,0.5}; do rm ./Rsq_zp_${i}/results_directional/data/*.txt -f; done
# for i in {0.03,0.05,0.1,0.15,0.3,0.4,0.5}; do rm ./Oligo_Rsq_zp_${i}/results_directional/data/*.txt -f; done

# for i in {0.03,0.05,0.1,0.15,0.3,0.4,0.5}; do rm ./Rsq_zp_${i}/std* -f ; rm ./Rsq_zp_${i}/*.log -f; done
# for i in {0.03,0.05,0.1,0.15,0.3,0.4,0.5}; do rm ./Oligo_Rsq_zp_${i}/std* -f ; rm ./Oligo_Rsq_zp_${i}/*.log -f; done





##
