
cd /shares/compbio/Group-Yang/a.xue/proj/addiction/ldsc/rg/

munge="/shares/compbio/Group-Yang/a.xue/app/ldsc/munge_sumstats.py"
ldsc="/shares/compbio/Group-Yang/a.xue/app/ldsc/ldsc.py"
panel="/shares/compbio/Group-Yang/a.xue/app/ldsc/eur_w_ld_chr/w_hm3.snplist"
ld="/shares/compbio/Group-Yang/a.xue/app/ldsc/eur_w_ld_chr/"

#t1=(AC,AC1,AC5,AUD)
#t2=(EA,HI)

for j in {CI,MCI,HCI,TI,MTI,HTI}

do

for i in {ASTHMA,ALLERGIC_RHINITIS,CARD,CANCER,DEPRESS,DIA2,DYSLIPID,HYPER,HEMORRHOIDS,HERNIA_ABDOMINOPELVIC,IRON_DEFICIENCY,IRRITABLE_BOWEL,OSTIOA,OSTIOP,PVD,PEPTIC_ULCERS,PSYCHIATRIC,VARICOSE_VEINS,CASES,SUM_OF_CASES}
do

$ldsc \
--rg ../${i}.sumstats.gz,../${j}.sumstats.gz \
--ref-ld-chr $ld \
--w-ld-chr $ld \
--out ./BMI/${i}_${j}_rg

    done
done



