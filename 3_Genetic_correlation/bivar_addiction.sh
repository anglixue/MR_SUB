
cd /shares/compbio/Group-Yang/a.xue/proj/addiction/ldsc/rg/

munge="/shares/compbio/Group-Yang/a.xue/app/ldsc/munge_sumstats.py"
ldsc="/shares/compbio/Group-Yang/a.xue/app/ldsc/ldsc.py"
panel="/shares/compbio/Group-Yang/a.xue/app/ldsc/eur_w_ld_chr/w_hm3.snplist"
ld="/shares/compbio/Group-Yang/a.xue/app/ldsc/eur_w_ld_chr/"

lables=(CI MCI HCI TI MTI HTI)
for i in `seq 0 5`; do
    for j in `seq 0 5`;do
        if [ "$j" -lt "$i" ]; then
#
$ldsc \
--rg ../${lables[$j]}.sumstats.gz,../${lables[$i]}.sumstats.gz \
--ref-ld-chr $ld \
--w-ld-chr $ld \
--out ./coffee_type/${lables[$j]}_${lables[$i]}_rg

        fi
    done
done



