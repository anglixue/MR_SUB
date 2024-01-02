library(data.table)
less=fread("../AC_less_log2_not_ascertained_removed_underreport_all/ukbEUR_AC_less_log2_not_ascertained_removed_underreport_all_cojo.txt",header=T);less=as.data.frame(less)
same=fread("../AC_same_log2_removed_underreport_all/ukbEUR_AC_same_log2_removed_underreport_all_cojo.txt",header=T);same=as.data.frame(same)
more=fread("../AC_more_log2_removed_underreport_all/ukbEUR_AC_more_log2_removed_underreport_all_cojo.txt",header=T);more=as.data.frame(more)

snp=Reduce(intersect,list(less$SNP,same$SNP,more$SNP))
less=less[match(snp,less$SNP),]
same=same[match(snp,same$SNP),]
more=more[match(snp,more$SNP),]

new=less
new$N=less$N+same$N+more$N
new$freq=(less$freq*less$N+same$freq*same$N+more$freq*more$N)/new$N

new$se=1/sqrt((1/less$se^2)+(1/same$se^2)+(1/more$se^2))
new$b=((less$b/less$se^2)+(same$b/same$se^2)+(more$b/more$se^2))/((1/less$se^2)+(1/same$se^2)+(1/more$se^2))
new$P=pchisq((new$b/new$se)^2,1,lower.tail=F)
new=na.omit(new)

print("Meta-analysis finished")
write.table(new,"ukbEUR_AC_not_ascertained_log2_adj_longitudinal_removed_underreport_all_cojo.txt",row.names=F,col.names=T,quote=F)

com=fread("../AC_less_log2_not_ascertained_removed_underreport_all/ukbEUR_AC_less_log2_not_ascertained_removed_underreport_all_common.txt",header=T)
com=as.data.frame(com)
com=com[match(new$SNP,com$SNP),]
res=cbind(com[,c(1,3)],new)
print("Writing common file")
write.table(res[,c(1,3,2,4,5,6,10,7:9)],"ukbEUR_AC_not_ascertained_log2_adj_longitudinal_removed_underreport_all_common.txt",row.names=F,col.names=T,quote=F)

print("Creating LDSC file")
res$Z=res$b/res$se
hm3=read.table("/shares/compbio/Group-Yang/a.xue/app/ldsc/eur_w_ld_chr/w_hm3.snplist",header=T)
snp=Reduce(intersect,list(hm3$SNP,res$SNP))

res=res[match(snp,res$SNP),]
write.table(res[,c("SNP","A1","A2","Z","N","P")],"ukbEUR_AC_not_ascertained_log2_adj_longitudinal_removed_underreport_all_ldsc.txt",row.names=F,col.names=T,quote=F)






