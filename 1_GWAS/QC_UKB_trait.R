########## 22 diseases GWAS summary QC ##########

## Read the file
PATH="/shares/compbio/Group-Yang/a.xue/proj/addiction/trait/"

library(data.table)

trait=fread(paste(PATH,"ukbEUR_MAF1e-4_trait.bolt",sep=""),header=T)
trait=as.data.frame(trait)
trait=na.omit(trait)

## Get sample size N
# phen=read.table("/shares/compbio/Group-Yang/a.xue/proj/addiction/Addiction_Substance_QCed.phen",header=T)
phen=read.table("/gpfs/gpfs01/polaris/Q0286/uqaxue/phen/addiction/trait_zscore_normal_adj_by_sex_age.phen",header=T)
co=read.table("/gpfs/gpfs01/polaris/Q0286/uqaxue/phen/covar_sex_age_10PCs.txt",header=T)

ind=Reduce(intersect,list(phen$IID,co$IID))
phen=phen[match(ind,phen$IID),]

if (length(table(phen$trait[phen$trait>=0]))==2){
N=nrow(phen[phen$trait>=0,])
} else {
N=sum(!is.na(phen$trait))
}

print(paste("The sample size is ",N),sep="")
## How many columns of P-valeu, keep the later one if there are two
if(ncol(trait)==12){trait=trait[,-12]}
trait=trait[,-4]

##Flip
snbuf=which(trait$A1FREQ>0.5)
t=trait[snbuf,"ALLELE1"]
trait[snbuf,"ALLELE1"]=trait[snbuf,"ALLELE0"]
trait[snbuf,"ALLELE0"]=t
trait[snbuf,"BETA"]=(-1)*trait[snbuf,"BETA"]
trait[snbuf,"A1FREQ"]=1-trait[snbuf,"A1FREQ"]

###### Determine if it is case control phenotype 

if (length(table(phen$trait[phen$trait>=0]))==2){
## Rename the colnames

K=table(phen$trait[phen$trait>=0])[2]/N
print("trait is a binary phenotype")
print(paste("The prevalence is ",K),sep="")

colnames(trait)=c("SNP","Chr","bp","A1","A2","Freq","N","b","se","p")

trait$N=N*(1-trait$N)
sample=trait[,c(1,7)]
trait=trait[,-7]

trait$K=rep(K,nrow(trait))
trait=trait[,c(2,1,3:10)]

## Convert
source("/gpfs/gpfs01/polaris/Q0286/uqaxue/transfer_tmp/qbi/script/OR_transformation_Luke_new.R")
new=LmToOddsRatio(trait,K,0)

## Clean
new$N=as.integer(sample$N)
new$BETA=log(new$OR)
new$SE=new$BETA*new$se/new$b

new$b=new$BETA
new$se=new$SE

colnames(new)[1:9]=c("CHR","SNP","BP","A1","A2","freq","b","se","P")

trait=na.omit(new)
trait=trait[,c("CHR","SNP","BP","A1","A2","freq","N","b","se","P")]
print("OR transformation completed.")

} else {
print("trait is a quantitative phenotype.")
## Rename the colnames
colnames(trait)=c("SNP","CHR","BP","A1","A2","freq","N","b","se","P")

trait$N=N*(1-trait$N)
trait=trait[,c(2,1,3:10)]


## Clean
trait$N=as.integer(trait$N)


}

#### !!! Need to be careful if there is any NA or Inf in the data !!!! ####
# print(paste(sum("Any NA for trait: ",rowSums(is.na(trait))>0),sep=""))

## Freq check

snbuf=which(trait$freq>=0.5)
t=trait[snbuf,"A1"]
trait[snbuf,"A1"]=trait[snbuf,"A2"]
trait[snbuf,"A2"]=t
trait[snbuf,"b"]=(-1)*trait[snbuf,"b"]
trait[snbuf,"freq"]=1-trait[snbuf,"freq"]

##### Finish cleaning summary
print("Start to save summary statistics......")
## Save common SNPs with MAF >= 0.01
# trait$b=formatC(trait$b,digits=6,format="f")
# trait$se=formatC(trait$se,digits=6,format="f")

write.table(trait[trait$freq>=0.01,],paste(PATH,"ukbEUR_trait_common.txt",sep=""),row.names=F,col.names=T,quote=F)

## Also save rare summary
write.table(trait[trait$freq<0.01,],paste(PATH,"ukbEUR_trait_rare.txt",sep=""),row.names=F,col.names=T,quote=F)

## Save COJO & SMR format file

write.table(trait[trait$freq>=0.01,c("SNP","A1","A2","freq","b","se","P","N")],paste(PATH,"ukbEUR_trait_cojo.txt",sep=""),row.names=F,col.names=T,quote=F)
## Save LDSC format file
hm3=read.table("/shares/compbio/Group-Yang/a.xue/app/ldsc/eur_w_ld_chr/w_hm3.snplist",header=T)
trait$Z=trait$b/trait$se
snp=Reduce(intersect,list(hm3$SNP,trait$SNP))
new=trait[match(snp,trait$SNP),]

write.table(new[,c("SNP","A1","A2","Z","N","P")],paste(PATH,"trait_ldsc.txt",sep=""),row.names=F,col.names=T,quote=F)


## The End
print("trait QC completed!")




