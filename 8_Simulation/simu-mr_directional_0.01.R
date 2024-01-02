#
source("./comfun_file.R")
source("./MR.R")
library(MASS)

args = commandArgs(trailingOnly=TRUE)
taskid = as.integer(args[1])    # 1-135
### tunning parameters
# m_xpys = c("400-0-100", "360-40-100", "320-80-100","280-120-100", "200-200-100")
m_xpys = c("200-0-100", "180-20-100", "160-40-100","140-60-100", "100-100-100")
b_xys = c(0, 0.1, -0.1)
cor_r=c(-0.3,-0.1,-0.05,0,0.05,0.1,0.3)

m_xpy = as.character(expand.grid(m_xpys, b_xys,cor_r)[taskid, 1])
b_xy = as.numeric(expand.grid(m_xpys, b_xys,cor_r)[taskid, 2])
cor_r = as.numeric(expand.grid(m_xpys, b_xys,cor_r)[taskid, 3])
m_x = as.numeric(strsplit(m_xpy, "-")[[1]][1])
m_p = as.numeric(strsplit(m_xpy, "-")[[1]][2])
m_y = as.numeric(strsplit(m_xpy, "-")[[1]][3])
m = m_x + m_p + m_y

### fixed parameter
b_cx=0.3; b_cy=0.3; sigmasq_c=1; Rsq_zx=0.1 ; Rsq_zp=0.01
n1a=2e4; n1b=2e4; n1=n1a+n1b; n2=5e3; n3=5e3; n=n1+n2+n3
rep=300; res=c()

### output
outfolder = "Rsq_zp_0.01/results_directional/data/"
#outfolder = "results2/data/"
mydircreate(outfolder)
out = paste(outfolder, "simu3_", m_xpy, "_bxy_", b_xy,"_cor_",cor_r, ".txt", sep="")

### function
SimGenoValue = function(m, genostd, sigmasq_b){
  if(m==0){
    b=NULL; gvalue=0
  }else{
    b = rnorm(m) * sqrt(sigmasq_b); gvalue=genostd %*% b
  }
  return(list(b=b, gvalue=gvalue))
}

    
for(repindex in 1:rep){
  
set.seed(taskid*repindex*906)
### genotype
geno = simSNP(m, n, 0.01, 0.5)
z = geno$geno012
### genotypic value
genovalue_x = SimGenoValue(m_x, geno$genostd[,1:m_x], sigmasq_b = Rsq_zx/(m_x+m_p))
genovalue_y = SimGenoValue(m_y, geno$genostd[,(m_x+m_p+1):m], sigmasq_b = Rsq_zx/(m_x+m_p))
### pleiotropy value
# cor_r=0.8 # directional pleio
Sigma <- matrix(c(1,cor_r,cor_r,1),2,2)

  if(m_p==0){
    b_px=NULL; gvalue_px=0
    b_py=NULL; gvalue_py=0
  }else{
    mb=mvrnorm(n=m_p, rep(0, 2), Sigma,empirical = T)
    b_px = mb[,1] * sqrt(Rsq_zp/(m_x+m_p)); gvalue_px = geno$genostd[,(m_x+1):(m_x+m_p)] %*% b_px
    b_py = mb[,2] * sqrt(Rsq_zp/(m_x+m_p)); gvalue_py = geno$genostd[,(m_x+1):(m_x+m_p)] %*% b_py
  }

genovalue_px = list(b=b_px, gvalue=gvalue_px)
genovalue_py = list(b=b_py, gvalue=gvalue_py)


### confounder
c = rnorm(n)
### exposure phenotype (x)
e_x = rnorm(n) * as.numeric(sqrt(1-var(genovalue_x$gvalue + genovalue_px$gvalue + b_cx*c)))
x = genovalue_x$gvalue + genovalue_px$gvalue + b_cx*c + e_x
### outcome phenotype (y)
e_y = rnorm(n) * as.numeric(sqrt(1 - var(x*b_xy + genovalue_py$gvalue + genovalue_y$gvalue + b_cy*c)))
y = x*b_xy + genovalue_py$gvalue + genovalue_y$gvalue + b_cy*c + e_y
### association in training dataset
ss_zx_train = t(apply(z[1:n1a,], 2, function(a, x) as.numeric(coefficients(summary(lm(x~a)))[2,]), x=x[1:n1a]))
ss_zy_train = t(apply(z[(n1a+1):n1,], 2, function(a, y) as.numeric(coefficients(summary(lm(y~a)))[2,]), y=y[(n1a+1):n1]))
maf=colMeans(geno$geno012)/2
ss_zx_train_sig = ss_zx_train[ss_zx_train[,4]<5e-8,]
ss_zy_train_sig = ss_zy_train[ss_zx_train[,4]<5e-8,]
index=which(ss_zx_train[,4]<5e-8)
nr_snp=length(index)
nr_causal=sum(index %in% c(1:m_x))

if(nr_snp<5){
	print(paste0("rep ",repindex," fail. Nr of SNPs < 5."))
	rep=rep+1
	next}

for(method in c("GSMR_heidi_v1","GSMR_heidi_v3_stepwise", "GSMR_noheidi", "IVW", "Median", "Mode", "Egger", "Robust", "Lasso", "RAPS", "PRESSO", "MRMix", "Con-Mix")){

res_mr = all_mr(bzx=ss_zx_train_sig[,1], bzx_se=ss_zx_train_sig[,2], bzx_pval=ss_zx_train_sig[,4],bzy=ss_zy_train_sig[,1], bzy_se=ss_zy_train_sig[,2], bzy_pval=ss_zy_train_sig[,4], method=method)

if(is.null(res_mr$se)){res_mr$se=NA}
if(is.null(res_mr$pleio)){res_mr$pleio=NA}

rg=cor(ss_zx_train[,1],ss_zy_train[,1])

pleio = length(res_mr$pleio)
pleio_true = sum(as.numeric(gsub("SNP","",res_mr$pleio))>m_x)

res = rbind(res,data.frame(repindex, m_xpy, b_xy, method, p=res_mr$pval, bxy=res_mr$bxy, se=res_mr$se,rg=rg,nr_snp=nr_snp,nr_causal=nr_causal,nr_pleio_rm=pleio,nr_pleio_true=pleio_true))

# write.table(res, out, sep="\t", quote=F, row.names=F, col.names=T)
        }
print(paste("taskid:",taskid,", rep:",repindex,", SNP: ",m_xpy,", bxy: ",b_xy,
	"cor: ",cor_r,sep=" "))
}

write.table(res, out, sep="\t", quote=F, row.names=F, col.names=T)


## END






