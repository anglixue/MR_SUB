####
# Simulation to prove differences in effect sizes in exposure-stratified regressions or 
# even stratified rg estimates are direct indications of non-linear effects

source("comfun_file.R")
source("MR.R")
library(MASS)
library(survey)
library(dplyr)
library(ggplot2)

# args = commandArgs(trailingOnly=TRUE)
# taskid = as.integer(args[1])    # 1-135
taskid = 53
### tunning parameters
m_xpys = c("200-0-100", "180-20-100", "160-40-100","140-60-100", "100-100-100")
b_xys = c(0, 0.2, -0.2)
cor_r = c(-0.5,-0.3,-0.1,0,0.1,0.3,0.5)

m_xpy = as.character(expand.grid(m_xpys, b_xys,cor_r)[taskid, 1])
b_xy = as.numeric(expand.grid(m_xpys, b_xys,cor_r)[taskid, 2])
cor_r = as.numeric(expand.grid(m_xpys, b_xys,cor_r)[taskid, 3])
m_x = as.numeric(strsplit(m_xpy, "-")[[1]][1])
m_p = as.numeric(strsplit(m_xpy, "-")[[1]][2])
m_y = as.numeric(strsplit(m_xpy, "-")[[1]][3])
m = m_x + m_p + m_y

### fixed parameter
Rsq_zx = 0.3; Rsq_zp = 0.1
n1a = 1e5; n1b = 1e5; n1 = n1a + n1b; n = n1
# n2=5e3; n3=5e3; n=n1+n2+n3
rep = 300; res = c()

### function
SimGenoValue = function(m, genostd, sigmasq_b){
  if(m==0){
    b=NULL; gvalue=0
  }else{
    b = rnorm(m) * sqrt(sigmasq_b); gvalue=genostd %*% b
  }
  return(list(b=b, gvalue=gvalue))
}

# res1 for non-linear causal effect
res1 = data.frame(repindex = 1:100, cor_rb = NA, cor_g1 = NA, cor_g2 = NA)
# res2 for linear causal effect (choose the 5th quantile as turning point)
res2 = data.frame(repindex = 1:100, cor_rb = NA, cor_g1 = NA, cor_g2 = NA)

for (repindex in 1:100){
#repindex=40
set.seed(taskid*repindex*906)
### genotype
geno = simSNP(m, n, 0.01, 0.5)
z = geno$geno012
### genotypic value
genovalue_x = SimGenoValue(m_x, geno$genostd[,1:m_x], sigmasq_b = Rsq_zx/(m_x+m_p))
genovalue_y = SimGenoValue(m_y, geno$genostd[,(m_x+m_p+1):m], sigmasq_b = Rsq_zx/(m_x+m_p))
### pleiotropy value
genovalue_px = SimGenoValue(m_p, geno$genostd[,(m_x+1):(m_x+m_p)], sigmasq_b = Rsq_zp/(m_x+m_p))
genovalue_py = SimGenoValue(m_p, geno$genostd[,(m_x+1):(m_x+m_p)], sigmasq_b = Rsq_zp/(m_x+m_p))

### exposure phenotype (x)
e_x = rnorm(n) * as.numeric(sqrt(1-var(genovalue_x$gvalue + genovalue_px$gvalue)))
x = genovalue_x$gvalue + genovalue_px$gvalue + e_x
### outcome phenotype (y)
e_y = rnorm(n) * as.numeric(sqrt(1 - var(x*b_xy + genovalue_py$gvalue + genovalue_y$gvalue )))
y = (x^2+x)*b_xy + genovalue_py$gvalue + genovalue_y$gvalue + e_y
y = scale(y)[,1]
### association in training dataset
ss_zx_train = t(apply(z[1:n1a,], 2, function(a, x) as.numeric(coefficients(summary(lm(x~a)))[2,]), x=x[1:n1a]))
ss_zy_train = t(apply(z[(n1a+1):n1,], 2, function(a, y) as.numeric(coefficients(summary(lm(y~a)))[2,]), y=y[(n1a+1):n1]))
maf = colMeans(geno$geno012)/2
ss_zx_train_sig = ss_zx_train[ss_zx_train[,4] < 5e-8,]
ss_zy_train_sig = ss_zy_train[ss_zx_train[,4] < 5e-8,]
index = which(ss_zx_train[,4] < 5e-8)
nr_snp = length(index)
nr_causal = sum(index %in% c(1:m_x))

##
res = cbind(x,y)
res = as.data.frame(res)
colnames(res) = c("exposure","outcome")
rownames(res) = 1:nrow(res)

# Find the turning point based on 10 quantiles
res$group = NA
nr_group = 10
prob = quantile(res$exposure, probs = seq(0, 1, 1/nr_group))
res$group = cut(res$exposure, breaks = prob, labels = 1:nr_group, include.lowest = T, right = T)

est = res %>% group_by(group) %>% summarize(count = n(), x = mean(exposure), y = mean(outcome), y_sd = sd(outcome))
est = as.data.frame(est)
est$y_se = est$y_sd / sqrt(est$count)

turn_p = which.min(est$y)

# Divide the cohort into two subgroups
res$group2 = as.numeric(as.numeric(res$group) >= turn_p) + 1

## First evaluate if the SNP effects will change direction

# Generate the control group
res[res$exposure < prob[1], "group2"] = 0

# association in group 1
index1 = which(res$group2 == 1)

ss_zx_train1 = t(apply(z[index1[index1<n1a],], 2, function(a, x) as.numeric(coefficients(summary(lm(x~a)))[2,]), x = x[index1[index1<n1a]]))
ss_zy_train1 = t(apply(z[index1[index1>=n1a],], 2, function(a, y) as.numeric(coefficients(summary(lm(y~a)))[2,]), y = y[index1[index1>=n1a]]))
ss_zx_train_sig1 = ss_zx_train1[ss_zx_train1[,4] < 5e-8,]
ss_zy_train_sig1 = ss_zy_train1[ss_zx_train1[,4] < 5e-8,]

# association in group 2
index2 = which(res$group2 == 2)

ss_zx_train2 = t(apply(z[index2[index2<n1a],], 2, function(a, x) as.numeric(coefficients(summary(lm(x~a)))[2,]), x = x[index2[index2<n1a]]))
ss_zy_train2 = t(apply(z[index2[index2>=n1a],], 2, function(a, y) as.numeric(coefficients(summary(lm(y~a)))[2,]), y = y[index2[index2>=n1a]]))
ss_zx_train_sig2 = ss_zx_train2[ss_zx_train2[,4] < 5e-8,]
ss_zy_train_sig2 = ss_zy_train2[ss_zx_train2[,4] < 5e-8,]

## Second evaluate if the rg will be in the opposite direction
res1$cor_rb[repindex] = cor(ss_zx_train[,1], ss_zy_train[,1])
res1$cor_g1[repindex] = cor(ss_zx_train1[,1], ss_zy_train1[,1])
res1$cor_g2[repindex] = cor(ss_zx_train2[,1], ss_zy_train2[,1])

# res1 = rbind(res1,est)

print(repindex)
}

for (repindex in 1:100){
  #repindex=40
  set.seed(taskid*repindex*906)
  ### genotype
  geno = simSNP(m, n, 0.01, 0.5)
  z = geno$geno012
  ### genotypic value
  genovalue_x = SimGenoValue(m_x, geno$genostd[,1:m_x], sigmasq_b = Rsq_zx/(m_x+m_p))
  genovalue_y = SimGenoValue(m_y, geno$genostd[,(m_x+m_p+1):m], sigmasq_b = Rsq_zx/(m_x+m_p))
  ### pleiotropy value
  genovalue_px = SimGenoValue(m_p, geno$genostd[,(m_x+1):(m_x+m_p)], sigmasq_b = Rsq_zp/(m_x+m_p))
  genovalue_py = SimGenoValue(m_p, geno$genostd[,(m_x+1):(m_x+m_p)], sigmasq_b = Rsq_zp/(m_x+m_p))
  
  ### exposure phenotype (x)
  e_x = rnorm(n) * as.numeric(sqrt(1-var(genovalue_x$gvalue + genovalue_px$gvalue)))
  x = genovalue_x$gvalue + genovalue_px$gvalue + e_x
  ### outcome phenotype (y)
  e_y = rnorm(n) * as.numeric(sqrt(1 - var(x*b_xy + genovalue_py$gvalue + genovalue_y$gvalue )))
  # linear effect
  y = x * b_xy + genovalue_py$gvalue + genovalue_y$gvalue + e_y
  y = scale(y)[,1]
  ### association in training dataset
  ss_zx_train = t(apply(z[1:n1a,], 2, function(a, x) as.numeric(coefficients(summary(lm(x~a)))[2,]), x=x[1:n1a]))
  ss_zy_train = t(apply(z[(n1a+1):n1,], 2, function(a, y) as.numeric(coefficients(summary(lm(y~a)))[2,]), y=y[(n1a+1):n1]))
  maf = colMeans(geno$geno012)/2
  ss_zx_train_sig = ss_zx_train[ss_zx_train[,4] < 5e-8,]
  ss_zy_train_sig = ss_zy_train[ss_zx_train[,4] < 5e-8,]
  index = which(ss_zx_train[,4] < 5e-8)
  nr_snp = length(index)
  nr_causal = sum(index %in% c(1:m_x))
  
  ##
  res = cbind(x,y)
  res = as.data.frame(res)
  colnames(res) = c("exposure","outcome")
  rownames(res) = 1:nrow(res)
  
  # Find the turning point based on 10 quantiles
  res$group = NA
  nr_group = 10
  prob = quantile(res$exposure, probs = seq(0, 1, 1/nr_group))
  res$group = cut(res$exposure, breaks = prob, labels = 1:nr_group, include.lowest = T, right = T)
  
  est = res %>% group_by(group) %>% summarize(count = n(), x = mean(exposure), y = mean(outcome), y_sd = sd(outcome))
  est = as.data.frame(est)
  est$y_se = est$y_sd / sqrt(est$count)
  
  # set the turning point at 50% quantile
  turn_p = 5
  
  # Divide the cohort into two subgroups
  res$group2 = as.numeric(as.numeric(res$group) > turn_p) + 1
  
  ## First evaluate if the SNP effects will change direction
  
  # Generate the control group
  res[res$exposure < prob[1], "group2"] = 0
  
  # association in group 1
  index1 = which(res$group2 == 1)
  
  ss_zx_train1 = t(apply(z[index1[index1<n1a],], 2, function(a, x) as.numeric(coefficients(summary(lm(x~a)))[2,]), x = x[index1[index1<n1a]]))
  ss_zy_train1 = t(apply(z[index1[index1>=n1a],], 2, function(a, y) as.numeric(coefficients(summary(lm(y~a)))[2,]), y = y[index1[index1>=n1a]]))
  ss_zx_train_sig1 = ss_zx_train1[ss_zx_train1[,4] < 5e-8,]
  ss_zy_train_sig1 = ss_zy_train1[ss_zx_train1[,4] < 5e-8,]
  
  # association in group 2
  index2 = which(res$group2 == 2)
  
  ss_zx_train2 = t(apply(z[index2[index2<n1a],], 2, function(a, x) as.numeric(coefficients(summary(lm(x~a)))[2,]), x = x[index2[index2<n1a]]))
  ss_zy_train2 = t(apply(z[index2[index2>=n1a],], 2, function(a, y) as.numeric(coefficients(summary(lm(y~a)))[2,]), y = y[index2[index2>=n1a]]))
  ss_zx_train_sig2 = ss_zx_train2[ss_zx_train2[,4] < 5e-8,]
  ss_zy_train_sig2 = ss_zy_train2[ss_zx_train2[,4] < 5e-8,]
  
  ## Second evaluate if the rg will be in the opposite direction
  res2$cor_rb[repindex] = cor(ss_zx_train[,1], ss_zy_train[,1])
  res2$cor_g1[repindex] = cor(ss_zx_train1[,1], ss_zy_train1[,1])
  res2$cor_g2[repindex] = cor(ss_zx_train2[,1], ss_zy_train2[,1])
  
  print(repindex)
}
# Save results

write.table(res1, "nonlinear_causal_strat_rg.txt", row.names = F, col.names = F, quote = F)
write.table(res2, "linear_causal_strat_rg.txt", row.names = F, col.names = F, quote = F)

