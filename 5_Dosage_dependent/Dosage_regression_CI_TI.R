####
library(ggplot2)
library(dplyr)
detach(package:plyr)
######### plot A #########
phen=read.table("~/Downloads/addiction_07MAR2018/phen/Addiction_Substance_QCed.phen",header=T)
raw=read.table("~/Downloads/addiction_07MAR2018/phen/raw_TI_CI.phen",header=T)
ukb=read.table("~/Downloads/addiction_07MAR2018/regression/ukb_v2_18_diseases_match_gera.phen",header=T)
co=read.table("~/Downloads/addiction_07MAR2018/phen/covar_sex_age_10PCs_alcohol_change_meal_SES.txt",header=T)
bmi=read.table("~/Downloads/addiction_07MAR2018/GSMR/BMI/9280_12505_UKBiobank_T2DBMI_240717.tab",header=T)

phen=phen[match(co$IID,phen$IID),]
ukb=ukb[match(co$IID,ukb$IID),]
bmi=bmi[match(co$IID,bmi$IID),]

## Replace values > 10 as 10
phen[phen$TI<0,"TI"]=NA
phen[phen$CI<0,"CI"]=NA
phen=phen[,1:4]
ukb$CASES=ukb$CASES2
ukb$SUM_OF_CASES=ukb$SUM_OF_CASES2
ukb=ukb[,c(-23,-24)]

## Combine the data sets
phen=cbind(phen,ukb[,c(-1,-2)])
back=phen
phen[!is.na(phen$TI) & phen$TI>10,"TI"]=10
phen[!is.na(phen$CI) & phen$CI>10,"CI"]=10


## Rough check about the phenotypic correlation
res1=phen %>% group_by(TI) %>% summarize(count=n(),risk=mean(SUM_OF_CASES,na.rm=T))
res1=as.data.frame(res1)

res2=phen %>% group_by(CI) %>% summarize(count=n(),risk=mean(SUM_OF_CASES,na.rm=T))
res2=as.data.frame(res2)

plot(res1$TI,res1$risk)
plot(res2$CI,res2$risk)

#### Main analysis for CI ####

## Add covariates, |sex|age|10 PCs|
main=cbind(phen,co[,c(3:14)])
main=main[!is.na(main$CI),]

## Control
con=main[main$CI==0,]

## Case groups
case=list()
for (i in 1:10){
  case[[i]]=main[main$CI==i,]
}

## results file
res=data.frame(CI=NA,beta=NA,se=NA,P=NA,Disease=NA)
res[1,1]=mean(con$CI,na.rm=T)
res[1,2]=0
res[1,3]=0
res[1,4]=1
res[1,5]=NA

## Disease name
dis=c("Asthma", "Allergic rhinitis", "Cardiovascular disease","Cancer", "Depressive disorder", 
      "Diabetes type 2", "Dyslipidemia", "Hypertensive disease", "Hemorrhoids", 
      "Hernia abdominopelvic cavity","Iron deficiency anemias", "Irritable bowel syndrome", 
      "Osteoarthritis", "Osteoporosis", "Peripheral vascular disease","Peptic ulcers", "Psychiatric disorder", 
      "Varicose veins of lower extremities")
res_all=res[0,]

for (j in 1:18){
  ## main analysis
  for (i in 1:10){
    tmp_con=con
    tmp_case=case[[i]]
    tmp_con$CI_group=0
    tmp_case$CI_group=1
    
    tmp=rbind(tmp_con,tmp_case)
    tmp$Sex=as.factor(tmp$Sex)
    tmp$Age=scale(tmp$Age)
    
    res[i+1,1]=mean(case[[i]]$CI,na.rm=T)
    
    res[i+1,c(2,3,4)]=coef(summary(glm(tmp[,j+4]~tmp$CI_group+
                                         factor(tmp$Sex)+.,data=tmp[,c(26:36)]),
                                   family=binomial(link='logit')))[2,c(1,2,4)]
    res[i+1,5]=dis[j]
    print(paste(dis[j],"vs",i,"cups",sep=" "))
  }
  
  res_all=rbind(res_all,res)
  
}
## assian names for controls results
res_all[seq(1,188,by=11),5]=dis
res_all[,5]=rep(dis,each=11)

res_CI=res_all

raw <- ggplot(res_CI, aes(x=as.factor(CI), y=exp(beta))) + 
  geom_hline(yintercept=1, linetype=2) +
  geom_point(size=0.5) +
  geom_errorbar(aes(ymax=exp(beta+1.96*se), ymin=exp(beta-1.96*se)), cex=0.3,width=0.2) + 
  facet_wrap(vars(Disease), nrow = 6, scales = "free") +
  # , scales = "free" +
  theme(strip.text.x = element_text(size = 8)) +
  #ggtitle("Cardiovascular diseases") +
  xlab("Coffee Intake (cups/day)") +
  ylab("Odds ratio (OR)") 
raw

ggsave("~/Documents/Study/Lab Work/Project3_addiction/2.2 addiction GSMR/Sup Figure/Sup_Figure_9_Coffee_intake_dosage_06OCT2020.png",
       raw,width=7,height=7)

# Threshold when OR change direction for 18 diseases
thresh=c(5,6,NA,8,6,7,8,NA,NA,NA,NA,NA,8,7,NA,6,8,NA)


#### Main analysis for TI ####

## Add covariates, |sex|age|10 PCs|
main=cbind(phen,co[,c(3:14)])
main=main[!is.na(main$TI),]

## Control
con=main[main$TI==0,]

## Case groups
case=list()
for (i in 1:10){
  case[[i]]=main[main$TI==i,]
}

## results file
res=data.frame(TI=NA,beta=NA,se=NA,P=NA,Disease=NA)
res[1,1]=mean(con$TI,na.rm=T)
res[1,2]=0
res[1,3]=0
res[1,4]=1
res[1,5]=NA

## Disease name
dis=c("Asthma", "Allergic rhinitis", "Cardiovascular disease","Cancer", "Depressive disorder", 
      "Diabetes type 2", "Dyslipidemia", "Hypertensive disease", "Hemorrhoids", 
      "Hernia abdominopelvic cavity","Iron deficiency anemias", "Irritable bowel syndrome", 
      "Osteoarthritis", "Osteoporosis", "Peripheral vascular disease","Peptic ulcers", "Psychiatric disorder", 
      "Varicose veins of lower extremities")
res_all=res[0,]

for (j in 1:18){
## main analysis
for (i in 1:10){
  tmp_con=con
  tmp_case=case[[i]]
  tmp_con$TI_group=0
  tmp_case$TI_group=1
  
  tmp=rbind(tmp_con,tmp_case)
  tmp$Sex=as.factor(tmp$Sex)
  tmp$Age=scale(tmp$Age)
  
  res[i+1,1]=mean(case[[i]]$TI,na.rm=T)
  
  res[i+1,c(2,3,4)]=coef(summary(glm(tmp[,j+4]~tmp$TI_group+
                                       factor(tmp$Sex)+.,data=tmp[,c(26:36)]),
                                 family=binomial(link='logit')))[2,c(1,2,4)]
  res[i+1,5]=dis[j]
  print(paste(dis[j],"vs",i,"cups",sep=" "))
}

  res_all=rbind(res_all,res)
  
}
## assian names for controls results
res_all[seq(1,188,by=11),5]=dis
res_all[,5]=rep(dis,each=11)

res_TI=res_all

raw <- ggplot(res_TI, aes(x=as.factor(TI), y=exp(beta))) + 
  geom_hline(yintercept=1, linetype=2) +
  geom_point(size=0.5) +
  geom_errorbar(aes(ymax=exp(beta+1.96*se), ymin=exp(beta-1.96*se)), cex=0.3,width=0.2) + 
  facet_wrap(vars(Disease), nrow = 6, scales = "free") +
  # , scales = "free" +
  theme(strip.text.x = element_text(size = 8)) +
  #ggtitle("Cardiovascular diseases") +
  xlab("Tea Intake (cups/day)") +
  ylab("Odds ratio (OR)") 
raw

ggsave("~/Documents/Study/Lab Work/Project3_addiction/2.2 addiction GSMR/Sup Figure/Sup_Figure_10_Tea_intake_dosage_06OCT2020.png",
       raw,width=7,height=7)

# Threshold when OR change direction for 18 diseases
thresh=c(NA,8,8,8,6,NA,NA,6,5,8,4,6,6,NA,6,NA,6,3)


#### CPD ####
####
library(ggplot2)
library(dplyr)
detach(package:plyr)
######### plot A #########
phen=read.table("~/Downloads/addiction_07MAR2018/phen/Addiction_Substance_QCed.phen",header=T)
raw=read.table("~/Downloads/addiction_07MAR2018/phen/raw_TI_CI.phen",header=T)
ukb=read.table("~/Downloads/addiction_07MAR2018/regression/ukb_v2_18_diseases_match_gera.phen",header=T)
co=read.table("~/Downloads/addiction_07MAR2018/phen/covar_sex_age_10PCs_alcohol_change_meal_SES.txt",header=T)
bmi=read.table("~/Downloads/addiction_07MAR2018/GSMR/BMI/9280_12505_UKBiobank_T2DBMI_240717.tab",header=T)
batch=read.table("~/Downloads/addiction_07MAR2018/phen/24101_12505_Genotyping_Batch.tab",header=T)

phen=phen[match(co$IID,phen$IID),]
ukb=ukb[match(co$IID,ukb$IID),]
bmi=bmi[match(co$IID,bmi$IID),]
batch=batch[match(co$IID,batch$f.eid),]

## Replace values > 10 as 10
phen[!is.na(phen$CPD) & phen$CPD<0,"CPD"]=NA
phen=phen[,c(1,2,6,9)]
ukb$CASES=ukb$CASES2
ukb$SUM_OF_CASES=ukb$SUM_OF_CASES2
ukb=ukb[,c(-23,-24)]

## Combine the data sets
phen=cbind(phen,ukb[,c(-1,-2)])
back=phen

## Rough check about the phenotypic correlation
res1=phen %>% group_by(CPD) %>% summarize(count=n(),risk=mean(SUM_OF_CASES,na.rm=T))
res1=as.data.frame(res1)

plot(res1$CPD,res1$risk)

#### Main analysis for CI ####

## Add covariates, |sex|age|10 PCs|
main=cbind(phen,co[,c(3:14)])

## Control
con=main[!is.na(main$SI) & main$SI==0,]
con=con[,-3]
main=main[!is.na(main$CPD),]
main=main[,-3]
## Case groups
case=list()
for (i in 1:6){
  case[[i]]=main[main$CPD>=(5*i-4) & main$CPD<=(5*i),]
}
case[[7]]=main[main$CPD>30,]

## results file
res=data.frame(CPD=NA,beta=NA,se=NA,P=NA,Disease=NA)
res[1,1]=0
res[1,2]=0
res[1,3]=0
res[1,4]=1
res[1,5]=NA

## Disease name
dis=c("Asthma", "Allergic rhinitis", "Cardiovascular disease","Cancer", "Depressive diorder", 
      "Diabetes type 2", "Dyslipidemia", "Hypertensive disease", "Hemorrhoids", 
      "Hernia abdominopelvic cavity","Iron deficiency anemias", "Irritable bowel syndrome", 
      "Osteoarthritis", "Osteoporosis", "Peripheral vascular disease","Peptic ulcers", "Psychiatric disorder", 
      "Varicose veins of lower extremities")
res_all=res[0,]

for (j in 1:18){
  ## main analysis
  for (i in 1:7){
    tmp_con=con
    tmp_case=case[[i]]
    tmp_con$CPD_group=0
    tmp_case$CPD_group=1
    
    tmp=rbind(tmp_con,tmp_case)
    tmp$Sex=as.factor(tmp$Sex)
    tmp$Age=scale(tmp$Age)
    
    res[i+1,1]=mean(case[[i]]$CPD,na.rm=T)
    
    res[i+1,c(2,3,4)]=coef(summary(glm(tmp[,j+4]~tmp$CPD_group+
                                         factor(tmp$Sex)+.,data=tmp[,c(26:36)]),
                                   family=binomial(link='logit')))[2,c(1,2,4)]
    res[i+1,5]=dis[j]
    print(paste(dis[j],"vs",i,"cigarettes",sep=" "))
  }
  
  res_all=rbind(res_all,res)
  
}
## assian names for controls results
res_all[seq(1,nrow(res_all),by=8),5]=dis
res_all[,5]=rep(dis,each=8)
res_all$CPD_group=rep(c("0","1~5","6~10","11~15","16~20","21~25","26~30",">30"),18)

res_all$CPD_group=as.factor(res_all$CPD_group)
res_all$CPD_group = factor(res_all$CPD_group, levels=levels(res_all$CPD_group)[c(2,3,8,4,5,6,7,1)])

res_CPD=res_all

raw <- ggplot(res_CPD, aes(x=as.factor(CPD_group), y=exp(beta))) + 
  geom_hline(yintercept=1, linetype=2) +
  geom_point(size=0.5) +
  geom_errorbar(aes(ymax=exp(beta+1.96*se), ymin=exp(beta-1.96*se)), cex=0.3,width=0.2) + 
  facet_wrap(vars(Disease), nrow = 6, scales = "free") +
  # , scales = "free" +
  theme(strip.text.x = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=6)) +
  xlab("Cigarettes per day (CPD)") +
  ylab("Odds ratio (OR)") 
raw

ggsave("~/Documents/Study/Lab Work/Project3_addiction/2.2 addiction GSMR/Sup Figure/Sup_Figure_3.5_CPD_dosage.png",
       raw,width=7,height=7)



####

phen$batch=batch$f.22000.0.0
phen=phen[!is.na(phen$batch),]
sum = phen %>% group_by(batch) %>% dplyr::summarize(count=n(),si=mean(SI,na.rm=T),cpd=mean(CPD,na.rm=T),
                                              asth=mean(ASTHMA,na.rm=T),aller=mean(ALLERGIC_RHINITIS,na.rm=T))

sum=as.data.frame(sum)


plot(sum$batch,sum$si)
plot(sum$batch,sum$cpd)
plot(sum$batch,sum$asth)
plot(sum$batch,sum$aller)










