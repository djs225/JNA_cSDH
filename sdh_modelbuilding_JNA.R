library(tidyverse)
library(rlist)
library(lattice)
library(survminer)
library(survival)
library(stargazer)
library(mice)
library(survival)
library(naniar)
library(gtools)
library(caret)
setwd("/Volumes/Untitled 1/Backup")
original<-read.csv("/Volumes/Untitled 1/Backup/sdh_final.csv")
RNGversion("3.5.3")#Allow reproducibility of imputation performed under previous version of R
set.seed(123)

 # to ensure reproducibility of imputation

###Variable tidying, create two data frames containing baseline patient variables (e.g. sex) 
# and then clinical data pertaining to the stay (e.g. complications on day 1 = any1)
orig_filt<-original %>% filter(!is.na(volatile)) %>% filter(op_time>0)
head(orig_filt)
orig_filt<-orig_filt %>% select(-height, -bmi)
orig_filt <-orig_filt %>% rename(mot=Initial.GCS..M) # Bit of renaming etc...
orig_filt<-orig_filt %>% rename(GCS=Initial.GCS) %>% rename(Sex=Gender)
orig_filt<-orig_filt %>% rename(wt=weight)
orig_filt <-orig_filt %>% rename(cog=cognitive_inpair)
orig_filt<-orig_filt %>% rename(reop=reop.y)
orig_filt<-orig_filt %>% mutate(asa=coalesce(asa, ASA)) # as available from two sources
orig_filt<-orig_filt %>% mutate(dos_comp=ifelse(is.na(sicker),tot0,tot0)) %>% mutate(sicker=ifelse(is.na(sicker),0,sicker))

baseline<-orig_filt %>% select(asa, dgh, Sex, GCS, mot, mRS, age, wt,cog,los, last_creat)
baseline2<-orig_filt %>% select(id, asa, dgh, Sex, GCS, mot, mRS, age, wt,cog, last_creat,reop) #need this as will need id to re-hash it back together once its imputed.
clinical<-orig_filt %>% select(id,dgh, any0,any1,any2,any3,any4,any5,any6,any7,any8,any9,any10,any11,any12,any13,tot0,tot1,tot2,tot3,tot4,tot5,tot6,tot7,tot8,tot9,tot10,tot11,tot12,tot13, MI, AKI, AKI_SENS, op_time,mactp,hypo_t,tcitp,co2tp,co2tn,fentanyl, morphine_iv,fentanyl_iv,alfentanil, morphine, volatile, tiva, airways, angina, anticoag, antiplatelet,heart_failure, ihd_htn,renal,smok_cess,los,low,lot,sicker,admission_comp,dos_comp,reop,lot)
md.pattern(baseline, rotate.names=TRUE)
orig_filt %>% filter(is.na(los))



#Recode predictor variables (e.g. group cvs disease, total time out of CO2 range)
clinical<-clinical %>% mutate(cvs_dx=ifelse(angina==1,1,ifelse(ihd_htn==1,1,0)))
clinical<-clinical %>% mutate(any_acoag=ifelse(anticoag==1,1,ifelse(antiplatelet==1,1,0)))
clinical<-clinical %>% mutate(co2_t=co2tn+co2tp)

#Morphine doses converted into fentanyl equivalents using equianalgesic doses from here: https://journals.sagepub.com/doi/pdf/10.1345/aph.1H421.  
#This gives a ratio of 50-100:1 for Morphine to Fentanyl.  Have chosen lower bound for this work.  This then presents the summary dose in mcg of fentanyl equivalents.  
# Rationale for this is that most commonly used opioid at our centre is fentanyl.  

head(clinical)
clinical_zero<-clinical %>%
  mutate_all(funs(ifelse(is.na(.), 0, .)))
clinical_zero<-clinical_zero%>%mutate(fent_equiv=fentanyl+fentanyl_iv+(morphine_iv/1000*50))
#intraop<-cbind(intraop,clinical_zero$fent_equiv)
#intraop<-intraop %>% rename(fent_equiv=`clinical_zero$fent_equiv`)

##Create outcome variable for end organ complications.
# MI = Trop > ULN 
# AKI = Creat rise at any point to >1.5x baseline
clinical_zero<-clinical_zero %>% mutate(end_organcomp=ifelse(MI==1,1,ifelse(AKI==1,1,0)))

##Now combine baseline data witih `clinical zero frame above to give complete case data framer`
#comp_case becomes main dataframe for complete case analyses and as base for MI, then requires tidying/change of format of some variables (e.g. characters)
comp_case<-baseline2 %>% left_join(clinical_zero, by="id")
comp_case<-comp_case %>% rename(reop=reop.y)
comp_case$cog<-as.character(comp_case$cog)
comp_case<-comp_case %>% mutate(cog_cat=(ifelse(str_detect(cog,"Yes, longstanding|Yes, temporary|UTA"),"Yes",cog))) 
comp_case<-comp_case %>% replace_with_na((replace=list(cog_cat="")))
comp_case<-comp_case %>% mutate(cog_cat=ifelse(cog_cat=="Yes",1,cog_cat)) %>% mutate(cog_cat=ifelse(cog_cat=="No",0,cog_cat))
comp_case$cog_cat<-as.integer(comp_case$cog_cat)
comp_case<-comp_case %>% mutate(m6=ifelse(mot>5,1,0))
comp_case<-comp_case %>% mutate(gcs15=ifelse(GCS>14,1,0))
comp_case<-comp_case %>% mutate(status=1)##Same imputed dataset was used for cox regression, therefore this had to be included
cc_km<-with(comp_case,Surv(los, status==1))


##########MULTIPLE IMPUTATION################
# This is to allow univariable screening and multivariable modelling.  Validation of the selected model
#from this process (using stopping rules of p<0.05 on pooled LRT) is then validated using a separate
# fold and impute strategy in file sdh_cv_loops_JNA.R.

for_mi<-comp_case %>% select(asa,dgh.x,Sex,age,m6,gcs15,cog_cat,mRS,last_creat,cvs_dx,any_acoag,admission_comp,dos_comp,sicker,airways,los,low,op_time,volatile,fent_equiv,hypo_t,co2_t,any3,status,end_organcomp, heart_failure,reop)
for_mi<-for_mi %>% mutate(nelson=nelsonaalen(for_mi,los,status)) # Need this in the imputation model to allow cox regression
for_mi<-for_mi %>% mutate(los_quant=quantcut(los,q=4,na.rm=FALSE))
for_mi<-for_mi %>% mutate(prol_los=ifelse(los_quant=="(10.6,161]",1,0))
for_mi <-for_mi %>% mutate(miss_mrs=ifelse(is.na(mRS),1,0)) %>% mutate(miss_asa=ifelse(is.na(asa),1,0)) %>% mutate(miss_creat=ifelse(is.na(last_creat),1,0)) %>% mutate(miss_cog=ifelse(is.na(cog_cat),1,0))
for_risk<-for_mi #Clones the dataframe that will be used for imputation to create a startingi point for the validation process in other file
vis_miss(for_mi)
gg_miss_upset(for_mi)
gg_miss_var(for_mi)

## Imputation model building

for_mi$asa<-as.factor(for_mi$asa)
for_mi$mRS<-as.factor(for_mi$mRS)
for_mi$cog_cat<-as.factor(for_mi$cog_cat)

head(for_mi)
for_mi<-for_mi %>% select(-los_quant)
inlist<-c("nelson", "status", "prol_los", "end_organcomp")#specify that outcomes must be inclued in imputation model
outlist<-c("miss_mrs", "miss_asa", "miss_cog", "miss_creat")
missing_var_name<-c("asa", "cog_cat", "mRS", "last_creat")
for_mi<-for_mi[,c(1,7,8,9,2,3,4,5,6,10:33)]
pred<-quickpred(for_mi, minpuc=0.5, include=inlist, exclude=outlist)
head(for_mi)
dim(for_mi)

methods<-c("polr","logreg","polr","pmm","","","","","","","","","","","","","","","","","","","","","","","","","","","","","")
sdh_mi<-mice(for_mi, method=methods, m=40, maxit=10, pred=pred) #WE have 37 percent missing data for mRS.  THerefore do m = 40 (m>/=miss%)
summary(sdh_mi)
#densityplot(sdh_mi)
#stripplot(sdh_mi)


impL<-complete(sdh_mi, "long", include = T)
head(impL)
i2<-impL %>% filter(.imp!=0)
table(i2$asa,i2$miss_asa)
table(i2$cog_cat,i2$miss_cog)
table(i2$mRS,i2$miss_mrs)

imp.labs<-factor(c(0, 1))
names(imp.labs)<-c("Observed", "Imputed")

asa<-ggplot(data=i2,aes(x=asa,group=miss_asa))+geom_bar(aes(y=..prop..,fill=asa), stat="count", position="dodge")+facet_grid(~miss_asa, labeller=labeller(miss_asa=imp.labs))+geom_text(aes(label=scales::percent(..prop..),y=..prop..),stat="count",vjust=-.5) +labs (y="Percent", fill="ASA")+scale_y_continuous(labels=scales::percent)
asa

cog<-ggplot(data=i2,aes(x=cog_cat,group=miss_cog))+geom_bar(aes(y=..prop..,fill=cog_cat), stat="count", position="dodge")+facet_grid(~miss_cog,labeller=labeller(miss_cog=imp.labs))+geom_text(aes(label=scales::percent(..prop..),y=..prop..),stat="count",vjust=-.5) +labs (y="Percent", fill="Cog Concern")+scale_y_continuous(labels=scales::percent)
cog

mrsr<-ggplot(data=i2,aes(x=mRS,group=miss_mrs))+geom_bar(aes(y=..prop..,fill=mRS), stat="count", position="dodge")+facet_grid(~miss_mrs, labeller=labeller(miss_mrs=imp.labs))+geom_text(aes(label=scales::percent(..prop..),y=..prop..),stat="count",vjust=-.5) +labs (y="Percent", fill="mRS")+scale_y_continuous(labels=scales::percent)
mrsr
densityplot(sdh_mi, xlab="Creatinine in micromols/litre")

head(sdh_mi)
stripplot(sdh_mi,asa+cog_cat+mRS~.imp, col=c("grey", mdc(2),pch=c(1,20)), hor=20, jitter=TRUE)
plot(sdh_mi)
head(sdh_mi$predictorMatrix, 4)
stargazer(head(sdh_mi$predictorMatrix, 4)[,1:13])
stargazer(head(sdh_mi$predictorMatrix, 4)[,14:28])

##Univariable screening in MI dataset first of all, output to LaTeX table via Stargazer package
univ_vars<-c("age", "Sex", "asa", "mRS", "dgh.x", "cog_cat", "m6", "gcs15", "last_creat", "cvs_dx", "any_acoag", "airways", "heart_failure", "admission_comp", "dos_comp", "sicker", "low", "op_time", "volatile", "fent_equiv", "hypo_t", "co2_t", "any3", "end_organcomp")
mi_modelsend<-list()
pooled_end<-list()
univ_vars_comp<-univ_vars[1:22]
univ_vars_comp<-c(univ_vars_comp)
for (var in univ_vars_comp){
  mi_modelsend[[var]]<-with(data=sdh_mi, glm(end_organcomp~get(var),family=binomial(link=logit)))
  pooled_end[[var]]<-summary(pool(mi_modelsend[[var]]))
}
pooled_flat_end<-list.rbind(pooled_end)
pooled_flat_end

pooled_flat_end<-pooled_flat_end %>% mutate(var=row.names(pooled_flat_end))
head(pooled_flat_end)
pooled_flat_end<-pooled_flat_end %>% mutate((or=exp(estimate)))
pooled_flat_end<-pooled_flat_end %>% mutate((lower=exp(estimate-1.96*std.error)))
pooled_flat_end<-pooled_flat_end %>% mutate((upper=exp(estimate+1.96*std.error)))
pooled_flat_end<-pooled_flat_end %>% select(-df,-estimate)
pooled_flat_end<-pooled_flat_end[,c(5,6,7,8,2,3, 4)]
head(pooled_flat_end)


colnames(pooled_flat_end) <- (c("Variable","OR", "95% Lower", "95% Upper", "se", "z", "p"))
head(pooled_flat_end)
pooled_flat_end$Variable
pooled_flat_end<-pooled_flat_end %>% filter(Variable!="los") %>% filter(Variable!="status") %>% filter(Variable!="nelson")
pooled_flat_end$Variable
pooled_flat_end<-pooled_flat_end[-c(1,3,5,10,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55),]
pooled_flat_end<-pooled_flat_end[,c(2,3,4,7)]
mi_variables_stargazer<-c("Age","Male", "ASA2", "ASA3", "ASA4", 
                          "ASA5", "mRS 1", "mRS 2", "mRS 3", "mRS 4" ,"Non CUH patient", "Cognitive Impairment", 
                          "M6 on admission", "GCS 15 on admission",  "Creatinine", "History of CVS Disease", "Anticoagulant on admission",
                          "Airways Disease",  "Heart Failure", "Admission EPOMS","DOS EPOMS", "Pre-op Deterioration", "Length of Wait","Op Time", "Volatile anaesthetic", 
                          "Fentanyl dose", "Time MAP not 80mmHg", "Time CO2 not 3 5kPa")
rownames(pooled_flat_end)<-mi_variables_stargazer
stargazer(pooled_flat_end, summary=FALSE)
head(pooled_flat_end)
### Multivariable Model Building in Multiply Imputed Datasets
##So select those with p<0.2 and build/test across MI datasets

pre_1_mi<-with(sdh_mi,glm(end_organcomp~Sex+as.numeric(asa)+dgh.x+cog_cat+m6+gcs15+last_creat+cvs_dx+any_acoag+airways+heart_failure+admission_comp, family=binomial))#Remove m6
sum1<-summary(pool(pre_1_mi), conf.int=TRUE, exponentiate = TRUE)#REMOVE M6
sum1

pre_2_mi<-with(sdh_mi,glm(end_organcomp~Sex+as.numeric(asa)+dgh.x+cog_cat+gcs15+last_creat+cvs_dx+any_acoag+airways+heart_failure+admission_comp, family=binomial)) #remove ccf
sum2<-summary(pool(pre_2_mi), conf.int=TRUE, exponentiate = TRUE)# REMOVE HEART FAILURE
sum2[-1]<-round(sum2[-1],5)
sum2
pool.compare(pre_1_mi,pre_2_mi,method="likelihood")

pre_3_mi<-with(sdh_mi,glm(end_organcomp~Sex+as.numeric(asa)+dgh.x+cog_cat+gcs15+last_creat+cvs_dx+any_acoag+airways+admission_comp, family=binomial)) 
sum3<-summary(pool(pre_3_mi), conf.int=TRUE, exponentiate = TRUE)# REMOVE GCS15
sum3[-1]<-round(sum3[-1],5)
sum3
pool.compare(pre_2_mi,pre_3_mi,method="likelihood")

pre_4_mi<-with(sdh_mi,glm(end_organcomp~Sex+as.numeric(asa)+dgh.x+cog_cat+airways+last_creat+cvs_dx+any_acoag+admission_comp, family=binomial)) 
sum4<-summary(pool(pre_4_mi), conf.int=TRUE, exponentiate = TRUE)# REMOVE AIRWAYS
sum4[-1]<-round(sum4[-1],5)
sum4
pool.compare(pre_3_mi,pre_4_mi,method="likelihood")

pre_5_mi<-with(sdh_mi,glm(end_organcomp~Sex+as.numeric(asa)+dgh.x+cog_cat+last_creat+cvs_dx+any_acoag+admission_comp, family=binomial)) 
sum5<-summary(pool(pre_5_mi), conf.int=TRUE, exponentiate = TRUE)# REMOVE Male
sum5[-1]<-round(sum5[-1],5)
sum5
pool.compare(pre_4_mi,pre_5_mi,method="likelihood")

pre_6_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+cog_cat+last_creat+cvs_dx+any_acoag+admission_comp, family=binomial))
sum6<-summary(pool(pre_6_mi), conf.int=TRUE, exponentiate = TRUE)#remove creat
sum6[-1]<-round(sum6[-1],5)
sum6
pool.compare(pre_5_mi,pre_6_mi,method="likelihood")

pre_7_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+cog_cat+cvs_dx+any_acoag+admission_comp, family=binomial)) 
sum7<-summary(pool(pre_7_mi), conf.int=TRUE, exponentiate = TRUE)#remove cog cat
sum7[-1]<-round(sum7[-1],5)
sum7
pool.compare(pre_6_mi,pre_7_mi,method="likelihood")

pre_8_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+cvs_dx+any_acoag+admission_comp, family=binomial))
sum8<-summary(pool(pre_8_mi), conf.int=TRUE, exponentiate = TRUE)#remove cvs dx
sum8[-1]<-round(sum8[-1],5)
sum8
pool.compare(pre_7_mi,pre_8_mi,method="likelihood")

pre_9_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+any_acoag+admission_comp, family=binomial))
sum9<-summary(pool(pre_9_mi), conf.int=TRUE, exponentiate = TRUE)
sum9[-1]<-round(sum9[-1],5)
sum9
pool.compare(pre_8_mi,pre_9_mi,method="likelihood")

#Now add in the postop/DOS variables

post_1_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+any_acoag+admission_comp+dos_comp+op_time+fent_equiv+hypo_t+co2_t+low, family=binomial))
sumpost1<-summary(pool(post_1_mi), conf.int=TRUE, exponentiate = TRUE)
sumpost1[-1]<-round(sumpost1[-1],5)#remove optime
sumpost1



post_2_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+any_acoag+admission_comp+dos_comp+fent_equiv+hypo_t+co2_t+low, family=binomial))
sumpost2<-summary(pool(post_2_mi), conf.int=TRUE, exponentiate = TRUE)
sumpost2[-1]<-round(sumpost2[-1],5)#remove hypotime
sumpost2
pool.compare(post_1_mi,post_2_mi,method="likelihood")

post_3_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+any_acoag+admission_comp+dos_comp+fent_equiv+co2_t+low, family=binomial))
sumpost3<-summary(pool(post_3_mi), conf.int=TRUE, exponentiate = TRUE)
sumpost3[-1]<-round(sumpost3[-1],5)#remove low
sumpost3
pool.compare(post_2_mi,post_3_mi,method="likelihood")

post_4_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+any_acoag+admission_comp+dos_comp+fent_equiv+co2_t, family=binomial))
sumpost4<-summary(pool(post_4_mi), conf.int=TRUE, exponentiate = TRUE)
sumpost4[-1]<-round(sumpost4[-1],5)
sumpost4#remove admission epoms
pool.compare(post_3_mi,post_4_mi,method="likelihood")

post_5_mi<-with(sdh_mi,glm(end_organcomp~as.numeric(asa)+dgh.x+any_acoag+dos_comp+fent_equiv+co2_t, family=binomial))
sumpost5<-summary(pool(post_5_mi), conf.int=TRUE, exponentiate = TRUE)
sumpost5[-1]<-round(sumpost5[-1],5)
sumpost5
pool.compare(post_4_mi,post_5_mi,method="likelihood")

######### COMPLETE CASE ANALYSIS ###############

univ_vars<-c("age", "Sex", "asa", "mRS", "dgh.x", "cog_cat", "m6", "gcs15", "last_creat", "cvs_dx", "any_acoag", "airways", "heart_failure", "admission_comp", "dos_comp", "sicker", "low", "op_time", "volatile", "fent_equiv", "hypo_t", "co2_t")

univ_formulas_endorgan<-sapply(univ_vars,function(x)as.formula(paste('end_organcomp~',x)))
univ_models_endorgan<-lapply(univ_formulas_endorgan,function(x){glm(x, data=for_mi,family=binomial)})
univ_models_endorgan
univ_results_endorgan <- lapply(univ_models_endorgan,function(x){return(exp(cbind(coef(x),confint(x))))})
univ_results_endorgan

res_uni_cc_log<-list.rbind(univ_results_endorgan)
se<-vector()
z<-vector()
p_val<-vector()
for (i in univ_models_endorgan){
  j<-summary(i)
  coef<-j$coefficients
  se<-c(se,coef[,2])
  z<-c(z,coef[,3])
  p_val<-c(p_val,coef[,4])
}

res_uni_cc_log<-cbind(res_uni_cc_log
                      ,se,z,p_val)

res_uni_cc_log ## For complete case data

res_uni_cc_log <-res_uni_cc_log[-c(1,3,5,10,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53),]
res_uni_cc_log<-res_uni_cc_log[,-c(4,5)]

rownames(res_uni_cc_log)<-mi_variables_stargazer
res_uni_cc_log
stargazer(res_uni_cc_log)


