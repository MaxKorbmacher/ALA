# Data prep for:
# Replication of Serum Alpha-Linolenic Acid and Long-Term Multiple Sclerosis Activity and Progression
#
# Max Korbmacher, 21 Jul 2025
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
#
print("Load packages")
pacman::p_load(haven,dplyr,reshape2)
#
print("Load data.")
#
print("All data will be merged into the (long) format of df.")
d = read.csv("/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv",sep=";")
df = read.csv("/Users/max/Documents/Local/MS/data/final_dat_subc.csv")
demo = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/Demographics.sav") # Demographics
ALA = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/fa.sas7bdat")
relapse = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/relapsereport.sas7bdat") # relapse info
demo10 = read.csv('/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv',sep = ";") # 10 years follow up demo & clin scores
#
print("prep and select")
df$TotalVol = df$TotalGrayVol + df$CerebralWhiteMatterVol + df$Right.Cerebellum.White.Matter + df$Left.Cerebellum.White.Matter
df$TotalWMVol = df$CerebralWhiteMatterVol + df$Right.Cerebellum.White.Matter + df$Left.Cerebellum.White.Matter
df = df %>% select(eid, session, sex, age, lesion_count, EstimatedTotalIntraCranialVol,TotalVol,TotalGrayVol,TotalWMVol)
#
#
d$eid = d$Patnr
#sum_T1_lesions = d %>% select(eid, MRI_T1Gd_sum)

#
# PASAT in long format
PASAT = d %>% select(eid, Treatment_OFAMS, BL_PASATcorrect,M12_PASATcorrect,M24_PASATcorrect, smoking_OFAMS)
PASAT = melt(PASAT,id.vars = c("eid","Treatment_OFAMS","smoking_OFAMS"))
PASAT = rename(PASAT, PASAT = value, session = variable)
PASAT$session = factor(PASAT$session)
levels(PASAT$session) = c(0,12,24)

# T1 lesions
T1Gd = d %>% select(eid,MRI_T1Gd_0 , MRI_T1Gd_12,MRI_T1Gd_24)
T1Gd = melt(T1Gd,id.vars = "eid")
T1Gd = rename(T1Gd, new_T1Gd_lesion = value, session = variable)
T1Gd$session = factor(T1Gd$session)
levels(T1Gd$session) = c(0,12,24)

# ALA
ALA = ALA %>% select(patno, visitno, FA18X3N3)
names(ALA) = c("eid","session","ALA")

# T2_lesions = read.csv("/Users/max/Documents/Local/MS/data/lesion_count.csv")
#
# Relapses
#relapsenb = relapse %>% filter(RELAPSENEW==1) %>% select(patno, relapseno) %>% group_by(patno) %>% summarise(max_relapseno = max(relapseno))
relapse$eid = relapse$patno
demo10$eid = demo10$Patnr
relapse = merge(relapse,demo10,by="eid",all=T)
# Convert to Date objects, specifying format for BL_DATEOFVISIT
DATEOFONSET_dt <- as.Date(relapse$DATEOFONSET)
BL_DATEOFVISIT_dt <- as.Date(relapse$BL_DATEOFVISIT, format="%m/%d/%Y")
# Calculate the difference in days
relapse$time_difference_days_r <- as.numeric(DATEOFONSET_dt - BL_DATEOFVISIT_dt, units = "days")
# Calculate the difference in years
relapse$time_difference_years_r <- as.numeric(DATEOFONSET_dt - BL_DATEOFVISIT_dt, units = "days")/365

relapse$session = ifelse(relapse$time_difference_years_r<1,1,0)+
  ifelse(relapse$time_difference_years_r>1 & relapse$time_difference_years_r<2,13,0)+
  ifelse(relapse$time_difference_years_r>2,25,0)-1
relapse$session = factor(relapse$session)
relapse$RELAPSENEW = ifelse(is.na(relapse$RELAPSENEW) == T, 0,relapse$RELAPSENEW)
relapse = relapse %>% select(eid, session, RELAPSENEW)

# EDSS
edss = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/edss.sav")
edss = melt(edss,id.vars = c("Patnr"))
names(edss) = c("eid","session","edss")
edss$session = factor(edss$session)
levels(edss$session) = c(0,6,12,18,24)

#brain_vol = read.csv("/Users/max/Documents/Local/MS/data/icv.csv")

#data which are not used
#MRI = read_spss("/Users/max/Documents/Local/MS/demographics/Jurate-statistikkfiler/MRI_A_E_D.sav")
#interrim = read.csv("/Users/max/Documents/Local/MS/results/interrim_data.csv")

# correspondence between T1w and T2w lesion count
#cor(d$MRI_T1Gd_12,d$MRI_T2_12,use = "complete.obs")
#cor(d$MRI_T1Gd_24,d$MRI_T2_24,use = "complete.obs")


# ---------------------- #
# MERGE DATA
# ---------------------- #
#
# Merge PASAT and brain vol
df = merge(PASAT,df,by=c("eid","session"),all=T)
# Further merge with EDSS
df = merge(edss,df,by=c("eid","session"),all=T)
# further merge with ALA
df = merge(ALA,df,by=c("eid","session"),all=T)
# further merge with new T1Gd lesions factor
df = merge(T1Gd,df,by=c("eid","session"),all=T)
# here is what the filtered version looks like
df %>% filter(session == 0 |session == 12 |session == 24)
df = merge(df,data.frame(eid = c(replicate(3,demo$Patnr)), 
                         session = c(replicate(nrow(demo),0),replicate(nrow(demo),12),replicate(nrow(demo),24)),
           DiseaseDuration = c(demo$DISEASE_DURATION, demo$DISEASE_DURATION+1,demo$DISEASE_DURATION+2)),
      by=c("eid","session"),all=T)

# further merge with new relapse
rel = merge(relapse,df,by=c("eid","session"),all=T)


# ---------------------- #
# SAVE DATA
# ---------------------- #

write.csv(df,file="/Users/max/Documents/Local/MS/ALA/data/clean.csv")
write.csv(rel,file="/Users/max/Documents/Local/MS/ALA/data/relapses_long.csv")
