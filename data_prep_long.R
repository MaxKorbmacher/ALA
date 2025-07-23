# Analysis prep based on previous script
# Max Korbmacher, Jul 2025
#
#
# clean up
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
savepath = "/Users/max/Documents/Local/MS/ALA/data/" # define results/save/oputput path
# 0. Prep ####
# load packages and install if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,lme4,lmerTest,effects,effectsize,interactions,gamm4,
               ggseg,ggtext,ggpubr,MuMIn,dplyr,ggplot2,standardize,haven,reshape2,
               simr, data.table)
# data
df = read.csv("/Users/max/Documents/Local/MS/data/final_dat_subc.csv")
T2lesions = read.csv("/Users/max/Documents/Local/MS/data/lesion_count.csv")
# standardize and join datasets
# df$scanner = ifelse(df$eid > 100 & df$eid < 200, "a", df$eid)
# df$scanner = ifelse(df$eid > 200 & df$eid < 300, "b", df$eid)
# df$scanner = ifelse(df$eid > 300 & df$eid < 400, "c", df$eid)
# df$scanner = ifelse(df$eid > 400 & df$eid < 500, "d", df$eid)
# df$scanner = ifelse(df$eid > 500 & df$eid < 600, "e", df$eid)
# df$scanner = ifelse(df$eid > 600 & df$eid < 700, "f", df$eid)
# df$scanner = ifelse(df$eid > 700 & df$eid < 800, "g", df$eid)
# df$scanner = ifelse(df$eid > 800 & df$eid < 900, "h", df$eid)
# df$scanner = ifelse(df$eid > 900 & df$eid < 1000, "i", df$eid)
# df$scanner = ifelse(df$eid > 1000 & df$eid < 1100, "j", df$eid)
# df$scanner = ifelse(df$eid > 1100 & df$eid < 1200, "k", df$eid)
# df$scanner = ifelse(df$eid > 1200 & df$eid < 1300, "l", df$eid)
# df$scanner = ifelse(df$eid > 1300 & df$eid < 1400, "m", df$eid)
# df$scanner = ifelse(df$eid > 1400 & df$eid < 1500, "n", df$eid)
# df$scanner = ifelse(df$eid > 1500 & df$eid < 1600, "o", df$eid)
# df$scanner = ifelse(df$eid > 1600, "p", df$eid)
# df$data = "OFAMS"
df$sex = factor(df$sex) # make sex a factor
levels(df$sex) = c("M","F")

df$TotalVol = df$TotalGrayVol + df$CerebralWhiteMatterVol + df$Right.Cerebellum.White.Matter + df$Left.Cerebellum.White.Matter
df$TotalWMVol = df$CerebralWhiteMatterVol + df$Right.Cerebellum.White.Matter + df$Left.Cerebellum.White.Matter

df = df %>% select(eid, session, sex, age, EstimatedTotalIntraCranialVol,TotalVol,TotalGrayVol,TotalWMVol)
df$session = df$session+1
#
#
# Add PASAT
demo10 = read.csv('/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv',sep = ";") # 10 years follow up demo & clin scores
pasat1 = demo10%>%dplyr::select(Patnr,BL_PASATcorrect,PASAT_24M,PASAT_OFAMS10)
pasat1 = reshape2::melt(pasat1, id.vars = "Patnr")
names(pasat1) = c("eid","session","PASAT")
pasat1$session = ifelse(pasat1$session == "BL_PASATcorrect",1,0)+ifelse(pasat1$session == "PASAT_24M",24,0)+ifelse(pasat1$session == "PASAT_OFAMS10",145,0)

df = merge(df,pasat1,by = c("eid","session"),all=T)

# add T1w lesions
t1wlesions = demo10 %>% select(Patnr, MRI_T1Gd_0, MRI_T1Gd_12, MRI_T1Gd_24, MRI_T1Gd_sum)
t1wlesions = reshape2::melt(t1wlesions, id.vars = "Patnr")
names(t1wlesions) = c("eid","session","T1wLesions")
t1wlesions$session = ifelse(t1wlesions$session == "MRI_T1Gd_0",1,0)+ifelse(t1wlesions$session == "MRI_T1Gd_12",24,0)+
  ifelse(t1wlesions$session == "MRI_T1Gd_24",24,0)+ifelse(t1wlesions$session == "MRI_T1Gd_sum",145,0)
df = merge(df,t1wlesions,by = c("eid","session"),all=T)

# add T2w lesions
t2wlesions = reshape2::melt(T2lesions, id.vars = "eid")
names(t2wlesions) = c("eid","session","t2wLesions")
t2wlesions$session = gsub("m","",t2wlesions$session)
t2wlesions$session = as.numeric(ifelse(t2wlesions$session == "baseline",0,t2wlesions$session))+1
t2wlesions$session = ifelse(t2wlesions$session == 121,145,t2wlesions$session)
df = merge(df,t2wlesions,by = c("eid","session"),all=T)


# Add EDSS
EDSS = demo10 %>% select(Patnr, edss_baseline, edss_month_12, edss_month18, edss_month_24, EDSS_score_10)
EDSS = reshape2::melt(EDSS, id.vars = "Patnr")
names(EDSS) = c("eid","session","EDSS")
EDSS$session = ifelse(EDSS$session == "edss_baseline",1,0)+ifelse(EDSS$session == "edss_month_12",24,0)+
  ifelse(EDSS$session == "edss_month18",18,0)+
  ifelse(EDSS$session == "edss_month_24",24,0)+ifelse(EDSS$session == "EDSS_score_10",145,0)
df = merge(df,EDSS,by = c("eid","session"),all=T)
df$EDSS = as.numeric(gsub(",",".",df$EDSS))

# remove duplicates (resulting from repeat scans)
X = as.data.table(df)
keys = df %>% select(eid,session,sex) %>% names
X1 = X %>% filter(session == 1)
X2 = X %>% filter(session == 24)
X3 = X %>% filter(session == 145)
X1 = X1[,lapply(.SD,mean),keys]
X2 = X2[,lapply(.SD,mean),keys]
X3 = X3[,lapply(.SD,mean),keys]
df = rbind(df %>% filter(session != 1 & session != 24 & session != 145),X1)
df = rbind(df,X2)
df = rbind(df,X3)
# write the long frame
write.csv(x = df, file = paste(savepath, "long_10yrs.csv",sep = ""))
# calculate rate of change
df = df %>% filter(session == 1 | session == 145)

df = df %>% 
  arrange(eid, t2wLesions, EstimatedTotalIntraCranialVol, TotalVol, 
          TotalGrayVol, TotalWMVol, PASAT, EDSS) %>% 
  group_by(eid) %>% 
  mutate(lesion_count_diff = c(diff(t2wLesions)/12, NA),
         EstimatedTotalIntraCranialVol_diff = c(diff(EstimatedTotalIntraCranialVol)/12, NA),
         TotalVol_diff = c(diff(TotalVol)/12, NA),
         TotalGrayVol_diff = c(diff(TotalGrayVol)/12, NA),
         TotalWMVol_diff = c(diff(TotalWMVol)/12, NA),
         PASAT_diff = c(diff(PASAT)/12, NA),
         EDSS_diff = c(diff(EDSS)/12, NA))
df_copy=df
df = df %>% filter(session == 1) %>% select(eid, sex, age, ends_with("_diff"))
df = merge(df,df_copy %>% filter(session == 145) %>% select(eid,T1wLesions),by="eid")

write.csv(x = df, file = paste(savepath, "rate_of_change_10yrs.csv",sep = ""))
