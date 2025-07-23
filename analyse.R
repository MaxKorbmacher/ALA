# Replication of Serum Alpha-Linolenic Acid and Long-Term Multiple Sclerosis Activity and Progression
#
# Max Korbmacher, Jul 2025
#
# ------------------------------------------------------------ #
# ---------------------- Contents ---------------------------- #
# ------------------------------------------------------------ #
# 0. Preparations -------------------------------------------- #
# 1. Descriptives & Validation ------------------------------- #
# 2. Analyse ------------------------------------------------- #
# 2.1 Mixed linear models ------------------------------------ #
# 2.2 Mixed negative binomial models ------------------------- #
# 3. Simple linear models ------------------------------------ #
# 3.1 T2w number of lesions  --------------------------------- #
# 3.2 Brain volume  ------------------------------------------ #
# 3.3 PASAT  ------------------------------------------------- #
# 3.4 EDSS  -------------------------------------------------- #
# 3.5 T1w new Lesions (Negative binomial models) ------------- #
# 4. Mediation analyses -------------------------------------- #
# 4.1 Brain vol > ALA > EDSS---------------------------------- #
# 4.2 T1w lesions > ALA > EDSS-------------------------------- #
# 5. Longitudinal predictions of ALA ------------------------- #
# ------------------------------------------------------------ #
# ------------------------------------------------------------ #
#
# 0. Preparations --------------------------------------------
# define data path
datapath = "/Users/max/Documents/Local/MS/ALA/data/"

# read packages
pacman::p_load(haven,dplyr,reshape2,lme4,lmerTest,VGAM,mediation, psych, MuMIn, reshape2)
# load data
df = read.csv(paste(datapath,"clean.csv",sep=""))
long = read.csv(paste(datapath,"rate_of_change_10yrs.csv",sep=""))
#
# 1. Descriptives & Validation -------------------------------
ICC(reshape(df%>%dplyr::select(eid,session,ALA), 
            idvar = "eid", timevar = "session", 
            direction = "wide") %>% dplyr::select(-eid))
# 2. Analyse ------------------------------------------------- 
# 2.1 Mixed linear models ------------------------------------ 
# 2.1.1 PASAT
m = lmer(PASAT ~ log(ALA) + (1|eid),df)
summary(m)
m = lmer(PASAT ~ log(ALA) + age + sex + Treatment_OFAMS + (1|eid),df)
summary(m)
# 2.1.2 brain volume ***
m = lmer(TotalVol ~ log(ALA) + EstimatedTotalIntraCranialVol + (1|eid),df)
summary(m)
effectsize::standardize_parameters(m)
r.squaredGLMM(m)

m = lmer(TotalVol ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol + (1|eid),df)
summary(m)
effectsize::standardize_parameters(m)
r.squaredGLMM(m)

# # 2.1.2.1 GM volume *
# m = lmer(TotalGrayVol ~ log(ALA) + EstimatedTotalIntraCranialVol + (1|eid),df)
# summary(m)
# effectsize::standardize_parameters(m)
# m = lmer(TotalGrayVol ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol + (1|eid),df)
# summary(m)
# effectsize::standardize_parameters(m)
# 
# # 2.1.2.2 WM volume
# m = lmer(TotalWMVol ~ log(ALA) + EstimatedTotalIntraCranialVol + (1|eid),df)
# summary(m)
# effectsize::standardize_parameters(m)
# m = lmer(TotalWMVol ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol + (1|eid),df)
# summary(m)
# effectsize::standardize_parameters(m)

# 2.1.3 EDSS
m = lmer(edss ~ log(ALA) + (1|eid),df)
summary(m)
effectsize::standardize_parameters(m)
r.squaredGLMM(m)
m = lmer(edss ~ log(ALA) + age + sex + Treatment_OFAMS + (1|eid),df)
summary(m)
effectsize::standardize_parameters(m)
r.squaredGLMM(m)

# 2.1.3 T2w number of lesions
m = lmer(lesion_count ~ log(ALA) + (1|eid),df)
summary(m)
effectsize::standardize_parameters(m)
# m = lmer(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS + (1|eid),df)
# summary(m)
# control also for intracranial volume
m = lmer(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol + (1|eid),df)
summary(m)
effectsize::standardize_parameters(m)

# 2.2 Mixed negative binomial models --------------------------
# 2.2.1 T1w new Lesions (Negative binomial models) ---------------------
m = glmer(new_T1Gd_lesion ~ log(ALA) + (1|eid),df,family = binomial(link = cloglog))
summary(m)
m = glmer(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS + (1|eid),df,family = binomial(link = cloglog))
summary(m)
# control also for intracranial volume
m = glmer(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol + (1|eid),df,family = binomial(link = cloglog))
summary(m)

# 3. Simple linear models ------------------------------------ 

# 3.1 T2w number of lesions  -----------------------------------------
# including ICV
m = lm(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==0))
summary(m)
effectsize::standardize_parameters(m)
m = lm(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==12))
summary(m)
m = lm(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==24))
summary(m)
# not including ICV
m = lm(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS ,df%>%filter(session==0))
summary(m)
m = lm(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==12))
summary(m)
m = lm(lesion_count ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==24))
summary(m)

# 3.2 Brain volume  -----------------------------------------
# including ICV
m = lm(TotalVol ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==0))
summary(m)
m = lm(TotalVol ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==12))
summary(m)
m = lm(TotalVol ~ log(ALA) + age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==24))
summary(m)
# not including ICV
m = lm(TotalVol ~ log(ALA) + age + sex + Treatment_OFAMS ,df%>%filter(session==0))
summary(m)
m = lm(TotalVol ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==12))
summary(m)
m = lm(TotalVol ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==24))
summary(m)


# 3.3 PASAT  -----------------------------------------
m = lm(PASAT ~ log(ALA),df%>%filter(session==0))
summary(m)
m = lm(PASAT ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==0))
summary(m)
m = lm(PASAT ~ log(ALA),df%>%filter(session==12))
summary(m)
m = lm(PASAT ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==12))
summary(m)
m = lm(PASAT ~ log(ALA),df%>%filter(session==24))
summary(m)
m = lm(PASAT ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==24))
summary(m)

# 3.4 EDSS  -----------------------------------------
m = lm(edss ~ log(ALA),df%>%filter(session==0))
summary(m)
m = lm(edss ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==0))
summary(m)
m = lm(edss ~ log(ALA),df%>%filter(session==12))
summary(m)
effectsize::standardize_parameters(m)
m = lm(edss ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==12))
summary(m)
m = lm(edss ~ log(ALA),df%>%filter(session==24))
summary(m)
effectsize::standardize_parameters(m)
m = lm(edss ~ log(ALA) + age + sex + Treatment_OFAMS,df%>%filter(session==24))
summary(m)
effectsize::standardize_parameters(m)
# 3.5 T1w new Lesions (Negative binomial models) ---------------------
# BL
m = glm(new_T1Gd_lesion ~ log(ALA),df%>%filter(session==0),family = binomial(link = cloglog))
summary(m)
m = glm(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS ,df%>%filter(session==0),family = binomial(link = cloglog))
summary(m)
# control also for intracranial volume
m = glm(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==0),family = binomial(link = cloglog))
summary(m)
# 12 months
m = glm(new_T1Gd_lesion ~ log(ALA),df%>%filter(session==12),family = binomial(link = cloglog))
summary(m)
m = glm(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS ,df%>%filter(session==12),family = binomial(link = cloglog))
summary(m)
# control also for intracranial volume
m = glm(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==12),family = binomial(link = cloglog))
summary(m)
# 24 months
m = glm(new_T1Gd_lesion ~ log(ALA),df%>%filter(session==24),family = binomial(link = cloglog))
summary(m)
m = glm(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS ,df%>%filter(session==24),family = binomial(link = cloglog))
summary(m)
# control also for intracranial volume
m = glm(new_T1Gd_lesion ~ log(ALA) +age + sex + Treatment_OFAMS + EstimatedTotalIntraCranialVol,df%>%filter(session==24),family = binomial(link = cloglog))
summary(m)

# 4. Mediation analyses --------------------------------------
detach(package:lmerTest,unload = T) # lmerTest needs to leave for these functions to work
# 4.1 Brain vol > ALA > EDSS----------------------------------
# select and log transform
df2 = df %>% dplyr::select(eid,age,sex,Treatment_OFAMS,ALA,edss,TotalVol,lesion_count) #%>% na.omit
df2$ALA = log(df2$ALA)
# models
fit.mediator = lmer(ALA ~ age + sex + Treatment_OFAMS + TotalVol + (1|eid),df2)
fit.dv = lmer(edss ~ age + sex + Treatment_OFAMS + ALA + TotalVol + (1|eid),df2)
# mediation
results <- mediation::mediate(fit.mediator, fit.dv, mediator='ALA', treat = 'TotalVol')
summary(results)

# 4.2 T1w lesions > ALA > EDSS--------------------------------
# models
fit.mediator = lmer(ALA ~ age + sex + Treatment_OFAMS + lesion_count + (1|eid),df2)
fit.dv = lmer(edss ~ age + sex + Treatment_OFAMS + ALA + lesion_count + (1|eid),df2)
# mediation
results <- mediation::mediate(fit.mediator, fit.dv, mediator='ALA', treat = 'lesion_count')
summary(results)

# 5. Longitudinal predictions of ALA -------------------------
# prep df
long.df = merge(df %>% dplyr::filter(session == 0)%>%dplyr::select(eid,ALA,EstimatedTotalIntraCranialVol,Treatment_OFAMS),long, id.vars = "eid")
# check associations of ALA with THE ANNUAL RATE OF CHANGE in ...
## EDSS
m=lm(EDSS_diff~ALA,data=long.df)
summary(m)
m=lm(EDSS_diff ~ALA+age + sex + Treatment_OFAMS,data=long.df)
summary(m)

## PASAT
m=lm(PASAT_diff~ALA,data=long.df)
summary(m)
m=lm(PASAT_diff ~ALA+age + sex + Treatment_OFAMS,data=long.df)
summary(m)

## T2w lesions
m=lm(lesion_count_diff~ALA,data=long.df)
summary(m)
m=lm(lesion_count_diff ~ALA+age + sex + Treatment_OFAMS+EstimatedTotalIntraCranialVol,data=long.df)
summary(m)

## Brain volume
m=lm(TotalVol_diff~ALA,data=long.df)
summary(m)
m=lm(TotalVol_diff ~ALA+age + sex + Treatment_OFAMS+EstimatedTotalIntraCranialVol,data=long.df)
summary(m)

## [[total]] T1w lesions
m=lm(T1wLesions~ALA,data=long.df)
summary(m)
m=lm(T1wLesions ~ALA+age + sex + Treatment_OFAMS+EstimatedTotalIntraCranialVol,data=long.df)
summary(m)
