rm(list=ls())
setwd("C:/Users/Marlen Fröhlich/Documents/R/MS Interaction engine")
library(lme4)
library(car)

# Subsetting data and transforming variables

test.data <- read.table ("data_repair.csv", header=TRUE, sep=",", stringsAsFactors=TRUE)
test.data <- subset(test.data, Context %in%  c("pl","jt","gr"))

test.data$age.dep=as.numeric(test.data$age==levels(test.data$age)[2])
test.data$age.imm=as.numeric(test.data$age==levels(test.data$age)[3])
test.data$sex.code=as.numeric(test.data$sex==levels(test.data$sex)[2])
test.data$kinship.mk=as.numeric(test.data$kinship==levels(test.data$kinship)[1])
test.data$kinship.mo=as.numeric(test.data$kinship==levels(test.data$kinship)[2])
test.data$agediff.ol=as.numeric(test.data$agediff==levels(test.data$agediff)[1])
test.data$agediff.yo=as.numeric(test.data$agediff==levels(test.data$agediff)[3])
test.data$context.pl=as.numeric(test.data$Context==levels(test.data$Context)[5])


per.sig=as.data.frame(na.omit(test.data[, c("Persist", "Signal", "species","age","age.dep","age.imm", "kinship.mo", "context.pl",  "sex.code", "group", "setting", "ID_coder", "ID_obs", "Dyad", "ID_sign","ID_rec", "Context")]))
per.sig2=as.data.frame(na.omit(test.data[, c("Persist", "ASO", "PerType", "Signal", "species", "age","age.dep","age.imm",  "kinship.mo", "context.pl",  "sex.code", "group", "setting", "ID_coder", "ID_obs", "Dyad", "ID_sign","ID_rec", "Context")]))

rep.sig=as.data.frame(na.omit(test.data[, c("RepeatASO2", "Signal", "species","age","age.dep","age.imm",  "kinship.mo", "context.pl",  "sex.code", "group", "setting", "ID_coder", "ID_obs", "Dyad", "ID_sign","ID_rec", "Context")]))
rep.sig2=as.data.frame(na.omit(test.data[, c("RepeatASO2","ASO", "Signal", "species","age","age.dep","age.imm",  "kinship.mo", "context.pl",  "sex.code", "group", "setting", "ID_coder", "ID_obs", "Dyad", "ID_sign","ID_rec")]))

ela.sig=as.data.frame(na.omit(test.data[, c("ElabP", "Signal", "species","age","age.dep","age.imm",  "kinship.mo", "context.pl",  "sex.code", "group", "setting", "ID_coder", "ID_obs", "Dyad", "ID_sign","ID_rec", "Context")]))
ela.sig2=as.data.frame(na.omit(test.data[, c("ElabP","ASO", "Signal", "species","age","age.dep","age.imm",  "kinship.mo", "context.pl",  "sex.code", "group", "setting", "ID_coder", "ID_obs", "Dyad", "ID_sign","ID_rec")]))

aso.rep.sig =subset(rep.sig2, RepeatASO2 ==1)
aso.ela.sig =subset(ela.sig2, ElabP ==1)
aso.per.sig =subset(per.sig2, Persist ==1)

contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000000))

# collinearity check: max vif = 1.3 (persist), 1.3 (elab), aso.rep (1.3), aso.ela (1.3)
vif(lm(rnorm(nrow(ela.sig)) ~ setting + species  + kinship.mo +  age + sex.code + context.pl, data = ela.sig))


#model 1: redoings

#full model containing essential random slopes
mod.per = glmer(formula = Persist ~ setting * species  + kinship.mo +  age.dep + age.imm + sex.code + context.pl + 
                  (0+age.dep|ID_rec)+ (0+age.imm|ID_rec)+ 
                  (0+context.pl|ID_sign)+  (0+context.pl|group)+
                  +(1|ID_obs)+(1|ID_rec) +(1|group:ID_sign), family = binomial, data = per.sig,  control=contr)

#null model containing essential random slopes
null.per = glmer(formula = Persist ~  age.dep +  age.imm + sex.code +   context.pl +
                   (0+age.dep|ID_rec)+ (0+age.imm|ID_rec)+ 
                   (0+context.pl|ID_sign)+  (0+context.pl|group)+
                   +(1|ID_obs)+(1|ID_rec) +(1|group:ID_sign), family = binomial, data = per.sig,  control=contr)

#Check number of residuals
length(residuals(mod.per)) #3869
length(residuals(null.per)) #3869

#Likelihood ratio test (full-null model comparison)
as.data.frame(anova(null.per, mod.per, test="Chisq"))

#model output
round(summary(mod.per)$coefficients, 3)

#Likelihood ratio test (individual fixed effects)
drop1(mod.per,  test ="Chisq")

#Postdhoc test of interaction terms
require(emmeans)
lsm <- lsmeans(mod.per, ~ setting * species, adjust ="sidak")
lsm1 <- lsmeans(mod.per, list(pairwise ~ setting|species, pairwise ~ species|setting))
lsm1

# model 2: elaboration

mod.ela = glmer(formula = ElabP ~ setting * species + kinship.mo + age.dep + age.imm + sex.code + context.pl +
                  +(0+context.pl|ID_rec)+
                  +(1|ID_obs)+(1|ID_rec), family = binomial, data = ela.sig,  control=contr)


null.ela = glmer(formula = ElabP ~  age.dep + age.imm + sex.code + context.pl +
                   +(0+context.pl|ID_rec)+
                   +(1|ID_obs)+(1|ID_rec), family = binomial, data = ela.sig,  control=contr)

length(residuals(mod.ela)) #1207
length(residuals(null.ela)) #1207

as.data.frame(anova(null.ela, mod.ela, test="Chisq"))


round(summary(mod.ela)$coefficients, 3)

drop1(mod.ela,  test ="Chisq")

require(lsmeans)
lsm <- lsmeans(mod.ela, ~ setting * species, adjust ="sidak")
lsm1 <- lsmeans(mod.ela, list(pairwise ~ setting|species, pairwise ~ species|setting))
lsm1


# model 3: effectiveness repetition 

mod.aso.rep = glmer(formula = ASO ~ setting + species + kinship.mo + age.dep + age.imm + sex.code + context.pl +
                         (0+kinship.mo|ID_sign) +
                         (0+sex.code|ID_rec)+
                         (0+context.pl|group)+
                         +(1|ID_obs) +(1|group), family = binomial, data = aso.rep.sig,  control=contr)

null.aso.rep = glmer(formula = ASO ~  age.dep + age.imm + sex.code + context.pl +
                       (0+kinship.mo|ID_sign) +
                       (0+sex.code|ID_rec)+
                       (0+context.pl|group)+
                       +(1|ID_obs) +(1|group), family = binomial, data = aso.rep.sig,  control=contr)

length(residuals(mod.aso.rep)) #727
length(residuals(null.aso.rep)) #727
as.data.frame(anova(null.aso.rep, mod.aso.rep, test="Chisq"))

drop1(mod.aso.rep,  test ="Chisq")


round(summary(mod.aso.rep)$coefficients, 3)


# model 4: effectiveness elaboration: no interaction term


mod.aso.ela = glmer(formula = ASO ~ setting + species + kinship.mo + age.dep + age.imm + sex.code + context.pl +
                         (0+kinship.mo|ID_sign) +
                         (0+sex.code|ID_rec), family = binomial, data = aso.ela.sig,  control=contr)

null.aso.ela = glmer(formula = ASO ~  age.dep + age.imm + sex.code + context.pl +
                       (0+kinship.mo|ID_sign) +
                       (0+sex.code|ID_rec), family = binomial, data = aso.ela.sig,  control=contr)

length(residuals(mod.aso.ela)) #476
length(residuals(null.aso.ela)) #476

round(summary(mod.aso.ela)$coefficients, 3)

drop1(mod.aso.ela,  test ="Chisq")


###########################################################################################################################
# Get individual means for four response variables

mod.per_kin=aggregate(x=per.sig$Persist, by=per.sig[, c("ID_sign", "species", "setting", "kinship.mo")], FUN=mean)
mod.ela_kin=aggregate(x=ela.sig$ElabP, by=ela.sig[, c("ID_sign", "species", "setting", "kinship.mo")], FUN=mean)

mod.aso.rep_kin=aggregate(x=aso.rep.sig$ASO, by=aso.rep.sig[, c("ID_sign", "species", "setting", "kinship.mo")], FUN=mean)
mod.aso.ela_kin=aggregate(x=aso.ela.sig$ASO, by=aso.ela.sig[, c("ID_sign", "species", "setting", "kinship.mo")], FUN=mean)


# plotting effects of setting, species and kinship

library(ggplot2)
theme_marlen_ss <- theme(panel.background = element_blank(),
                         panel.border =element_rect(colour="black", fill=NA),
                         
                         plot.background = element_blank(),
                         panel.grid = element_blank(),
                         axis.line = element_line(colour ="black"),
                         axis.text.x = element_text (size = 12,colour= "black", family="sans"),
                         axis.text.y = element_text (size = 12,colour= "black", family="sans"),
                         axis.ticks.y = element_line(colour="black"),
                         axis.ticks.x = element_line(colour=NA),
                         axis.title.x = element_text(size = 12, vjust = -0.5, family="sans"),
                         axis.title.y = element_text(size = 12, vjust = 2, family="sans"),
                         legend.text=  element_text(size = 11, family="sans", margin = margin(t = 10)),
                         legend.key = element_blank(),
                         legend.position = "right",
                         legend.spacing.x = unit(0.2, 'cm'),
                         strip.text = element_text(size = 12))

levels(mod.per_kin$setting) <- c("Zoo", "Wild")
mod.per_kin$kinship.mo=as.factor(mod.per_kin$kinship.mo)
mod.per_kin$kinship.mo <- relevel(mod.per_kin$kinship.mo, "1")

levels(mod.ela_kin$setting) <- c("Zoo", "Wild")
mod.ela_kin$kinship.mo=as.factor(mod.ela_kin$kinship.mo)
mod.ela_kin$kinship.mo <- relevel(mod.ela_kin$kinship.mo, "1")

levels(mod.aso.ela_kin$setting) <- c("Zoo", "Wild")
mod.aso.ela_kin$kinship.mo=as.factor(mod.aso.ela_kin$kinship.mo)
mod.aso.ela_kin$kinship.mo <- relevel(mod.aso.ela_kin$kinship.mo, "1")

levels(mod.aso.rep_kin$setting) <- c("Zoo", "Wild")
mod.aso.rep_kin$kinship.mo=as.factor(mod.aso.rep_kin$kinship.mo)
mod.aso.rep_kin$kinship.mo <- relevel(mod.aso.rep_kin$kinship.mo, "1")


# Boxplot example: Number of gestural pursuits (Fig. S1)

dodge.posn <- position_dodge(.9)
mod.site <- ggplot(count.per_kin, aes(x = species, y = x))
mod.site + geom_boxplot(aes(fill = kinship.mo), width = 0.9) +
  geom_point(aes(fill = kinship.mo),position= dodge.posn, shape = 1, colour = "black", alpha = 0.5) +
  theme_marlen_ss +
  scale_y_continuous("Mean number of gestural redoings after communicative failure") +
  scale_x_discrete("Orang-utan species",
                   limits = c("Bor", "Sum"),
                   labels = c("Bornean", "Sumatran"))+
  scale_fill_manual(values=c("grey70", "seagreen3"),name="Interaction dyad",
                    breaks=c("1","0"),
                    labels=c("mother-offspring", "other dyad"))+
  facet_wrap(~setting)+
  stat_summary(fun=mean, geom="point",shape =23, fill ="black",aes(group=kinship.mo), position=position_dodge(.9), 
               color="black", size=3)

# model stability

source("C:/Users/Marlen Fröhlich/Documents/R/lm_course/glmm_stability.r")
mod.per.stab=glmm.model.stab(model.res=mod.per, contr=contr, ind.cases=F, para=T)
                              )
save.image("C:/Users/Marlen Fröhlich/Documents/R/orangutan_repair.RData")
