# Code for the analysis and visualization of data from the manuscript "Habitat alteration impacts predation risk in an aposematic amphibian"
# Code Authors: Doriane Hagnier and Carolin Dittrich 
# data:  forest_div.csv and attacks_data.csv
# R version 4.4.1 

#### Libraries to load####
library(MASS)
library(AICcmodavg)
library(lme4)
library(DHARMa)
library(reshape2)
library(effsize)
library(vegan)
library(gridExtra)
library(MuMIn)
library(tidyverse)

####forest data as csv file####

# data include parameters on location name, forest inventory number, core/managed zone, number of strata, percent per strata
# mean age per strata, shannon and simpson index, number of tree species, number of attacks by birds and mammals, all per forest patch
envi<-read.csv("forest_div.csv",header=TRUE, sep=";")
head(envi)
str(envi)


#### Welsh two way t-tests for differences between forest zones ####

t.test(envi$nspecies~envi$core) #significant with p=0.025

t.test(envi$simp~envi$core) # not significant, only trend p-value = 0.1124, mean in group core 0.2477286  mean in group managed 0.3967104 

t.test(envi$shan~envi$core)# not significant, trend p-value = 0.08424, mean in group core 0.4316020 mean in group managed 0.7335435

t.test(envi$strata1age~envi$core) #not significant, p=0.321, mean in group core 94.92857 mean in group managed 78.00000 

#### Hedges g as effect size (Cohens d with hedges correction) ####

cohen.d(nspecies~core, data = envi, hedges.correction=TRUE) # large -1.05 (CI 1-93 to -0.17)

cohen.d(simp~core, data = envi, hedges.correction=TRUE) # medium 0.65 (CI -1.5 to 0.20)

cohen.d(shan~core, data = envi, hedges.correction=TRUE) # medium -0.74 (CI -1.59 to 0.12)

cohen.d(strata1age~core, data = envi, hedges.correction=TRUE) # small 0.36 (CI -0.47 to 1.2)


#### Figure 3 with ggplot2####
# start with number of tree species and assign names to the plots for printing 

Fig3a<- ggplot(envi, aes(x = nspecies, y = core)) +
  geom_jitter(height= 0.05, width= 0.05, size=3.5, shape=21, aes(fill= core)) +
  stat_summary(fun.data = mean_cl_boot,
               geom= "pointrange",
               size= 2.0,
               shape= 21,
               fill = "white")+
  scale_x_continuous(breaks = seq(1,7))+
  scale_y_discrete(labels=c("protected", "managed"))+
  xlab("number of tree species")+
  ylab("forest zone")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position = "none",
        plot.margin = unit(c(1, 0.5, 1, 1), "lines"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))
# check figure
Fig3a

# plot shannon index per forest patch in the two forest zones
Fig3b<- ggplot(envi, aes(x = shan, y = core)) +
  geom_jitter(height= 0.01, width= 0.01, size=3.5, shape=21, aes(fill= core)) +
  stat_summary(fun.data = mean_cl_boot,
               geom= "pointrange",
               size= 2.0,
               shape= 21,
               fill = "white")+
  scale_y_discrete(labels=c("protected", "managed"))+
  xlab("Shannon index")+
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=22),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 0.5), "lines"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))
# check figure
Fig3b

# plot simpson index per forest patch in the two forest zones
Fig3c<-ggplot(envi, aes(x = simp, y = core)) +
  geom_jitter(height= 0.01, width= 0.01, size=3.5, shape=21, aes(fill= core)) +
  stat_summary(fun.data = mean_cl_boot,
               geom= "pointrange",
               size= 2.0,
               shape= 21,
               fill = "white")+
  scale_y_discrete(labels=c("protected", "managed"))+
  xlab("Simpson index of diversity")+
  ylab("forest zone")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position = "none",
        plot.margin = unit(c(1, 0.5, 1, 1), "lines"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))
# check figure
Fig3c

# plot mean age of first strata per forest patch in the two forest zones
Fig3d<-ggplot(envi, aes(x = strata1age, y = core)) +
  geom_jitter(height= 0.05, width= 0.05, size=3.5, shape=21, aes(fill= core)) +
  stat_summary(fun.data = mean_cl_boot,
               geom= "pointrange",
               size= 2.0,
               shape= 21,
               fill = "white")+
  scale_x_continuous(breaks = seq(10,230, by=20))+
  scale_y_discrete(labels=c("protected", "managed"))+
  xlab("age of the first tree strata")+
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=22),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 0.5), "lines"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))
#check figure
Fig3d

##### print/save figures####
# combine all four panels into one figure 
Fig3<- grid.arrange(Fig3a, Fig3b, Fig3c,Fig3d, nrow = 2)
ggsave("Figure3.png", plot = Fig3, device ="png", width = 25, height =17, units = c("cm"), dpi = 300)


####Attacks data an salamander plasticine models, cleaned for missing models, rodent attacks included####
attacks<-read.csv("attacks_data.csv", sep=";", h=T)
head(attacks)
str(attacks)

# Considering treatment, predator type, core as factors
attacks$treatm<-as.factor(attacks$treatm)
attacks$predtype<-as.factor(attacks$predtype)
attacks$core<-as.factor(attacks$core)

summary(attacks)

# In this study, we had many rodent attacks that were most probably not due to predation
# We decided to remove them for the analyses and provide a rational here that inclusion would not change the main finding of this study

count(attacks,core,rodent,treatm)
# rodents did "attack" models with small and large markings equally for both forest zones
count(attacks,rodent)
#### Generalized linear mixed models ##### 

# We lack the statistical power to model both predator types together, thus we build the models separately
# using cbind to account for the time (days) the models have been deployed, this varies from 4-5 days

##### GLMM Models:bird attacks probability as response, Rodent as non events####

# full model with interaction between core and treatment, location as random effect
mod1<-glmer(cbind(birdatt,time)~core*treatm+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod1)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod1, plot = T)

# core  + treatment, location as random effect
mod2<-glmer(cbind(birdatt,time)~core+treatm+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod2)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod2, plot = T)

# just core, location as random effect
mod3<-glmer(cbind(birdatt,time)~core+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod3)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod3, plot = T)

# just treatment, location as random effect
mod4<-glmer(cbind(birdatt,time)~treatm+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod4)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod4, plot = T)

# just random effect location (null model)
mod5<- glmer(cbind(birdatt, time)~1+ (1|location), data= attacks, family = binomial)
summary(mod5)
# effect of the transect more than the area => 0 in location 
simulateResiduals(fittedModel = mod5, plot = T)


# Check the difference in AIC to choose best model
models <- list(mod1, mod2, mod3, mod4, mod5)
aictab(cand.set = models)


##### Models:mammal attack probability, Rodent as non events####
# full model with interaction between core and treatment, location as random effect
mod6<-glmer(cbind(mammatt,time)~core*treatm+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod6)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod6, plot = T)

# core  + treatment, location as random effect
mod7<-glmer(cbind(mammatt,time)~core+treatm+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod7)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod7, plot = T)

# just core, random factor location
mod8<-glmer(cbind(mammatt,time)~core+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod8)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod8, plot = T)

# just treatm, random factor location
mod9<-glmer(cbind(mammatt,time)~treatm+(1|location), na.action=na.omit, data=attacks, family=binomial)
summary(mod8)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod8, plot = T)

# just random factor location
mod10<- glmer(cbind(mammatt, time)~1+ (1|location), data= attacks, family = binomial)
summary(mod10)
# check models visually with the DHARMa package
simulateResiduals(fittedModel = mod10, plot = T)


# Check the difference in AIC
models <- list(mod6, mod7, mod8, mod9, mod10)
aictab(cand.set = models)


####Permutations####
# Here is a function allowing to permute conditions for binary variables 
# The function will compare our observed result (difference in proportion) to a simulated null model to increase statistical power

SuperPermProp <- function(fixval, cond.permut, your.data) {
  obs_tab <- table(fixval, cond.permut)  # observed contingency table
  
  obs_diff <- obs_tab[2, 1] / sum(obs_tab[, 1]) - obs_tab[2, 2] / sum(obs_tab[, 2])  # observed difference
  
  dal <- numeric(10000)  # permutation box
  for (i in 1:9999) {
    your.data$al <- sample(cond.permut)
    perm_tab <- table(fixval, your.data$al)  # contingency table for permutation
    
    dal[i] <- perm_tab[2, 1] / sum(perm_tab[, 1]) - perm_tab[2, 2] / sum(perm_tab[, 2])
  }  # permutation of the condition habitat or morph
  
  dal[10000] <- obs_diff  # putting the observed difference at the end of the permutation box
  
  # counts how many times the difference in your.data$cond.permutation box is >= to the observed difference
  tap <- table(abs(dal) >= abs(obs_diff))
  
  pvalue <- tap[2] / 10000
  pvalue
  
  limsum <- quantile(dal, c(0.025, 0.975), na.rm = TRUE)
  limsum
  
  dal[10000]
  
  par(mfrow = c(1, 1))
  hist(dal, breaks = max(1, min(30, length(unique(dal)))), main = "Histogram of p.value")
  abline(v = limsum[1], col = "red")
  abline(v = limsum[2], col = "red")
  abline(v = dal[10000], lty = 2, col = "blue")
  
  return(pvalue)
}

# For bird attacks 
# core
SuperPermProp(attacks$birdatt, attacks$core, attacks)
# treatm
SuperPermProp(attacks$birdatt, attacks$treatm, attacks)

# For mammal attacks
# Rodent as non events
# core
SuperPermProp(attacks$mammatt, attacks$core, attacks)
#treatm
SuperPermProp(attacks$mammatt, attacks$treatm, attacks)

# for headattacks
# rodents as non event 
# core
SuperPermProp(attacks$headatt, attacks$core, attacks)
#treatm
SuperPermProp(attacks$headatt, attacks$treatm, attacks)
# predator type 
SuperPermProp(attacks$headatt, attacks$predtype, attacks)

#### Figure 4 ####
# Percentage of attacks by different predators (birds = green, mammals = brown) in function of forest zone (managed vs protected)
# combine mammal and bird attacks
mambird<-subset(attacks, mammatt == 1 | birdatt ==1) # we make a dataset for just mammal and bird attacks
mambird <- mambird %>% 
  mutate(core = recode(core, `0` = "managed", `1` = "protected"))

# plot mammal and bird attacks by environment
# filter out the two various observations 
filtered_df <- mambird %>% filter(predtype != "various")

Fig4<-ggplot(filtered_df)+
  geom_bar(aes(x=core, fill= predtype), position="dodge")+
  scale_fill_manual(values = c("#22521A","#B2592C")) +
  theme_bw() +
  labs(y = "Percentage of attacks", x = "forest zone", fill = "Predator type") +
  coord_cartesian(ylim = c(0, 40)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"))+
  theme(text = element_text(size = 18),
        legend.position = "none")

ggsave("Figure4.png", plot = Fig4, device ="png", width = 15, height =10, units = c("cm"), dpi = 300)

#### Generalized models of number of bird attacks per forest patch and diversity index####

# linear model 1 with core and shannon index as predictors
lm1a<-lm(n_attack_bird ~ core * shan, data=envi)
summary(lm1a)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(lm1a, n = 1000, plot = T)

lm1b<-lm(n_attack_bird ~ core + shan, data=envi)
summary(lm1b)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(lm1b, n = 1000, plot = T)

# compare model with and without interaction term 
AIC(lm1a, lm1b) # better without interaction
anova(lm1a,lm1b) # almost no difference
summary(lm1a)$adj.r.squared
summary(lm1b)$adj.r.squared # better without interaction

# use negative binomial distribution for the response with shannon index and forest zone
# with interaction
nb1a<- MASS::glm.nb(n_attack_bird~shan*core, data=envi)
summary(nb1a)
simulation_output1a <- simulateResiduals(fittedModel = nb1a)
# Plot residual diagnostics
plot(simulation_output1a) # looks god

#without interaction
nb1b<- MASS::glm.nb(n_attack_bird~shan+core, data=envi)
summary(nb1b)
simulation_output1b <- simulateResiduals(fittedModel = nb1b)
# Plot residual diagnostics
plot(simulation_output1b) # looks even better

# compare model fit 
AIC(nb1a, nb1b) #AIC slightly better without interaction interaction 
anova(nb1a, nb1b)

# Calculate 95% confidence intervals for the coefficients and estimeates
est <- cbind(Estimate = coef(nb1b), confint(nb1b))
# use exponents of estimate and CI tbe used as incident rate ratios (IRR)
exp(est)



####simpson 
# use linear distribution for the response with simpson index and forest zone
# with interaction 
lm2a<-lm(n_attack_bird ~ core * simp, data=envi)
summary(lm2a)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(lm2a, n = 1000, plot = T)

#without interaction
lm2b<-lm(n_attack_bird ~ core + simp, data=envi)
summary(lm2b)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(lm2b, n = 1000, plot = T)

# compare model with and without interaction term 
AIC(lm2a, lm2b) # better without interaction
anova(lm2a,lm2b) # almost no difference
summary(lm2a)$adj.r.squared
summary(lm2b)$adj.r.squared # better 

# model without interaction is better for bird attack and diversity indices
# 
# use negative binomial distribution for the response with simpson index and forest zone
# with interaction
nb2a<- MASS::glm.nb(n_attack_bird~simp*core, data=envi)
summary(nb2a)
simulation_output2a <- simulateResiduals(fittedModel = nb2a)
# Plot residual diagnostics
plot(simulation_output2a) # looks god

#without interaction
nb2b<- MASS::glm.nb(n_attack_bird~simp+core, data=envi)
summary(nb2b)
simulation_output2b <- simulateResiduals(fittedModel = nb2b)
# Plot residual diagnostics
plot(simulation_output2b) # looks even better

# compare model fit 
AIC(nb2a, nb2b) #AIC slightly better without interaction interaction 
anova(nb2a, nb2b)

# Calculate 95% confidence intervals for the coefficients and estimeates
(est2 <- cbind(Estimate = coef(nb2b), confint(nb2b)))
# use exponents of estimate and CI tbe used as incident rate ratios (IRR)
exp(est2)

####Figure 5 bird attacks and diversity indices of forest patches in the two forest zones####

Fig5a<-ggplot(envi, aes(x = shan, y = n_attack_bird))+
  geom_point(size=3.5, shape=21, aes(fill= coref))+
  geom_smooth(method = glm.nb, formula= y~x, se=T, aes(colour=coref))+
  scale_y_continuous(limits= c(0,11),breaks=seq(0,11))+
  xlab("Shannon index")+
  ylab("Number of bird attacks")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position = "none",
        plot.margin = unit(c(1, 0.5, 1, 1), "lines"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))+
  scale_colour_manual(values = c("#FFB90F", "#104E8B"))

#check
Fig5a

Fig5b<-ggplot(envi, aes(x = simp, y = n_attack_bird))+
  geom_point(size=3.5, shape=21, aes(fill= coref))+
  geom_smooth(method = glm.nb, formula= y~x, se=T, aes(colour=coref))+
  scale_y_continuous(limits= c(0,11),breaks=seq(0,11, by=1))+
  xlab("Simpson index of diversity")+
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 0.5), "lines"))+
  scale_colour_manual(values = c("#FFB90F", "#104E8B"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))

#check
Fig5b


#### Linear models of number of mammal attacks per forest patch and diversity index####

# model 3 with core and shannon index as predictors
m3a<-lm(n_attack_mam ~ core * shan, data=envi)
summary(m3a)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(m3a, n = 1000, plot = T) # not a good fit 
simulateResiduals(m3a, quantreg=T, plot = T)
plotResiduals(m3a) # not good fit
testDispersion(m3a)

m3b<-lm(n_attack_mam ~ core + shan, data=envi)
summary(m3b)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(m3b, n = 1000, plot = T) # not a good fit
simulateResiduals(m3b, quantreg=T, plot = T)
plotResiduals(m3b) # not good fit
testDispersion(m3b)

# compare model with and without interaction term 
AIC(m3a, m3b) # better without interaction
anova(m3a,m3b) # 
summary(m3a)$adj.r.squared # better
summary(m3b)$adj.r.squared  

#negative binomial model with MASS package
nb3a<-MASS::glm.nb(n_attack_mam~shan*core, data=envi)
summary(nb3a)
simulation_output3a <- simulateResiduals(fittedModel = nb3a)
# Plot residual diagnostics
plot(simulation_output3a) # looks way better to fit the data

#negative binomial model with MASS package
nb3b<-MASS::glm.nb(n_attack_mam~shan+core, data=envi)
summary(nb3b)
simulation_output3b <- simulateResiduals(fittedModel = nb3b)
# Plot residual diagnostics
plot(simulation_output3b) #fits less, interaction term seems to be important 

AIC(nb3a, nb3b) #AIC slightly better for interaction 
anova(nb3a, nb3b) # interaction seems better fit

# Calculate 95% confidence intervals for the coefficients and estimeates
(est3 <- cbind(Estimate = coef(nb3a), confint(nb3a)))
# use exponents of estimate and CI tbe used as incident rate ratios (IRR)
exp(est3)


# linear model 4 with core and simpson index as predictors and interaction 
lm4a<-lm(n_attack_mam ~ core * simp, data=envi)
summary(lm4a)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(lm4a, n = 1000, plot = T) # not a good fit, shows deviation from model assumptions

# without interaction
lm4b<-lm(n_attack_mam ~ core + simp, data=envi)
summary(lm4b)
# check with DHARMa and 1000 simulations to produce more stable results
simulateResiduals(lm4b, n = 1000, plot = T) # not a good fit, shows deviation from model assumptions

# compare model with and without interaction term for better fit
AIC(lm4a, lm4b) # better without interaction
anova(lm4a,lm4b) # 
summary(lm4a)$adj.r.squared # better with interaction
summary(lm4b)$adj.r.squared 

#linear models shows deviation and poisson overdispersion, use negative binomial model with MASS package
nb4a<-MASS::glm.nb(n_attack_mam~simp*core, data=envi)
summary(nb4a)
simulation_output4a <- simulateResiduals(fittedModel = nb4a)
# Plot residual diagnostics
plot(simulation_output4a) # better but deviation detected

#negative binomial model with MASS package
nb4b<-MASS::glm.nb(n_attack_mam~simp+core, data=envi)
summary(nb4b)
simulation_output4b <- simulateResiduals(fittedModel = nb4b)
# Plot residual diagnostics
plot(simulation_output4b) #fits better, no problems detected 

# test models for better fit
AIC(nb4a, nb4b) #AIC slightly better for interaction 
anova(nb4a, nb4b) # interaction seems better

# Calculate 95% confidence intervals for the coefficients and estimeates
(est4 <- cbind(Estimate = coef(nb4a), confint(nb4a)))
# use exponents of estimate and CI tbe used as incident rate ratios (IRR)
exp(est4)

####Figure 5c-d mammal attacks and diversity indices of forest patches in the two forest zones####

Fig5c<-ggplot(envi, aes(x = shan, y = n_attack_mam))+
  geom_point(size=3.5, shape=21, aes(fill= coref))+
  geom_smooth(method = glm.nb, formula= y~x, se=T, aes(colour=coref))+
  scale_y_continuous(breaks=seq(0,11))+
  xlab("Shannon index")+
  ylab("Number of mammal attacks")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position = "none",
        plot.margin = unit(c(1, 0.5, 1, 1), "lines"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))+
  scale_colour_manual(values = c("#FFB90F", "#104E8B"))

#check
Fig5c

Fig5d<-ggplot(envi, aes(x = simp, y = n_attack_mam))+
  geom_point(size=3.5, shape=21, aes(fill= coref))+
  geom_smooth(method = glm.nb, formula= y~x, se=T, aes(colour=coref))+
  scale_y_continuous(limits=c(0,11), breaks=seq(0,11))+
  xlab("Simpson index of diversity")+
  ylab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=22),
        axis.title.y = element_text(size=22),
        legend.position = "none",
        plot.margin = unit(c(1, 1, 1, 0.5), "lines"))+
  scale_colour_manual(values = c("#FFB90F", "#104E8B"))+
  scale_fill_manual(values = c("#FFB90F", "#104E8B"))

#check
Fig5d

##### print/save figures####
# figure 5 combination of all 
Fig5<- grid.arrange(Fig5a, Fig5b, Fig5c, Fig5d, nrow=2)
ggsave("Figure5.png", plot = Fig5, device ="png", width = 25, height =20, units = c("cm"), dpi = 300)
