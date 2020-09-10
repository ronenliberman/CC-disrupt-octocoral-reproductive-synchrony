require(coxme)
library(purrr)
library(tidyr)
library(tidyverse)
library(survminer)
library(survival)
library(ggsci)
library(cowplot)

d <- read.csv("trial_best.csv")
str(d)
head(d)

##using this method to transfrom my counted data into a bimal data with dead planula=1
#and living ones at day 7 = 0

library(splitstackshape)
oo <- expandRows(d, "count")
View(oo)
str(oo)

#check the data has transformed well
fit_oo<- survfit(Surv(day, status) ~ treatment, data=oo)

summary(fit_oo)
#visualize the data 

survplot <- ggsurvplot_facet(fit_oo, oo, facet.by = "source",
                 palette = "jco", pval = TRUE, legend = "bottom",
  xlab="Days post release", legend.labs = c("Ambeint","RCP 4.5", "RCP 8.5"),panel.labs = list(source = c("RSS-derived embryos", "Reef-derived embryos")))+
  theme_classic() %+replace% 
  theme(legend.position = "bottom",
        axis.line = element_line(colour = "black"),legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor.y=element_line(colour = "gray"),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        axis.title =element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(size=12))

survplot

ggsave(plot= survplot, "survplot.tiff", height = 8, width = 16, units = "cm", dpi = 600)

#model with fixed == source +treatment Random = aquaria 
simp_model_all <- coxph(Surv(day, status) ~ source+treatment,data=oo)
summary(simp_model_all) 

##model with all the data and random effect
 mix_eff<- coxme(Surv(day, status) ~ treatment+source+(1|tank),data=oo)
summary(mix_eff)

#mix_eff_2<- coxme(Surv(day, status) ~ source+(1|treatment)+(1|tank),data=oo)
#summary(mix_eff_2)
#anova(simp_model_all,mix_eff, test = 'Chisq')

#Random effects
#Group Variable Std Dev Variance   r.brood Intercept 0.0199992574 0.0003999703   r.worker Intercept 0.0199992574 0.0003999703

#to see if removal is overall dependend upon treatment i use the package multcomp as following with the second warning message:
library(multcomp )
coxme.glht<-glht(mix_eff,linfct=mcp(source="Tukey"))
summary(coxme.glht,test=Chisqtest())

#coxme.glht_source<-glht(mix_eff_2,linfct=mcp(source="Tukey"))
#summary(coxme.glht_source,test=Chisqtest())

coxme.source.glht<-glht(all_rand,linfct=mcp(source="Tukey"))
summary(coxme.source.glht,test=Chisqtest())
#Is it the effect sizes?

fit_oo <- survfit(formula = Surv(day, status) ~ treatment+source,data=oo, type = "kaplan-meier", 
        conf.type = "log")
summary(fit_oo)


#Kaplan estimator for each sourch separetly 
#now, look first only on the planulae from the experiment 

#check difference only in aquaria
d_aquaria <-  oo%>% 
  filter(source == "aquaria")

fit_aquaria <- survfit(formula = Surv(day, status) ~ treatment,data=d_aquaria, type = "kaplan-meier", 
                  conf.type = "log")
#summary of survival in the aquaria 
summary(fit_aquaria)
#check of significance using KM test
s_diff_aquaria <- survdiff(Surv(day, status) ~ treatment, data=d_aquaria)
1 - pchisq(s_diff_aquaria $chisq, length(s_diff_aquaria $n) - 1)

res <- pairwise_survdiff(Surv(day, status) ~ treatment,
                         data = d_aquaria)
res
#now, look first only on the planulae from the experiment 

#check difference only in aquaria
d_sea <-  oo%>% 
  filter(source == "sea")

fit_sea <- survfit(formula = Surv(day, status) ~ treatment,data=d_sea, type = "kaplan-meier", 
                       conf.type = "log")
summary(fit_sea)
s_diff_sea<- survdiff(Surv(day, status) ~ treatment, data=d_sea)
1 - pchisq(s_diff_sea $chisq, length(s_diff_sea $n) - 1)

res <- pairwise_survdiff(Surv(day, status) ~ treatment,
                         data = d_sea)
res


###from here i dont think i need anymore 

g_simple <- ggadjustedcurves(m1,data=d_aquaria,method = "conditional", variable = "treatment",
                 fun = "events" , size= 2, 
                 palette = "jco",legend = "bottom",
                 legend.title = "Treatment",xlab="Day post release",
                 legend.labs = c("Control",
                                 "RCP 4.5", "RCP 8.5"))+
  theme_bw() %+replace% 
  theme(legend.position = "none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor.y=element_line(colour = "gray"),
        panel.background = element_blank(),
        axis.title =element_text(size=10),
        axis.text=element_text(size=10),
        axis.text.x = element_text(size=10))

g_simple


# model selection 
#simple one

m1 <- coxph(Surv(day, status) ~ treatment,data=d_aquaria)
m1          
res.cox <- coxph(Surv(day, status) ~ treatment,data=oo)
ggcoxfunctional(res.cox,  data = oo, point.col = "blue", point.alpha = 0.5)

# random effect model
m2 <- coxme(Surv(day, status) ~ treatment  + (1 | tank),data=d_aquaria)
m2

#nested model - not taken, not better log liklihood or AIC with 1 higher df
m_nest <-  coxme(Surv(day, status) ~ treatment  + (1 | tank/treatment),data=d_aquaria)
print(m_nest)

#nested2 model - not taken - here AIC is uch higher than m2. Although more varience
#is explained by the nested the AIC is too high.. also, p in high = 0.059 - slightly not sign.
m_nest2 <-  coxme(Surv(day, status) ~ treatment  + (1 | falcon/colony_id),data=d_aquaria)
print(m_nest2)


##### embryos derived from the sea 

d_sea <-  oo%>% 
  filter(source == "sea")

#d$day <- as.factor(as.character(d$day))

view(d_sea)

###using cox test with random effect for falcon

#the null model 
m10 <- coxph(Surv(day, status) ~ treatment+falcon,data=d_sea)
m10
#ggforest(m1)

g_simple_sea <- ggadjustedcurves(m10,data=d_sea,method = "conditional", variable = "treatment",
                             fun = "events" , size= 2, 
                             palette = "jco",legend = "bottom",
                             legend.title = "Treatment",xlab="Day post release")+
  theme_bw() %+replace% 
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), panel.grid.minor.y=element_line(colour = "gray"),
    panel.background = element_blank(),
    axis.title =element_text(size=10),
    axis.text=element_text(size=10),
    axis.text.x = element_text(size=10))

g_simple_sea


# model selection 
#simple one

m10 <- coxph(Surv(day, status) ~ treatment,data=d_sea)
m10          
# random effect model
m20 <- coxme(Surv(day, status) ~ treatment  + (1 | tank),data=d_sea)
m20

#nested model - not taken, not better log liklihood or AIC with 1 higher df
m_nest_sea <-  coxme(Surv(day, status) ~ treatment  + (1 | tank/treatment),data=d_sea)
print(m_nest_sea)

#nested2 model - not taken - here AIC is uch higher than m2. Although more varience
#is explained by the nested the AIC is too high.. also, p in high = 0.059 - slightly not sign.
#m_nest2_sea <-  coxme(Surv(day, status) ~ treatment  + (1 | falcon/colony_id),data=d_sea)
#print(m_nest2)

##sohw these plots - i need to make similar y axis
library(scales)
#scale_y_continuous(labels = scales::percent)+
gg_sea <- g_simple_sea+ scale_y_continuous(limits  = c(0, 1), breaks = seq(0, 1, by = 0.2)) 
gg_aquaria <- g_simple+ scale_y_continuous(limits  = c(0, 1), breaks = seq(0, 1, by = 0.2))

library(cowplot)
Surv_together <- plot_grid(gg_sea, gg_aquaria, labels = c('A', 'B'), label_size = 12)
Surv_together

ggsave(plot= Surv_together, "combined survival.tiff", height = 8, width = 16, units = "cm", dpi = 600)

##statisti between survial of A and B - here i using a simple model with association between 
#the source pf the embyos-sea/ aquaria- and  
str(oo)
#the simpe model - 

simp<- coxph(Surv(day, status) ~ source,data=oo)
 summary(simp)

 #model with association 
simp_model_all <- coxph(Surv(day, status) ~ source:treatment,data=oo)
summary(simp_model_all) 

##model with all the data and random effect
all_rand <- coxme(Surv(day, status) ~ treatment+source+(1|tank),data=oo)
summary(all_rand)
stem(exp(ranef(all_rand[[1]])))
         
anova(simp_model_all, test = 'Chisq')

#Random effects
#Group Variable Std Dev Variance   r.brood Intercept 0.0199992574 0.0003999703   r.worker Intercept 0.0199992574 0.0003999703

#to see if removal is overall dependend upon treatment i use the package multcomp as following with the second warning message:
library(multcomp )
 coxme.glht<-glht(all_rand,linfct=mcp(treatment="Tukey"))
 summary(coxme.glht,test=Chisqtest())

 coxme.source.glht<-glht(all_rand,linfct=mcp(source="Tukey"))
 summary(coxme.source.glht,test=Chisqtest())
 #Is it the effect sizes?

#What is diffeent is the treatments RCP 4.5 and 8.5 of theembryos derived from the system 

##risk table for my graphs to show N.                