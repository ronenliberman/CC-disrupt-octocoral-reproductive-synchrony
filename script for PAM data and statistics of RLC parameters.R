library(tidyverse)
library(lubridate) # for working with dates
library(ggplot2)  # for creating graphs
install.packages("phytotools")
library(phytotools)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
rm(list=ls())
d

#calculating RLC parameters using phytotools (pratt 1980)
d <- read_csv('my_Retr.csv')

rlc.data <- d[c(2,3,1)]
rlc.data$Colony <- as.factor(rlc.data$Colony)

rlc.data<- setNames(rlc.data, c("par","etr", "id"))
head(rlc.data)

rlc.data$etr <- na_if(rlc.data$etr, 0)

#rlc.data <- data.frame(stringsAsFactors=FALSE)

ncurves <- length(unique(rlc.data$id)) # number of unique ids in the data 
ids <- unique(rlc.data$id) # store the unique ids 

str(rlc.data)

rlc.parameters <- data.frame(
  id = ids, 
  alpha = 0, 
  beta = 0, 
  ETRmax = 0, 
  Ek = 0, 
  ps = 0
)


for (i in 1:ncurves){
  
  temp.id = ids[i] # extract the id of the curve to be fitted
  
  print(paste("Now fitting curve ", as.character(temp.id))) # to keep track what's happening if the data has many curves
  
  temp.rlc.data <- rlc.data[rlc.data$id==temp.id,] # extract the the data of a single curve into a temporary variable
  PAR = temp.rlc.data$par 
  ETR = temp.rlc.data$etr
  
  fit = fitPGH(PAR, ETR, fitmethod = "Port") # for more options and explanation see package phytotools manual
  
  # store the fitted RLC values into temporary variables
  alpha.rlc = fit$alpha[1]
  beta.rlc = fit$beta[1]
  ps.rlc = fit$ps[1]
  
  # store the parameters
  rlc.parameters$id[i] <- temp.id
  rlc.parameters$alpha[i] <- alpha.rlc
  rlc.parameters$beta[i] <- beta.rlc
  rlc.parameters$ps[i] <- ps.rlc
  
  # calculate ETRmax and Ek for the PGH model (see e.g.Ralph & Gademann 2005 Aquatic Botany 82 (3): 222 - 237). 
  # Note that the equation depends on the model fitted, the code below applies only to the PGH model! 
  # Model equations are documented in the phytotools package code examples (and in the original papers): https://cran.r-project.org/web/packages/phytotools/phytotools.pdf
  
  ETRmax = ps.rlc*(alpha.rlc/(alpha.rlc + beta.rlc))*(beta.rlc/(alpha.rlc+beta.rlc))^(beta.rlc/alpha.rlc)
  Ek = ETRmax/alpha.rlc 
  
  # store the variables
  rlc.parameters$ETRmax[i] <- ETRmax
  rlc.parameters$Ek[i] <- Ek
  
  #plotting the curve and fitted model into a tiff file. By default the file name is the id of the curve. 
  tiff(file=paste0(temp.id, ".tiff"), compression="lzw")
  
  # plot the data, 
  plot(x=PAR, y=ETR, main=temp.id) 
  
  # plot the model fit
  with(fit, {
    P <- ps.rlc*(1-exp(-1*alpha.rlc*PAR/ps.rlc))*exp(-1*beta.rlc*PAR/ps.rlc) # the PGH model equation
    lines(PAR,P)
  }
  ) # end of with
  dev.off() #close the plotting devide. if this is not done, the next run of the loop will override the plot. 
  
}

# now the data frame rlc.parameters contains the fitted values for each curve. Tiff plots should be in current working directory. 
rlc.parameters
View(rlc.parameters)


warnings()

rlc.parameters_first <- rlc.parameters
View(rlc.parameters_first)

write.csv(file="rlc.parameters.csv",rlc.parameters)


####check pam equations difference
library(lme4)
library(nlme)
library(lmerTest)
library(MASS)
library(car)
install.packages("predictmeans")
library(predictmeans)

rm(rlc.parameters)
#FvFm 
d_PAM <- read.csv(file="all PAM_data.csv")
View(d_PAM)
FVFM <- subset(d_PAM, d_PAM$PAR=='1')
View(FVFM)
names(FVFM)[8]<-paste("FVFM")

###manualy copy the FVFM to the rlc parameter data  

##manually calculate maxNPQ for each colony 

np <- read.csv(file = "all PAM_data.csv")
names(np)[1]<-paste("id")
head(np)
str(np)


npq_max <- ddply(np, .(id,treatment), summarise,
      max = max(NPQ, na.rm = F))
view(npq_max)

d <- read_csv('rlc.parameters_NEW.csv')

head(d)
View(d)
d$treatment <- as.factor(d$treatment)
d$Aquaria <- as.factor(d$Aquaria)
d$id <- as.factor(d$id)

etrMAX <- lmer(ETRmax~treatment + (1|Aquaria), data=d)
summary(etrMAX)
anova(etrMAX)

summary(glht(etrMAX, linfct=mcp(treatment="Tukey")))

#normality assumptions
shapiro.test(d$ETRmax)
leveneTest(ETRmax~treatment,d=d)

###alpha
Alpha <- lmer(alpha~treatment + (1|Aquaria), data=d)
anova(Alpha)

summary(glht(Alpha, linfct=mcp(treatment="Tukey")))

#normality assumptions
shapiro.test(d$alpha)
leveneTest(alpha~treatment,d=d)

###eK
head(d)
eK <- lmer(Ek~treatment + (1|Aquaria), data=d)
summary(eK)
anova(eK)

summary(glht(eK, linfct=mcp(treatment="Tukey")))

#normality assumptions
shapiro.test(d$Ek)
leveneTest(Ek~treatment,d=d)

###FV/FM
head(d)
fv <- lmer(FvFm~treatment + (1|Aquaria), data=d)
anova(fv)

summary(glht(fv, linfct=mcp(treatment="Tukey")))

#normality assumptions
shapiro.test(d$FvFm)
leveneTest(FvFm~treatment,d=d)

#npqMAX
head(d)
NPQ<- lmer(npqMAX~treatment + (1|Aquaria), data=d)
anova(NPQ)

summary(glht(NPQ, linfct=mcp(treatment="Tukey")))

#normality assumptions
shapiro.test(d$npqMAX) 
leveneTest(npqMAX~treatment,d=d)

## summary of results 
library(plotrix)
d <- read_csv('rlc.parameters_NEW.csv')
Sum_all <- d %>% 
  group_by(treatment) %>% 
  summarise_each(funs(mean(., na.rm=T), n = sum(!is.na(.)),
                      se = sd(., na.rm=T)/sqrt(sum(!is.na(.)))), alpha:npqMAX)

View(Sum_all)

###graphics 

d <- read.csv("all PAM_data.csv")
head(d)
names(d)[8]<-paste("effective")
#d$treatment[grepl('ambient', d$treatment)] <- 'Ambient'
#View(d)

p_rETR2 <-  ggplot(d, aes(y=rETR*0.5, x=PAR, colour=treatment,group=treatment)) + 
  #geom_point(size=2,alpha=0.5)+
  geom_smooth(size=2,alpha=0.5,na.rm = FALSE,level=0.95)+
  #geom_errorbar(data=Sum_all, aes(ymin=ETRmax_mean-ETRmax_se, ymax=ETRmax_mean+ETRmax_se)
  scale_color_jco() +
  theme_bw() %+replace% 
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor.y=element_blank(), 
        panel.background = element_blank(),
        axis.title =element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(size=12))+
  #scale_x_continuous("PAR", labels = as.numeric(PAR), breaks = PAR)+
  xlab(bquote('�mol photons '~ m^-2~s^-1*''))+ylab('Relative ETR')
# labs(x=parse(text='Âµmol photons, m^-2s^-1'),y=parse(text='rETR'))+
#guides(fill=guide_legend(title="Treatments"))
p_rETR2

p_effective <-  ggplot(d, aes(y=effective, x=PAR, colour=treatment,group=treatment)) + 
  #geom_line(size=0.5)+#geom_point(aes(shape=treatment, color=treatment))+
  #geom_errorbar(aes(ymin = mean_eff-se, ymax= mean_eff+se),  
  #width=.5, position=position_dodge(0.05))+
  #geom_point(aes(colour= treatment),size=1,alpha=0.5)+
  geom_smooth(size=2,alpha=0.5,na.rm = FALSE,level=0.95)+
  scale_color_jco() +
  theme_bw() %+replace% 
  theme(legend.position="none",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor.y=element_blank(), 
        panel.background = element_blank(),
        axis.title =element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(size=12))+
  #xlab(bquote('Âµmol photons' ('m^-2 s^-1')))
  xlab(bquote('�mol photons '~ m^-2~s^-1*''))+ylab('Effective quantum yield')
#labs(x=expression(paste("Âµmol photons ", m^{-2})),y='Effective quantum yield')
p_effective

p_NPQ2 <-  ggplot(d, aes(y=NPQ, x=PAR, colour=treatment,group=treatment)) + 
  #  geom_line(size=0.5)+#geom_point(aes(shape=treatment, color=treatment))+
  # geom_errorbar(aes(ymin = mean_eff-se, ymax= mean_eff+se),  
  #              width=.5, position=position_dodge(0.05))+
  #geom_point(aes(colour= treatment),size=2,alpha=0.5)+
  geom_smooth(size=2,alpha=0.5,na.rm = FALSE,level=0.95)+
  scale_color_jco() +
  theme_bw() %+replace% 
  theme(legend.position="left",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor.y=element_blank(), 
        panel.background = element_blank(),
        axis.title =element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(size=12))+
  xlab(bquote('�mol photons '~ m^-2~s^-1*''))+ylab('Non-photochemical quenching')+
  # scale_x_continuous(limits  = c(21, 32), breaks = seq(21, 32, by = 1)) 
  #labs(x=expression(paste("Âµmol photons ", m^{-2})),
  #      y='Non-photochemical quenching')+
  # labs(x=parse(text='Âµmol photons, m^-2s^-1'),y=parse(text='rETR'))+
  guides(color=guide_legend(override.aes=list(fill=NA)))
p_NPQ2

ggsave("for legend.tiff", plot = p_NPQ2, width = 8.5 , height = 8 ,dpi=600, 
       units = "cm")

one_grid <- plot_grid(p_rETR2, p_effective, p_NPQ2, labels = c('A', 'B', 'C'),
                      # A negative rel_height shrinks space between elements
                      label_x = 0.15, label_y = 0.985, label_size = 14,ncol = 1)
one_grid
ggsave("RLC_June2020.tiff", plot = one_grid, width = 8.5 , height = 21 ,dpi=600, 
       units = "cm")
