library(glmmTMB)
library(splitstackshape)
library(ggplot)
library(scales)
library(ggsci)
set <- read.csv('advanced_merged.csv')
View(set)
set <- set[c(1:4)]

###making my data only 1'ns
sett<- expandRows(set, "count")
View(sett)


 ggplot(set,aes(x=count, fill= stage)) +
   geom_histogram() + facet_wrap(~time)

  
 plot_data <- set %>% 
   count(SizeOfTumor, Age, RemovalSurgery) %>% 
   group_by(Age, SizeOfTumor) %>% 
   mutate(percent = n/sum(n))
 table( set$stage)
 
# set$stage <- ordered( set$stage, levels = c( "Pre-metamorphosed", "Early-metamorphosed ", "Advanced-metamorphosed1", "Advanced-metamorphosed2"))
 levels(set$time)
set$time <- factor(set$time, labels= c("10 days polyp", "20 days polyp"))
 head(set)
 
 p <- ggplot(set, aes(x = treatment, y = count, fill = stage)) + 
   geom_col(position = "fill") + 
   scale_fill_d3()+
  # geom_label(aes(label = percent(percent)), position = "fill", color = "white", vjust = 1, show.legend = FALSE) +
   scale_y_continuous(labels = percent) +
   facet_grid(.~time)+  theme_classic2() %+replace% 
   theme(legend.position="bottom",
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), panel.grid.minor.y=element_blank(), 
         panel.background = element_blank(),
         axis.title =element_text(size=12),
         axis.text=element_text(size=12),
         axis.text.x = element_text(size=10))+
   #scale_x_continuous("PAR", labels = as.numeric(PAR), breaks = PAR)+
   xlab('')+ylab('Development of primary polyp')+ labs(fill= "")
 p
 
 ggsave("development polyp.tiff", plot = p, width = 17 , height = 8.5 ,dpi=600, 
        units = "cm")
 
 hist(set$count)
 
 glmm.fit1=glm(count~ 0+time*stage+treatment ,family="poisson", data = set) 
 summary(glmm.fit1)
Anova(glmm.fit1)

glht(glmm.fit1, mcp(treatment="Tukey"))
 
plot(glmm.fit1$residuals)

### using ordinal logistic regression 
set <- read.csv('Tachles_LEVELS.csv')
#View(set)
set <- set[c(1:4)]
#head(set)

###making my data only 1'ns
sett<- expandRows(set, "count")
View(sett)
#Ordering the dependent variable
sett$stage = factor(sett$stage, levels = c("Early-metamorphosed", "Advanced-metamorphosed1", "Advanced-metamorphosed2"), ordered = TRUE) 
sett$time = factor(sett$time, levels = c("1", "2"), ordered = FALSE) 
sett$treatment = factor(sett$treatment, levels = c("Ambient ", "RCP 4.5","RCP 8.5"), ordered = FALSE)  

head(sett)

str(sett)
#Exploratory data analysis 
#Summarizing the data
summary(sett)
#Making frequency table
table(sett$stage, sett$time, sett$treatment)

#Dividing data into training and test set
#Random sampling 
samplesize = 0.60*nrow(sett)
set.seed(100)
index = sample(seq_len(nrow(sett)), size = samplesize)
#Creating training and test set 
datatrain = sett[index,]
datatest = sett[-index,]
#Build ordinal logistic regression model
mod= polr(stage ~ treatment*time  , data = datatrain, Hess = TRUE)
summary(mod)

#Compute confusion table and misclassification error
predict_develpment = predict(mod,datatest)
table(datatest$stage, predict_develpment)
mean(as.character(datatest$stage) != as.character(predict_develpment))
#Plotting the effects 
install.packages("effects")
library("effects")

Effect(focal.predictors = "treatment",mod)
plot(Effect(focal.predictors = "time",mod))
plot(Effect(focal.predictors = c("treatment", "time"),mod))


########################another model ##########
ftable(xtabs(~ treatment + time + stage, data =sett))

mod2= polr(stage ~ treatment+time  , data = datatrain, Hess = TRUE)
summary(mod2)

# 
## view a summary of the model
(ctable <- coef(summary(mod2)))
p2 <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p2))

(ci <- confint(mod2)) # default method gives profiled CIs

confint.default(mod2) #### assuming normality - not using this

## odds ratios
exp(coef(mod2))

## OR and CI
exp(cbind(OR = coef(mod2), ci))


install.packages("Hmisc")
library(Hmisc)

sf <- function(y) {
  c('Y>=1' = qlogis(mean(y >= 1)),
    'Y>=2' = qlogis(mean(y >= 2)),
    'Y>=3' = qlogis(mean(y >= 3)))
}

(s <- with(sett, summary(as.numeric(stage) ~ treatment:time, fun=sf)))

glm(I(as.numeric(stage) >= 2) ~ treatment, family="binomial", data = sett)

plot(s, which=1:3, pch=1:3, xlab='logit', main=' ', xlim=range(s[,3:4]))

newdat <- data.frame(
  treatment= rep(1:3, 100),
  time = rep(1:2, each = 150))

newdat <- cbind(newdat, predict(mod2, newdat, type = "probs"))


########################
### Visualizing the data with stack plot 
d <- read.csv('settelment data for ordinal.csv')
View(set)
d <- d[c(1:4)]

#head(set)
d$time <- factor(d$time, labels= c("18 days polyp", "28 days polyp"))
head(sett)
levels(d$stage)

d$stage <- ordered(d$stage, levels = c( "Early-metamorphosed", "Advanced-metamorphosed", "Fully-metamorphosed"))
p <- ggplot(d, aes(x = treatment, y = count, fill = stage)) + 
  geom_col(position = "fill") + 
  scale_fill_grey()+
  # geom_label(aes(label = percent(percent)), position = "fill", color = "white", vjust = 1, show.legend = FALSE) +
  scale_y_continuous(labels = percent) +
  facet_grid(.~time)+  theme_classic2() %+replace% 
  theme(legend.position="bottom",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor.y=element_blank(), 
        panel.background = element_blank(),
        axis.title =element_text(size=12),
        axis.text=element_text(size=12),
        axis.text.x = element_text(size=10))+
  #scale_x_continuous("PAR", labels = as.numeric(PAR), breaks = PAR)+
  xlab('')+ylab('Development of primary polyp')+ labs(fill= "")
p

ggsave("development polyp_June.tiff", plot = p, width = 17 , height = 8.5 ,dpi=600, 
       units = "cm")




## cumulative link model with clm #########
###making my data only 1'ns
sett<- expandRows(d, "count")
# set$stage <- ordered( set$stage, levels = c( "Pre-metamorphosed", "Early-metamorphosed ", "Advanced-metamorphosed1", "Advanced-metamorphosed2"))
levels(sett$time)


install.packages("ordinal")
library(ordinal)

fm1 <- clm(stage ~ treatment * time, data=sett)
summary(fm1)
anova(fm1, type="III")

fm.nom <- clm(stage ~ treatment, nominal = ~ time, data=sett)
summary(fm.nom)
anova(fm1, fm.nom)

#Modelling scale effects
#To allow the scale of the latent variable distribution to depend on explanatory variables
fm.sca <- clm(stage ~ treatment + time, scale = ~ treatment, data=sett)
summary(fm.sca)

#visualizing the model data
pr1 <- profile(fm1, alpha=1e-4)
plot(pr1)

#Assessment of model convergence
#Likelihood slices

slice.fm1 <- slice(fm1, lambda = 5)
par(mfrow = c(2, 3))
plot(slice.fm1)

#Parameter accuracy
convergence(fm1)

