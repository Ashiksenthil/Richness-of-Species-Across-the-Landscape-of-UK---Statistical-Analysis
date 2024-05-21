if(!is.null(dev.list())) dev.off()  # clear out the past 
rm(list = ls())
cat("\014")
install.packages("DescTools")
install.packages("reshape2")
install.packages("moments")
install.packages("dplyr")
install.packages("tidyr")
library(DescTools)
library(dplyr) # enables pipes %>% and more
library(tidyr) # for spliting on the period see below
library(moments)
library(reshape2)
par(mfrow=c(1, 1)) 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# CHANGE TO THE DIRECTORY WHERE YOU PUT THE CSV FILE
setwd("D:/R/R codes/MA334- DA AND STATS/MA 334 - project/CSV Data")
Proj_data_all_11 <-  read.csv("proportional_species_richness_NAs_removed.csv")

# ASSAIGNING THE 5 GROUPS ASSIGNED  
eco_selected_names <- c("Bird","Butterflies","Carabids","Grasshoppers_._Crickets","Vascular_plants")

#Calculating biodiversity measure for as mean for the 5 group given
mean_selected <- rowMeans(Proj_data_all_11[,eco_selected_names]) # mean the 5 columns 
Proj_data <- Proj_data_all_11%>%select("Location",eco_selected_names,
                                       "Easting","Northing","dominantLandClass",      
                                       "ecologicalStatus","period")%>%
  mutate(eco_status_5=mean_selected)
names(Proj_data) 

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Proj_data$period <- as.factor(Proj_data$period) # must set categorical vars
Proj_data$dominantLandClass <- as.factor(Proj_data$dominantLandClass)

# Selecting the dominantclass of Region England
Proj_data<- Proj_data%>%filter(grepl("e",dominantLandClass))
View(Proj_data)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Creating a Winzorized Mean for the selected data
Winzorized_Mean <- function(var1){
  x <- sort(var1)
  len<-floor(0.2*length(x)) 
  L<-length(x)
  for (i in 1:len){
    x[i]<-x[len+1]
    }
  for (i in ((L-len+1):L)){
    x[i]<- x[L-len]
    }
 return(mean(x))
}


# this creates a table of univariate stats  
table <- data.frame()
for(i in 2:6){
  table <- rbind(table,c(names(Proj_data)[i],
                   min(Proj_data[,i],na.rm = TRUE),
                   quantile(Proj_data[,i],probs = 0.25),
                   round(mean(Proj_data[,i],na.rm = TRUE),digits = 2),
                   median(Proj_data[,i],na.rm = TRUE),
                   quantile(Proj_data[,i],probs = 0.75),
                   max(Proj_data[,i],na.rm = TRUE),
                   round(Winzorized_Mean(Proj_data[,i]),digits = 2)))
}

colnames(table) <- c("taxi_group","min","1st Q","mean","median","3rd Q","Max","Winsorized Mean")
View(table)
write.csv(table,"Univariate Tables")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Extend data exploration with correlations between continuous variables
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

names(Proj_data)
cor_values <- Proj_data%>%select(c(2:6)) 
names(cor_values)
#creating a correlation table.
correlation_table <- round(x = cor(cor_values,use="pairwise.complete.obs"), digits = 2)
View(correlation_table)
write.csv(correlation_table,"Correlation Table") 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# boxplot for only one variable in BD5
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plot(Proj_data$Butterflies~Proj_data$period, xlab = "Eco Period", ylab = "Butterflies")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#hypothesis 1 - KOLOMOGROV-SMIROV TEST
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# comparing the two by CULMULATIVE DISTRIBUTIONS  
par(mfrow=c(1, 1))  # divide graph area in 1 columns
qqplot(Proj_data$eco_status_5,Proj_data$ecologicalStatus)
abline(0,1,col="red")
# both cdfs together  and do a kolmogorov test H0: distributions are the same
BD5_cdf <- ecdf(Proj_data$eco_status_5)
BD11_cdf <- ecdf(Proj_data$ecologicalStatus)
plot(BD11_cdf,col="red") #plotting both the cdfs to compare.
lines(BD5_cdf,col="green")
ks.test(Proj_data$eco_status_5,Proj_data$ecologicalStatus)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#hypothesis2 T-TEST HYPOTHESIS 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

names(Proj_data)
Proj_data_split <- Proj_data%>%select(Location,period,eco_status_5)%>%
  pivot_wider(names_from =period,values_from=eco_status_5)%>%
  mutate(BD5_change=Y00-Y70)
View(Proj_data_split)
hist(Proj_data_split$BD5_change)  # the distribution of the BD5 change 
BD5_change <- Proj_data_split%>%pull(BD5_change)
t.test(BD5_change,mu=0)  # t test with H0: mu=0

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# here inner join  the dataframes for BD5 and BD11

Proj_data_split <- Proj_data%>%select(Location,period,eco_status_5)%>%
  pivot_wider(names_from =period,values_from=eco_status_5)%>%
  mutate(BD5_change=Y00-Y70)

Proj_data_all_11_split <- Proj_data%>%select(Location,period,ecologicalStatus)%>%
  pivot_wider(names_from =period,values_from=ecologicalStatus)%>%
  mutate(BD11_change=Y00-Y70)

Eco_change_BD11 <- Proj_data_all_11_split%>%select(Location,BD11_change)
Eco_change_BD5 <- Proj_data_split%>%select(Location,BD5_change)
Both_eco_change <- inner_join(Eco_change_BD11,Eco_change_BD5,by="Location")
View(Both_eco_change)

# here add two columns for BD5up and BD11up 
Both_eco_change <- Both_eco_change%>%
  mutate(BD11up=ifelse(Both_eco_change$BD11_change>0,1,0))%>%
  mutate(BD5up=ifelse(Both_eco_change$BD5_change>0,1,0))
table(Both_eco_change$BD11up)  # distribution of BD11up
table(Both_eco_change$BD5up)   # distribution of BD5up

# now the joint distribution , a contingency table to interpret 
Table_up_down <- table(Both_eco_change$BD11up,Both_eco_change$BD5up) # contingency table for interpretation 
colnames(Table_up_down) <- c("down","up");rownames(Table_up_down) <- c("down","up")
Table_up_down
View(Table_up_down)
GTest(Table_up_down) # log likelihood ratio test 

#Independent Table
Independent_table <- round(outer(rowSums(Table_up_down),colSums(Table_up_down))/sum(Table_up_down))
summary(Table_up_down) # summary also gives the chi squared test (similar p value)
Independent_table

#odds Ratio 
(OddsRatio(Table_up_down))
#or below method can be used.
p_bd11_up <- Table_up_down[1, 1]/Table_up_down[1, 2]
p_bd5_up <- Table_up_down[2, 1]/ Table_up_down[2,2]
(odd_ratio <- p_bd11_up/ p_bd5_up)

#Sensitivity
(sensitivity <- Table_up_down[1, 1]/colSums(Table_up_down)[1])
#Specificity
(specificity <- Table_up_down[2, 2]/colSums(Table_up_down)[2])
#youden's index
(youden_index <- sensitivity + specificity - 1)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SIMPLE LINEAR REGRESSION
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

par(mfrow = c(1,1))
 
# filtering the data for my dominant land class to run a linear model
BD1 <- Proj_data_all_11 %>% filter(grepl("e",dominantLandClass))

# ploting for the linear model
plot(BD1$Isopods~Proj_data$eco_status_5)
abline(0,1,col="red")

lin_mod <- lm(BD1$Isopods~Proj_data$eco_status_5)
summary(lin_mod)#model summary
abline(lin_mod,col="green")

# some diagnostics can be done here. 
plot(jitter(fitted(lin_mod)),residuals(lin_mod),xlab="Fitted",ylab="Residuals")
abline(h=0,col="blue")

qqnorm(residuals(lin_mod))
qqline(residuals(lin_mod),col="red")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# INITIAL MULTIPLE LINEAR MODEL
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# linear model for the five vaaibles 
lmMod <- lm(BD1$Isopods~.,
            data=Proj_data[c(eco_selected_names)],y=TRUE)
summary (lmMod)  # model summary
cor(lmMod$fitted.values,lmMod$y) # corelation with the data 
AIC(lmMod)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# REDUCED MULTIPLE LINEAR MODEL FROM THE INITIAL MODEL
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# removing butterflies from the initial model.
lmMod_reduced <- lm(BD1$Isopods~.,
                    data=Proj_data[c("Bird","Carabids","Grasshoppers_._Crickets","Vascular_plants")],y=TRUE)
summary(lmMod_reduced)
AIC(lmMod_reduced,lmMod) # here lmMod is preferred by p and AIC criteria

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# INTERACTION MULTIPLE LINEAR MODEL
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# interacting Vascular_plants with Grasshoppers_._Crickets
lmMod_interaction <- lm(BD1$Isopods~
                          Bird+Butterflies+Carabids+Grasshoppers_._Crickets+Vascular_plants
                        +Vascular_plants*Grasshoppers_._Crickets,   
                        data=Proj_data,y=TRUE)
summary(lmMod_interaction )
AIC(lmMod,lmMod_reduced,lmMod_interaction) # model with interataction prefered 
cor(lmMod_interaction$fitted.values,lmMod_interaction$y)# corelation slightly improved  
interaction_data <- BD1$Isopods-lmMod_interaction$fitted.values
plot(interaction_data~lmMod_interaction$fitted.values) # look for unwanted pattern in residuals
abline(0,0,col="red")
qqnorm(interaction_data) # check for normality of residuals in prediction
qqline(interaction_data,col="red")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ONE PERIOD AS TRAINING SET AND ANOTHER ONE AS TEST SET
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# now use one period as the training set and one as the test set 

table(Proj_data$period)
nrow(Proj_data)
Proj_data_Y70 <- Proj_data_all_11%>%filter(period=="Y70") # training set
Proj_data_Y00 <- Proj_data_all_11%>%filter(period=="Y00") # test set
nrow(Proj_data_Y00);nrow(Proj_data_Y00)

lmMod_70 <- lm(Proj_data_Y70$Isopods~.,
               data=Proj_data_Y70[c(eco_selected_names)],y=TRUE)
summary(lmMod_70)
qqnorm(lmMod_70$residuals);qqline(lmMod_70$residuals,col="red")

Predict_00 <- predict(lmMod_70,Proj_data_Y00)

mean((Proj_data_Y70$Isopods-lmMod_70$fitted.values)^2)  # MSE on train data set 
mean((Proj_data_Y00$Isopods-Predict_00)^2)  # MSE on test data (higher)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# OPEN ANALYSIS 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

analysis_data <- Proj_data_all_11%>%select("Location",eco_selected_names,
                                       "Easting","Northing","dominantLandClass",      
                                       "ecologicalStatus","period")%>%
  mutate(eco_status_5=mean_selected)

analysis_data <- analysis_data %>% filter((grepl('SK', Location)|grepl('SP', Location)|grepl('SJ', Location)|grepl('SO', Location))&(grepl("e",dominantLandClass)))
analysis_data <- analysis_data %>% filter(period=='Y00')
analysis_outcome <- lm(eco_status_5
                    ~ dominantLandClass, data = analysis_data, y=T)
summary (analysis_outcome) # Summary of the Linear Model. cor(open_analysis$fitted.values,open_analysis$y)

cor(analysis_outcome$fitted.values,analysis_outcome$y)

qqnorm(analysis_outcome$residuals)
qqline(analysis_outcome$residuals, col= 'red')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++