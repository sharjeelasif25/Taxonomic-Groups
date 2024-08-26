if(!is.null(dev.list())) dev.off()  # clear out the past 
rm(list = ls())
cat("\014")
install.packages('pander')

library(DescTools)
library(dplyr) 
library(tidyr) 
library(moments)
library(reshape2)
par(mfrow=c(1, 1)) 
library(DescTools)
library(ggplot2)
library(reshape2)

# Setting the working directory and reading the dataset

setwd("D:\\MSc Work\\Autumn\\Data Analytics with R\\Final Project")
Proj_data_all_11 <-  read.csv("proportional_species_richness_NAs_removed.csv")

# Selecting my 5 taxonomic groups

eco_selected_names <- c("Bird","Bryophytes","Hoverflies","Isopods","Ladybirds")

# Calculating the biodiversity measure as the mean of 5 groups 

mean_of_5_Varaibles <- rowMeans(Proj_data_all_11[,eco_selected_names]) 
Proj_data <- Proj_data_all_11%>%select("Location",eco_selected_names,"Easting","Northing","dominantLandClass","ecologicalStatus","period")%>%mutate(eco_status_5=mean_of_5_Varaibles)
names(Proj_data)



# Setting the categorical varialbes Periods and Dominant land Class

Proj_data$period <- as.factor(Proj_data$period) 
Proj_data$period
Proj_data$dominantLandClass <- as.factor(Proj_data$dominantLandClass)
Proj_data$dominantLandClass

#-------------------------------------------------------------------------------------------------#

#Question 1a & 1b (Univariate analysis)

table <- data.frame()
for(i in c(2:6)){
  table <- rbind(table,
                   c(names(Proj_data)[i],
                   round(mean(Proj_data[,i],na.rm = TRUE),digits = 2),
                   round(sd(Proj_data[,i],na.rm = TRUE),digits = 2),
                   round(median(Proj_data[,i],na.rm = TRUE),digits = 2),
                   round(min(Proj_data[,i],na.rm = TRUE),digits = 2),
                   round(max(Proj_data[,i],na.rm = TRUE),digits = 2),
                   round(mean(Winsorize(Proj_data[,i],probs= c(0.10,0.90))), digits = 2),
                   round(quantile(Proj_data[,i],na.rm = TRUE,probs = 0.25),digits = 2),
                   round(quantile(Proj_data[,i],na.rm = TRUE,probs = 0.75),digits = 2)
))}
colnames(table) <- c("Taxonomic Group","Mean","SD","Median","Min","Max","Winsorized Mean", "First Quantile", 'Third Quantile')
print(table)
pandoc.table(table)
pandoc.table(table, style = "rmarkdown")
formattable(table)

# ------------------------------------------------------------------------------- #

# Question 2 (Correlation matrix)

cont_vars <- Proj_data%>%select(c(2:6))
cormat <- round(x = cor(cont_vars,use="pairwise.complete.obs"), digits = 2)
cormat

get_lower_tri<-function(cormat){ # Getting the lower triangle values
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

get_upper_tri <- function(cormat){ # Getting the lower triangle values
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat) 
melted_cormat <- melt(upper_tri, na.rm = TRUE)

reorder_cormat <- function(cormat){ # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  }

cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
melted_cormat <- melt(upper_tri, na.rm = TRUE) # Melt the correlation matrix

# Create a ggheatmap

ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+ geom_tile(color = "white")+ scale_fill_gradient2(low = "white", high = "darkblue", mid = "lightblue", midpoint = 0, limit = c(-1,1), space = "Lab", name="Correlation\nMatrix") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()

ggheatmap + geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = c(0.5, 0.7),
      legend.direction = "horizontal") + 
      guides(fill = guide_colorbar(barwidth = 6, barheight = 1, title.position = "top", title.hjust = 0.5))

# ----------------------------------------------------------------------------------- #

# Question 3 (Box Plot of Isopods)

eco_period <- Proj_data%>%pull(period)
boxplot(Proj_data$Isopods~eco_period, xlab = 'Ecological Period', ylab = 'Isopods Diversity', col = c('lightgrey','lightblue'), at = c(1,2), names = c('Year 2000','Year 1970'))

# ----------------------------------------------------------------------------------- #

# Hypothesis Test 1 (T - Test)

Proj_data_split <- Proj_data%>%select(Location,period,eco_status_5)%>% pivot_wider(names_from =period,values_from=eco_status_5)%>% mutate(BD5_change=Y00-Y70)
hist(Proj_data_split$BD5_change)
BD5_change <- Proj_data_split%>%pull(BD5_change)
t.test(BD5_change,mu=0)  

# ---------------------------------------------------------------------------------- #

# Hypothesis Test 2 (KS - Test)

BD5_cdf <- ecdf(Proj_data$eco_status_5)
BD11_cdf <- ecdf(Proj_data$ecologicalStatus)
plot(BD11_cdf,col="red")
lines(BD5_cdf,col="green")
ks.test(Proj_data$eco_status_5,Proj_data$ecologicalStatus)

# ---------------------------------------------------------------------------------- #

# Contingency Table

Eco_change_BD11 <- Proj_data_all_11_split%>%select(Location,BD11_change)
Eco_change_BD5 <- Proj_data_split%>%select(Location,BD5_change)
Both_eco_change <- inner_join(Eco_change_BD11,Eco_change_BD5,by="Location")

Both_eco_change <- Both_eco_change%>%                                         # Adding the two columns for BD5up and BD11up
  mutate(BD11up=ifelse(Both_eco_change$BD11_change>0,"Increase","Decrease"))%>%
  mutate(BD5up=ifelse(Both_eco_change$BD5_change>0,"Increase","Decrease"))
table(Both_eco_change$BD11up)  # distribution of BD11up
table(Both_eco_change$BD5up)   # distribution of BD5up

Contingency_Table <- table(Both_eco_change$BD11up,Both_eco_change$BD5up) # contingency table of BD11 and BD5
Contingency_Table

# Test of independence

Chisq_test <- chisq.test(Contingency_Table)
Chisq_test

# likelihood-ratio statistic

likelihood_ratio_stat <- GTest(Contingency_Table) 
likelihood_ratio_stat

# Odds_Ratio 

odds_ratio <- (Contingency_Table[1, 1] * Contingency_Table[2, 2]) / (Contingency_Table[1, 2] * Contingency_Table[2, 1])
odds_ratio

# Sensitivity 

Sensitivity <- Contingency_Table[2, 2] / sum(Contingency_Table[2, ])
Sensitivity

# Specificity 

Specificity <- Contingency_Table[1, 1] / sum(Contingency_Table[1, ])
Specificity

# Youden's Index

Youden_Index <- Sensitivity + Specificity - 1
Youden_Index


# ---------------------------------------------------------------------------------- #

# Simple Linear Regression

plot(Proj_data_all_11$Carabids~Proj_data$eco_status_5)
abline(0,1,col="red")
lin_mod <- lm(Proj_data_all_11$Carabids~Proj_data$eco_status_5)
abline(lin_mod,col="green")
summary(lin_mod)

# Residual Plots 

plot(jitter(fitted(lin_mod)),residuals(lin_mod),xlab="Fitted",ylab="Residuals") 
abline(h=0,col="blue")
qqnorm(residuals(lin_mod))
qqline(residuals(lin_mod),col="red")

# ---------------------------------------------------------------------------------- #

# Multiple Linear Regression

lmMod <- lm(Proj_data_all_11$Carabids~., data=Proj_data[c(eco_selected_names)],y=TRUE)
summary (lmMod)  # model summary
cor(lmMod$fitted.values,lmMod$y) # corelation with the data 
plot(Proj_data_all_11$Carabids~lmMod$fitted.values)
abline(0,1,col="red")

# Residuals for the training model fit to the test data set 

mis_fit_to_Data <- Proj_data_all_11$Carabids-lmMod$fitted.values

# Unwanted patterns in residuals

plot(mis_fit_to_Data~lmMod$fitted.values) 
abline(0,0,col="red")

# check for normality of residuals in prediction

qqnorm(mis_fit_to_Data) 
qqline(mis_fit_to_Data,col="red")

# AIC and Model Summary

AIC(lmMod)
summary(lmMod)

# Creating a reduced model after dropping the insignificant group

lmMod_reduced <- lm(Proj_data_all_11$Carabids~.,data=Proj_data[c("Bryophytes","Hoverflies","Isopods","Ladybirds")],y=TRUE)
summary(lmMod_reduced)
AIC(lmMod_reduced,lmMod) # Reduced model is giving slightly better results and is a better fit

# Creating an interaction model

lmMod_interaction <- lm(Proj_data_all_11$Carabids~
                          Bird+Hoverflies+Isopods+Bryophytes+Ladybirds
                        +Ladybirds*Bryophytes,  # Ladybirds and Bryophytes as interaction terms
                        data=Proj_data,y=TRUE)
summary(lmMod_interaction )
AIC(lmMod,lmMod_reduced,lmMod_interaction) # Interaction model gives the best results and is the best fit
cor(lmMod_interaction$fitted.values,lmMod_interaction$y) # Correlation improved 

# Creating Y00 as test set and Y70 as training set

table(Proj_data$period)
nrow(Proj_data)
Proj_data_Y70 <- Proj_data_all_11%>%filter(period=="Y70") # training set
Proj_data_Y00 <- Proj_data_all_11%>%filter(period=="Y00") # test set
nrow(Proj_data_Y00)
nrow(Proj_data_Y00)

#Running the model on the training set

lmMod_70 <- lm(Proj_data_Y70$Carabids~.,data=Proj_data_Y70[c(eco_selected_names)],y=TRUE)
qqnorm(lmMod_70$residuals);qqline(lmMod_70$residuals,col="red")
plot(lmMod_70$residuals~lmMod_70$fitted.values) # look for unwanted pattern in residuals
abline(0,0,col="red")

#Running the model on the test set

Predict_00 <- predict(lmMod_70,Proj_data_Y00)
plot(Predict_00~Proj_data_Y00$Carabids)
abline(0,1,col="red")

# Getting the mean square error of both data sets

mean((Proj_data_Y70$Carabids-lmMod_70$fitted.values)^2)  # MSE on train data set 
mean((Proj_data_Y00$Carabids-Predict_00)^2)  # MSE on test data (higher)


# ---------------------------------------------------------------------------------- #

# Open Analysis - Comparing the biodiversity change in coastal and mountainous regions of all of UK

# Getting the mountainous dominant land classes in Y70

Mountains70 <- Proj_data %>%
  filter(dominantLandClass == c('18e','22e','23e','18s','19s','21s','22s','23s','24s','17w1','17w2','17w3','18w'), period == "Y70") %>%
  group_by(eco_status_5, dominantLandClass, period) %>% summarize()

# Getting the mountainous dominant land classes in Y700

Mountains00 <- Proj_data %>%
  filter(dominantLandClass ==  c('18e','22e','23e','18s','19s','21s','22s','23s','24s','17w1','17w2','17w3','18w'), period == "Y00") %>%
  group_by(eco_status_5, dominantLandClass, period) %>% summarize()

# Plotting the biodiversity change in mountainous regions as a histogram

Mountains <- rbind(data.frame(period = "Y00", eco_status_5 = Mountains00$eco_status_5),
                   data.frame(period = "Y70", eco_status_5 = Mountains70$eco_status_5))
ggplot(Mountains, aes(x=eco_status_5, fill=period)) +
  geom_histogram(binwidth = 0.01, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("#B1B6BD", "#4D8DDE")) +
  labs(title = "BD5 in Mountaineous Regions", x = "BD5 Group", y = "Occurrences") + theme(plot.title = element_text(hjust = 0.5))

# Getting the coastal dominant land classes in Y70

Coast70 <- Proj_data %>%
  filter(dominantLandClass == c('7e','8e','7s','29s','30s','31s','32s','7w'), period == "Y70") %>%
  group_by(eco_status_5, dominantLandClass, period) %>% summarize()

# Getting the coastal dominant land classes in Y70

Coast00 <- Proj_data %>%
  filter(dominantLandClass ==  c('7e','8e','7s','29s','30s','31s','32s','7w'), period == "Y00") %>%
  group_by(eco_status_5, dominantLandClass, period) %>% summarize()

# Plotting the biodiversity change in coastal regions as a histogram

Coast <- rbind(data.frame(period = "Y00", eco_status_5 = Coast00$eco_status_5),
               data.frame(period = "Y70", eco_status_5 = Coast70$eco_status_5))
ggplot(Coast, aes(x=eco_status_5, fill=period)) +
  geom_histogram(binwidth = 0.01, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = c("#B1B6BD", "#4D8DDE")) +
  labs(title = "BD5 in Costal Regions", x = "BD5 Group", y = "Occurrences") + theme(plot.title = element_text(hjust = 0.5))

# As per the two histograms, biodiversity measure has a more observable change in coastal areas as compared to mountainous regions of the UK
# The probably causes and reasons of this variation in effects is discussed in the pdf file attached alongwith this code.