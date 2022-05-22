## Oculina Heterotrophy Experiment Data Analysis ####
# Author: Hanny Rivera
# Script associated with manuscript Rivera et al. 2022: TITLE 

## Setup ####

#setwd("/Users/hannyrivera/Documents/BU/DaviesLab/_Data/Oculina_Exp/clean_repo")
# Reading in the master data frame with all the raw experimental values 
all_data<-read.table("input_files/Master_data.txt", header=T)
all_data$treatment<-as.factor(all_data$treatment)

library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(ggbiplot)
library(xts)
library(zoo)
library(TTR)
library(scales)
library(signal)
library(stargazer)
library(Rmisc)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(vegan)
library(ggvenn)
library(ggpubr)
library(multcomp)
library(lme4)
library(emmeans)

## Temperature Data #### 
# load temperature file from tank loggers
temp_unfed_hot<-xts(zoo(read.table("input_files/loggers/sys5_unfed_hot.txt", header=TRUE)$temp,seq.POSIXt(ISOdate(2019,8,9,13,0,0),ISOdate(2019,9,9,17,0,0), "30 min", tz="GMT")))
temp_unfed_ctrl<-xts(zoo(read.table("input_files/loggers//sys7_unfed_ctrl.txt", header=TRUE)$temp,seq.POSIXt(ISOdate(2019,8,9,13,0,0),ISOdate(2019,9,9,17,0,0), "30 min", tz="GMT")))
temp_fed_hot<-xts(zoo(read.table("input_files/loggers//sys8_fed_hot.txt", header=TRUE)$temp,seq.POSIXt(ISOdate(2019,8,9,13,0,0),ISOdate(2019,9,9,17,0,0), "30 min", tz="GMT")))
temp_fed_ctrl<-xts(zoo(read.table("input_files/loggers//sys6_fed_ctrl.txt", header=TRUE)$temp,seq.POSIXt(ISOdate(2019,8,9,13,0,0),ISOdate(2019,9,9,17,0,0), "30 min", tz="GMT")))

# Merge series into one object
alltemps<-merge(temp_unfed_hot, temp_unfed_ctrl, temp_fed_hot, temp_fed_ctrl, join="left", fill=NA)

alltemps<-alltemps["2019-08-13 12:00:00/"]



# Name columms
colnames(alltemps)<-c("unfed_hot", "unfed_ctrl", "fed_hot", "fed_ctrl")
# Convert to data frame for ggplot plotting
alltemps_gg<-as.data.frame(alltemps)
# renames rows as increasing numbers
row.names(alltemps_gg)<-as.character(seq(1,length(alltemps_gg[,1])))
# add in date time values as variable
alltemps_gg$datetime<-time(alltemps)
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
alltemps_gg<-melt(alltemps_gg, id.vars = "datetime", measure.vars = c("unfed_hot", "unfed_ctrl", "fed_hot", "fed_ctrl"))
# fix column headers again
colnames(alltemps_gg)<-c("datetime", "treatment", "temp")

# Plot temperature time series
# You'll have to color of the lines for the systems to match the treatment colors from the PAM file. 
# See the metadata file to check what system was what treatment. 

(temp_plot<-ggplot(alltemps_gg,aes(x=datetime, y=temp, colour=treatment))+geom_line(size=.6)+theme_minimal()+
  scale_y_continuous(breaks=seq(23, 33, by=1),limits = c(22.8,33.2))+
  scale_x_datetime(date_labels=c("28", "7", "14", "21"), date_breaks="days")+
  ggtitle("")+ylab("Temperature (°C)")+xlab("")+
  theme_bw()+
  theme(axis.text.x = element_text(size=10, face="bold", colour="black"),
        axis.text.y = element_text(size=10, face="bold", colour="black"),
        axis.title.y = element_text(size=12, face="bold", colour="black"),
        panel.grid.major.y = element_line(size=.75),
        panel.grid.major.x = element_line(size=.75), 
        panel.grid.minor.x = element_line(size=.3), 
        panel.grid.minor.y = element_blank(),
        legend.position = "top", 
        legend.background = element_rect(fill="white", colour="grey60"),
        legend.text = element_text(size=10, face="bold"),
        legend.title = element_text(size=10, face="bold"),
        plot.margin = margin(15,15,0,0,"pt"))+
  scale_color_manual(limits=c("unfed_ctrl", "fed_ctrl", "unfed_hot",  "fed_hot"),
                     values=c("navyblue", "cyan3", "firebrick", "darkgoldenrod2"),
                     labels=c("Unfed, Ambient", "Fed, Ambient", "Unfed, Heated", "Fed, Heated"))+
    labs(colour="Treatment"))

ggsave(filename = "Loggers/FigureX_logger_temps.png", height= 5, width=7, units="in", dpi=300)

## PAM ####
PAM<-read.table("input_files/physio/PAM.txt", header=TRUE)

PAM %>% group_by(timepoint, treatment) %>% summarise(yield=mean(fvfm),ster=sd(fvfm)/sqrt(length(fvfm)))->pam_sum

#add in the temps when they start diverging
colMeans(alltemps["2019-08-27",]) # 25.5 C for hot treatments 
colMeans(alltemps["2019-09-03",]) # 31C for hot treatment 
colMeans(alltemps["2019-09-09",]) # 33C for hot treatment

 
PAM %>% group_by(timepoint, coral_id, treatment) %>% summarise(avg_yield=mean(fvfm))->pam_avg
pam_avg$treatment<-as.factor(pam_avg$treatment)
pam_avg$timepoint<-as.factor(pam_avg$timepoint)

#add in pam data to all_data for PCA later 
# re-organize pam data for wgcna matrix 
pivot_wider(dplyr::select(PAM,coral_id,timepoint,fvfm), names_from =timepoint, values_from =fvfm)->pam_wide
colnames(pam_wide)[2:5]<-c("pam1", "pam2", "pam3", "pam4")
# combine into all_data dataframe
all_data<-left_join(all_data, pam_wide, by="coral_id")


# test for differences in PAM over time. 
pam_lm <- lmer(avg_yield~treatment*timepoint+(1|coral_id), data = pam_avg) 
anova(pam_lm)
emmeans(pam_lm, list(pairwise ~ treatment|timepoint), adjust = "tukey")
# $`pairwise differences of treatment | timepoint`
# timepoint = 1:
#   2                      estimate     SE  df t.ratio p.value
# fed_ctrl - fed_hot      0.04835 0.0400 172   1.208  0.6228
# fed_ctrl - unfed_ctrl   0.01763 0.0400 172   0.440  0.9714
# fed_ctrl - unfed_hot    0.00885 0.0393 172   0.225  0.9960
# fed_hot - unfed_ctrl   -0.03072 0.0400 172  -0.767  0.8691
# fed_hot - unfed_hot    -0.03950 0.0393 172  -1.005  0.7468
# unfed_ctrl - unfed_hot -0.00878 0.0393 172  -0.223  0.9961
# 
# timepoint = 2:
#   2                      estimate     SE  df t.ratio p.value
# fed_ctrl - fed_hot      0.01303 0.0400 172   0.325  0.9881
# fed_ctrl - unfed_ctrl   0.01165 0.0400 172   0.291  0.9914
# fed_ctrl - unfed_hot    0.01785 0.0393 172   0.454  0.9687
# fed_hot - unfed_ctrl   -0.00137 0.0400 172  -0.034  1.0000
# fed_hot - unfed_hot     0.00483 0.0393 172   0.123  0.9993
# unfed_ctrl - unfed_hot  0.00620 0.0393 172   0.158  0.9986
# 
# timepoint = 3:
#   2                      estimate     SE  df t.ratio p.value
# fed_ctrl - fed_hot      0.14990 0.0400 172   3.744  0.0014
# fed_ctrl - unfed_ctrl  -0.00440 0.0400 172  -0.110  0.9995
# fed_ctrl - unfed_hot    0.16063 0.0393 172   4.086  0.0004
# fed_hot - unfed_ctrl   -0.15429 0.0400 172  -3.854  0.0009
# fed_hot - unfed_hot     0.01073 0.0393 172   0.273  0.9929
# unfed_ctrl - unfed_hot  0.16502 0.0393 172   4.198  0.0003
# 
# timepoint = 4:
#   2                      estimate     SE  df t.ratio p.value
# fed_ctrl - fed_hot      0.27947 0.0400 172   6.981  <.0001
# fed_ctrl - unfed_ctrl   0.02054 0.0400 172   0.513  0.9559
# fed_ctrl - unfed_hot    0.25576 0.0393 172   6.506  <.0001
# fed_hot - unfed_ctrl   -0.25894 0.0400 172  -6.468  <.0001
# fed_hot - unfed_hot    -0.02372 0.0393 172  -0.603  0.9309
# unfed_ctrl - unfed_hot  0.23522 0.0393 172   5.983  <.0001
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 

(a<-ggplot(pam_sum, aes(x=timepoint, y=yield, color=treatment))+
  geom_line(position=position_dodge(width=0.2))+
  geom_errorbar(aes(ymin=yield-ster, ymax=yield+ster), position=position_dodge(width=0.2), width=0)+
  ylab(expression(Photosynthetic~Yield~(f[v]/f[m])))+
  xlab("Experiment Day")+
  scale_x_discrete(labels=c("9", "14", "22", "28"), 
                   limits=c(1,2,3,4))+
  scale_colour_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                    values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                    labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
  labs(colour="Treatment")+theme_bw()+
    theme(axis.text.x = element_text(size=10, face="bold", colour="black"),
          axis.text.y = element_text(size=10, face="bold", colour="black"),
          axis.title.y = element_text(size=12, face="bold", colour="black"),
          axis.title.x = element_text(size=12, face="bold", colour="black"),
          panel.grid.major.y = element_line(size=.75),
          panel.grid.major.x = element_line(size=.75), 
          panel.grid.minor.x = element_line(size=.3), 
          panel.grid.minor.y = element_blank(),
          legend.position = "top", 
          legend.background = element_rect(fill="white", colour="grey60"),
          legend.text = element_text(size=10, face="bold"),
          legend.title = element_text(size=10, face="bold"))+
    annotate("text", x=0.9, y=0.65, label="italic(PAM)", size=5, parse=T)+
    annotate("text", x=2, y=0.53, label="25.5°C", size=4)+
    annotate("text", x=3, y=0.39, label="31°C", size=4)+
    annotate("text", x=4, y=0.26, label="33°C", size=4)+
    annotate("text", x=3, y=0.54, label="B", size=4)+
    annotate("text", x=3, y=0.66, label="A", size=4)+
    annotate("text", x=4, y=0.4, label="D", size=4)+
    annotate("text", x=4, y=0.64, label="C", size=4))
  


## Symbiont Counts #### 
#these data are already in the master data sheet
ols_sym<-lm(dil_corr_symdens_percm2~treatment, data=all_data)
summary(ols_sym)
# Call:
#   lm(formula = dil_corr_symdens_percm2 ~ treatment, data = all_data)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2578695 -1200246  -393141   881782  5883657 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)          2781893     566071   4.914 1.17e-05 ***
#   treatmentfed_hot    -1203103     769139  -1.564   0.1246    
# treatmentunfed_ctrl    94204     769139   0.122   0.9031    
# treatmentunfed_hot  -1868634     769139  -2.430   0.0191 *  
# 

# Post-hoc pairwise comparisons
summary(glht(ols_sym, linfct = mcp(treatment = "Tukey")))
# Multiple Comparisons of Means: Tukey Contrasts
# Fit: lm(formula = dil_corr_symdens_percm2 ~ treatment, data = all_data)
# Linear Hypotheses:
# Estimate Std. Error t value Pr(>|t|)  
# fed_hot - fed_ctrl == 0     -1203103     769139  -1.564   0.4085  
# unfed_ctrl - fed_ctrl == 0     94204     769139   0.122   0.9993  
# unfed_hot - fed_ctrl == 0   -1868634     769139  -2.430   0.0857 .
# unfed_ctrl - fed_hot == 0    1297306     736394   1.762   0.3045  
# unfed_hot - fed_hot == 0     -665531     736394  -0.904   0.8028  
# unfed_hot - unfed_ctrl == 0 -1962838     736394  -2.665   0.0499 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)

# only unfed, hot vs unfed, ctrl is sig, will reflect that in plot 

(b<-ggplot(all_data, aes(x=treatment, y=dil_corr_symdens_percm2))+
    geom_jitter(colour="grey30", width = 0.2, height = 0)+
    geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA) +
    xlab("") + ylab(expression(Symbiont~Density~(million~cells/cm^2)))+theme_bw()+
    theme(axis.text.x = element_blank())+
    scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                     labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+ 
    scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                      values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                      labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
    guides(colour="none", fill="none")+
    scale_y_continuous(breaks=c(0,2000000,4000000,6000000), 
                       limits=c(0,7000000),
                       labels=c("0", "2","4","6"))+
    annotate("text", x=1.2, y=6900000, label="italic(Sym~Density)", size=5, parse=T)+
    annotate("text", x=1, y=6000000, label="A,B", size=4)+
    annotate("text", x=2, y=6000000, label="A", size=4)+
    annotate("text", x=3, y=6000000, label="A,B", size=4)+
    annotate("text", x=4, y=6000000, label="B", size=4))


## Chlorophyll ####
chl<-read.table("input_files/physio/chl_data.txt", header=T)

mutate(chl, avg=rowMeans(dplyr::select(chl, rep1,rep2, rep3), na.rm = T))->chl
pivot_wider(chl, names_from = wavelength, values_from = avg, id_cols = coral_id)->chl

#These equations come from Jeffrey and Haxo 1968 for 90% acetone extraction of dinoflagellates
chl$chlA <- 13.31*chl$a663 - 0.27*chl$a630
chl$chlC2 <- -8.37*chl$a663 + 51.72*chl$a630

# a few samples (C1, Q5, and O11) didn't seem to work properly and have negative values, 
# there's also one strangely high sample change those to NA
bad_samps<-c("C1", "Q5", "O11", "J13")
chl[chl$coral_id %in% bad_samps,4:5]<-NA
#combine chl data with all_data df 

all_data<-left_join(all_data, dplyr::select(chl, coral_id, chlA, chlC2), by="coral_id")

# results are in ug/ml convert to need to get to ug/cm2 
# The factor of 20 is becuase the original equations assume 1 ml of extract and 200 ul was used here
all_data$chlA_total<-(all_data$chlA*20)/all_data$surface_area_cm2
all_data$chlC2_total<-(all_data$chlC2*20)/all_data$surface_area_cm2


ols_chla<-lm(chlA_total~treatment, data=all_data)
summary(ols_chla)
# Call:
#   lm(formula = chlA_total ~ treatment, data = all_data)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -14.6282  -3.9054  -0.9282   2.1793  20.3911 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           17.298      2.295   7.537 3.97e-09 ***
# treatmentfed_hot     -10.149      3.246  -3.127  0.00333 ** 
# treatmentunfed_ctrl   -2.650      3.118  -0.850  0.40058    
# treatmentunfed_hot   -12.005      3.537  -3.394  0.00159 ** 

summary(glht(ols_chla, linfct = mcp(treatment = "Tukey")))
# Linear Hypotheses:
# Estimate Std. Error t value Pr(>|t|)   
# fed_hot - fed_ctrl == 0      -10.149      3.246  -3.127  0.01669 * 
# unfed_ctrl - fed_ctrl == 0    -2.650      3.118  -0.850  0.82974   
# unfed_hot - fed_ctrl == 0    -12.005      3.537  -3.394  0.00812 **
# unfed_ctrl - fed_hot == 0      7.499      3.118   2.405  0.09266 . 
# unfed_hot - fed_hot == 0      -1.857      3.537  -0.525  0.95240   
# unfed_hot - unfed_ctrl == 0   -9.355      3.420  -2.735  0.04414 * 


marginal  = emmeans(ols_chla, ~ treatment)
cld(marginal, Letters=letters)
# treatment  emmean   SE df lower.CL upper.CL .group
# unfed_hot    5.29 2.69 39   -0.151     10.7  a    
# fed_hot      7.15 2.30 39    2.507     11.8  ab   
# unfed_ctrl  14.65 2.11 39   10.378     18.9   bc  
# fed_ctrl    17.30 2.30 39   12.656     21.9    c  


(c<-ggplot(all_data, aes(x=treatment, y=chlA_total))+
    geom_jitter(colour="grey30",width = 0.2, height = 0)+
    geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA)+
    xlab("") + ylab(expression(Chl[A]~Concentration~(ug/cm^2)))+ theme_bw()+
    theme(axis.text.x = element_blank())+
    scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                     labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
    scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                      values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                      labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
    guides(colour="none", fill="none")+
    scale_y_continuous(breaks=c(0,15,30,45))+
    annotate("text", x=1, y=50, label="italic(Chl[A])", size=5, parse=T)+
    annotate("text", x=1, y=44, label="A", size=4)+
    annotate("text", x=2, y=44, label="A,B", size=4)+
    annotate("text", x=3, y=44, label="B,C", size=4)+
    annotate("text", x=4, y=44, label="C", size=4))


ols_chlc<-lm(chlC2_total~treatment, data=all_data)
summary(ols_chlc)
summary(glht(ols_chlc, linfct = mcp(treatment = "Tukey")))
# Linear Hypotheses:
# Estimate Std. Error t value Pr(>|t|)  
# fed_hot - fed_ctrl == 0      -6.1756     3.5233  -1.753   0.3103  
# unfed_ctrl - fed_ctrl == 0    0.4642     3.3851   0.137   0.9991  
# unfed_hot - fed_ctrl == 0    -9.1501     3.8395  -2.383   0.0967 .
# unfed_ctrl - fed_hot == 0     6.6399     3.3851   1.961   0.2194  
# unfed_hot - fed_hot == 0     -2.9745     3.8395  -0.775   0.8650  
# unfed_hot - unfed_ctrl == 0  -9.6144     3.7130  -2.589   0.0615 .


marginal  = emmeans(ols_chlc, ~ treatment)
cld(marginal, Letters=letters)
# treatment  emmean   SE df lower.CL upper.CL .group
# unfed_hot    4.17 2.92 39    -1.74     10.1  a    
# fed_hot      7.15 2.49 39     2.11     12.2  a    
# fed_ctrl    13.32 2.49 39     8.28     18.4  a    
# unfed_ctrl  13.79 2.29 39     9.15     18.4  a    


(d<-ggplot(all_data, aes(x=treatment, y=chlC2_total))+
    geom_jitter(colour="grey30",width = 0.2, height = 0)+
    geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA) +
    xlab("") +theme_bw()+
    theme(axis.text.x = element_blank())+
    ylab(expression(Chl[C2]~Concentration~(ug/cm^2)))+
    scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                   labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
    scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                    values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                    labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
    guides(colour="none", fill="none")+
    scale_y_continuous(limits=c(0,34), breaks=c(0,10,20,30))+
    annotate("text", x=1, y=33, label="italic(Chl[C2])", size=5, parse=T))


## Protein ####

#See the website below for more information about how Thermo Fisher 
#suggests analyzing Bradford Assay data: 
#http://tools.thermofisher.com/content/sfs/brochures/TR0057-Read-std-curves.pdf

# Load the standards reading for each plate
d <- data.frame(read.csv('input_files/physio/standard_curves.csv', header = T))
d$plate<-as.factor(d$plate)

coral_prot<-read.table("input_files/physio/coral_protein_data.txt",header=T)

#calculates the average absorbency for each row 
mutate(d, avg_abs=(abs1+abs2+abs3)/3)->d
# this adds a new column to the data frame that has the average fluorescence value for all three replicates

#plots of standard curves
for (plateID in levels(d$plate)){
  p<-ggplot(data=subset(d, plate==paste(plateID)), aes(x=avg_abs, y=conc))+
    theme_bw()+
    ggtitle(paste(plateID))+
    geom_point()+
    geom_smooth(method = "lm", formula = y ~ stats::poly(x, 3, raw=TRUE))+
    ylab("Known standard concentration")+
    xlab("Average absorbance (n=3 replicates)")
  print(p)
} 

# Generate polynomial curves and save model outputs 
# Fit a two order polynomial to the data to generate equation for calculating protein concentration from absorbency values 
# This code loops through the plate IDs runs the model and then extracts the intercept and coefficients into a data frame. 
# You can then load a dataframe that has your samples and their corresponding plate IDs and merge the two (done below) 
# so that you can calculate protein values 

# Initializes and empty data frame to hold values
# modify this to however many plate you have!
plates=6
model_results<-data.frame(plateID=levels(d$plate), intercept=c(rep(0,plates)), coef1=c(rep(0,plates)), coef2=c(rep(0,plates)), coef3=c(rep(0,plates)))

# The individual models are saved as PlateID_model 
# The coefficients for each are then grabbed and put into model_results 

# Uses a polynomial model (order 3), as this was found to be the best fit during initial excel exploration for all the plates
for (plateID in levels(d$plate)){
  assign(paste(plateID, "_model", sep=""), lm(conc ~ stats::poly(avg_abs,3, raw=TRUE), data=subset(d, plate==plateID)))
  model_results[model_results$plateID==plateID,]$intercept<-get(paste(plateID, "_model", sep=""))$coefficients[1]
  model_results[model_results$plateID==plateID,]$coef1<-get(paste(plateID, "_model", sep=""))$coefficients[2]
  model_results[model_results$plateID==plateID,]$coef2<-get(paste(plateID, "_model", sep=""))$coefficients[3]
  model_results[model_results$plateID==plateID,]$coef3<-get(paste(plateID, "_model", sep=""))$coefficients[4]
}

# Convert sample to protein values 

# Now use the equations to convert you sample values to protein concentrations 
# 
# Step 1: Convert raw values from assay to normalized values based on standard curve (this will be plate specific)
# From Step 1 we will get a concentration of carb for each sample in ug/uL
# 
# Step 2: Then we want to convert that to TOTAL carb in that sample, based on the slurry volume. This will require multiplying the slurry volume by 1000 to convert from ml to ul and then also dividing the carb concentration by 100 to convert from ug to mg
# 
# Step 3: Convert the total carb into a carb concentration/cm^2 using the surface area of each sample (this final value should be in mg/cm2)

# Read in sample data
# You will need to have a column that identifies what plate that sample was read on. 
# Since both host and symbiont protein was assessed for each sample, we need a column for Host_Plate_ID and one for Sym_Plate_ID

# Then you can just join your model results with your sample data frame so that you 
# can easily convert the sample absorbency to protein concentrations based on the curve/model for each plate. 
# 
# So now we need to match the host/sym plate ID for each sample in order to calculate the corrected absorbency values 
# We will combine the model results with the prot_metadata data frame

model_results$prot_host_plateID<-model_results$plateID
coral_prot<-left_join(coral_prot, dplyr::select(model_results, -plateID), by="prot_host_plateID")

model_results$prot_sym_plateID<-model_results$plateID
coral_prot<-left_join(coral_prot, dplyr::select(model_results, -prot_host_plateID, -plateID), by="prot_sym_plateID")

# this will make it so that the columns labeled 
# intercept.x, coef1.x, coef2.x, coef3.x correspond to the host curve parameters 
# and those that are .y correspond to the sym curve parameters these .x and .y are 
# unrelated to the curve equation below. 

# Then we can use the equtions to solve for absorbency 
# y = ax^3+ bx^2 + cx + d
# 
# y = protein concentration
# x = raw absorbency from Bradford Assay

# convert host sample absorbency to protein concentration
coral_prot$host_prot1_cor<- (coral_prot$host_prot1)^3 *(coral_prot$coef3.x) + (coral_prot$host_prot1)^2 *(coral_prot$coef2.x) + (coral_prot$host_prot1)*(coral_prot$coef1.x) + (coral_prot$intercept.x)
coral_prot$host_prot2_cor<- (coral_prot$host_prot2)^3 *(coral_prot$coef3.x) + (coral_prot$host_prot2)^2 *(coral_prot$coef2.x) + (coral_prot$host_prot2)*(coral_prot$coef1.x) + (coral_prot$intercept.x)
coral_prot$host_prot3_cor<- (coral_prot$host_prot3)^3 *(coral_prot$coef3.x) + (coral_prot$host_prot3)^2 *(coral_prot$coef2.x) + (coral_prot$host_prot3)*(coral_prot$coef1.x) + (coral_prot$intercept.x)
# these columns should be in ug/ul 

# convert sym sample absorbency to protein concentration}
coral_prot$sym_prot1_cor <- (coral_prot$sym_prot1)^3 *(coral_prot$coef3.y) + (coral_prot$sym_prot1)^2 *(coral_prot$coef2.y) + (coral_prot$sym_prot1)*(coral_prot$coef1.y) + (coral_prot$intercept.y)
coral_prot$sym_prot2_cor <- (coral_prot$sym_prot2)^3 *(coral_prot$coef3.y) + (coral_prot$sym_prot2)^2 *(coral_prot$coef2.y) + (coral_prot$sym_prot2)*(coral_prot$coef1.y) + (coral_prot$intercept.y)
coral_prot$sym_prot3_cor <- (coral_prot$sym_prot3)^3 *(coral_prot$coef3.y) + (coral_prot$sym_prot3)^2 *(coral_prot$coef2.y) + (coral_prot$sym_prot3)*(coral_prot$coef1.y) + (coral_prot$intercept.y)
# these columns should be in ug/ul 

#remove any nonsense negative values 
replace(coral_prot[,18:23], coral_prot[,18:23]<0, NA)->fixed
coral_prot<-cbind(coral_prot[,1:17], fixed)

#average corrected values
mutate(coral_prot, 
       host_cor_prot_avg=apply(coral_prot[,18:20], MARGIN = 1, function(x) mean(x,na.rm=TRUE)), 
       sym_cor_prot_avg=apply(coral_prot[,21:23], MARGIN = 1, function(x) mean(x,na.rm=TRUE)))->coral_prot

#add in the average, standard corrected values to the all_data dataframe 
all_data<-left_join(all_data, dplyr::select(coral_prot,coral_id, host_cor_prot_avg, sym_cor_prot_avg), by="coral_id")

# Convert values to final units of mg/cm2
# multiply slurry volume (in ul) with sym_cor_avg and host_cor_avg to get total protein for sample, 
# then divide by 1000 to get total protein in mg instead of ug (these cancel) so then
# divide by surface area
# For symbiont date the slurry volume is only 5 ml as the symbionts are pelleted and then resuspended in 5 ml
all_data$host_prot_concentration <- (all_data$host_cor_prot_avg*all_data$slurry_vol_ml)/(all_data$surface_area_cm2)
all_data$sym_prot_concentration <- (all_data$sym_cor_prot_avg*5)/(all_data$surface_area_cm2)
# gives final protein concentration in mg/cm2


ols_hostp<-lm(host_prot_concentration~treatment, data=all_data)
summary(ols_hostp)
summary(glht(ols_hostp, linfct = mcp(treatment = "Tukey")))
# Linear Hypotheses:
# Estimate Std. Error t value Pr(>|t|)    
# fed_hot - fed_ctrl == 0     -0.006004   0.001655  -3.628  0.00352 ** 
# unfed_ctrl - fed_ctrl == 0  -0.004340   0.001655  -2.623  0.05503 .  
# unfed_hot - fed_ctrl == 0   -0.008377   0.001655  -5.062  < 0.001 ***
# unfed_ctrl - fed_hot == 0    0.001665   0.001584   1.051  0.72065    
# unfed_hot - fed_hot == 0    -0.002373   0.001584  -1.498  0.44711    
# unfed_hot - unfed_ctrl == 0 -0.004037   0.001584  -2.548  0.06570 .  

marginal  = emmeans(ols_hostp, ~ treatment)
cld(marginal, Letters=letters)
# treatment   emmean      SE df lower.CL upper.CL .group
# unfed_hot  0.00512 0.00112 46  0.00286  0.00737  a    
# fed_hot    0.00749 0.00112 46  0.00523  0.00974  a    
# unfed_ctrl 0.00915 0.00112 46  0.00690  0.01141  ab   
# fed_ctrl   0.01349 0.00122 46  0.01104  0.01595   b   


# Plot final sample data 
 (e<-ggplot(all_data, aes(x=treatment, y=host_prot_concentration))+
    geom_jitter(colour="grey30",width = 0.2, height = 0)+
    geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA) +
    xlab("") + 
    ylab(expression(Host~Protein~(mg/cm^2)))+
    theme_bw()+theme(axis.text.x = element_blank())+
    scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                     labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
    scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                      values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                      labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
    guides(colour="none", fill="none")+
    scale_y_continuous(breaks=c(0,0.006,0.012,0.018), limits=c(0,0.023))+
    annotate("text", x=1.2, y=0.023, label="italic(Host~Protein)", size=5, parse=T)+
     annotate("text", x=1, y=0.021, label="A", size=4)+
     annotate("text", x=2, y=0.021, label="A,B", size=4)+
     annotate("text", x=3, y=0.021, label="B", size=4)+
     annotate("text", x=4, y=0.021, label="B", size=4))
  
ols_symp<-lm(sym_prot_concentration~treatment, data=all_data)
summary(ols_symp)
summary(glht(ols_symp, linfct = mcp(treatment = "Tukey")))

# Linear Hypotheses:
#   Estimate Std. Error t value Pr(>|t|)
# fed_hot - fed_ctrl == 0     -5.209e-05  4.953e-04  -0.105    1.000
# unfed_ctrl - fed_ctrl == 0  -6.573e-04  4.953e-04  -1.327    0.551
# unfed_hot - fed_ctrl == 0   -7.732e-04  5.047e-04  -1.532    0.427
# unfed_ctrl - fed_hot == 0   -6.052e-04  4.742e-04  -1.276    0.582
# unfed_hot - fed_hot == 0    -7.211e-04  4.840e-04  -1.490    0.452
# unfed_hot - unfed_ctrl == 0 -1.159e-04  4.840e-04  -0.239    0.995
# (Adjusted p values reported -- single-step method)

marginal  = emmeans(ols_symp, ~ treatment)
cld(marginal, Letters=letters)
# treatment   emmean       SE df lower.CL upper.CL .group
# unfed_hot  0.00147 0.000349 45 0.000770  0.00218  a    
# unfed_ctrl 0.00159 0.000335 45 0.000914  0.00226  a    
# fed_hot    0.00219 0.000335 45 0.001519  0.00287  a    
# fed_ctrl   0.00225 0.000365 45 0.001512  0.00298  a   

(f<-ggplot(all_data, aes(x=treatment, y=sym_prot_concentration))+
    geom_jitter(colour="grey30",width = 0.2, height = 0)+
    geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA) +
    xlab("") + 
    ylab(expression(Sym~Protein~(mg/cm^2)))+
    theme_bw()+theme(axis.text.x = element_blank())+
    scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                     labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
    scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                      values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                      labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
    guides(colour="none", fill="none")+
    scale_y_continuous(breaks=c(0,0.002,0.004,0.006), limits=c(0,0.0063))+
    annotate("text", x=1.1, y=0.0062, label="italic(Sym~Protein)", size=5, parse=T))
  

## Carbohydrates ####
carbs_stdcurves <- data.frame(read.csv('input_files/physio/carb_std_curves.csv', header = T))
carbs_stdcurves$plate<-as.factor(carbs_stdcurves$plate)

coral_carb<-read.table("input_files/physio/coral_carb_data.txt", header=T)

# calculate the average absorbency
# apply calculates the mean across columns 3 to 8 over each row, ignoring the missing NA values when they occur
mutate(carbs_stdcurves, 
       avg_abs=apply(carbs_stdcurves[,3:8], MARGIN = 1, function(x) mean(x,na.rm=TRUE)))->carbs_stdcurves

#plot the average values with linear fits
for (plateID in levels(carbs_stdcurves$plate)){
  p<-ggplot(data=subset(carbs_stdcurves, plate==paste(plateID)), aes(x=avg_abs, y=conc))+
    theme_bw()+
    ggtitle(paste(plateID))+
    geom_point()+
    geom_smooth(method = "lm", formula = y ~ stats::poly(x, 1, raw=TRUE))+
    ylab("Known standard concentration")+
    xlab("Average absorbance (n=3-6 replicates")
  print(p)
} 

# Initializes and empty data frame to hold values
# modify this to however many plate you have!
plates=7
model_results_carbs<-data.frame(plateID=levels(carbs_stdcurves$plate), intercept=c(rep(0,plates)), coef1=c(rep(0,plates)))

# The individual models are saved as PlateID_model 
# The coefficients for each are then grabbed and put into model_results 

# Uses a polynomial model (order 3), as this was found to be the best fit during initial excel exploration for all the plates
for (plateID in levels(carbs_stdcurves$plate)){
  assign(paste(plateID, "_model", sep=""), lm(conc ~ stats::poly(avg_abs,1, raw=TRUE), data=subset(carbs_stdcurves, plate==plateID)))
  model_results_carbs[model_results_carbs$plateID==plateID,]$intercept<-get(paste(plateID, "_model", sep=""))$coefficients[1]
  model_results_carbs[model_results_carbs$plateID==plateID,]$coef1<-get(paste(plateID, "_model", sep=""))$coefficients[2]
}

model_results_carbs$carb_host_plateID<-model_results_carbs$plateID
coral_carb<-left_join(coral_carb, dplyr::select(model_results_carbs, -plateID), by="carb_host_plateID")

model_results_carbs$carb_sym_plateID<-model_results_carbs$plateID
coral_carb<-left_join(coral_carb, dplyr::select(model_results_carbs, -carb_host_plateID, -plateID), by="carb_sym_plateID")

# this will make it so that the columns labeled 
# intercept.x, coef1.x, correspond to the host curve parameters 
# and those that are .y correspond to the sym curve parameters these .x and .y are 
# unrelated to the curve equation below. 

# Then we can use the equations to solve for absorbency 
# y = ax + b
# y = carbohydrate concentration
# x = raw absorbency from assay

# Now use the equations to convert sample values to carb concentrations 
# 
# Step 1: Convert raw values from assay to normalized values based on standard curve (this will be plate specific)
# From Step 1 we will get a concentration of protein for each sample in mg/ml
# 
# Step 2: Then we want to convert that to TOTAL protein in that sample, based on the slurry volume. 
# This will require multiplying the slurry volume 
# 
# Step 3: Convert the total protein into a protein concentration/cm^2 using the surface area of each sample (this final value should be in mg/cm2)
# 

# convert host sample absorbency to carb concentration
coral_carb$host_carb1_cor<- (coral_carb$rawcarbhost_1)*(coral_carb$coef1.x) + (coral_carb$intercept.x)
coral_carb$host_carb2_cor<- (coral_carb$rawcarbhost_2)*(coral_carb$coef1.x) + (coral_carb$intercept.x)
coral_carb$host_carb3_cor<- (coral_carb$rawcarbhost_3)*(coral_carb$coef1.x) + (coral_carb$intercept.x)
# these columns should be in mg/ml 

# convert sym sample absorbency to carb concentration}
coral_carb$sym_carb1_cor <- (coral_carb$rawcarbsym_1)*(coral_carb$coef1.y) + (coral_carb$intercept.y)
coral_carb$sym_carb2_cor <- (coral_carb$rawcarbsym_2)*(coral_carb$coef1.y) + (coral_carb$intercept.y)
coral_carb$sym_carb3_cor <- (coral_carb$rawcarbsym_3)*(coral_carb$coef1.y) + (coral_carb$intercept.y)
# these columns should be in mg/ml 

#average corrected values
mutate(coral_carb, 
       host_cor_carb_avg=apply(coral_carb[,14:16], MARGIN = 1, function(x) mean(x,na.rm=TRUE)), 
       sym_cor_carb_avg=apply(coral_carb[,17:19], MARGIN = 1, function(x) mean(x,na.rm=TRUE)))->coral_carb

#add in the average, standard corrected values to the all_data dataframe 
all_data<-left_join(all_data, dplyr::select(coral_carb, coral_id,host_cor_carb_avg,sym_cor_carb_avg), by="coral_id")

# Convert values to final units of mg/cm2
# multiply slurry volume  with sym_cor_avg and host_cor_avg to get total carb for sample, 
all_data$host_carb_concentration <- (all_data$host_cor_carb_avg*all_data$slurry_vol_ml)/(all_data$surface_area_cm2)
all_data$sym_carb_concentration <- (all_data$sym_cor_carb_avg*5)/(all_data$surface_area_cm2)
#there's a really high sym carb outlier remove: 
all_data$sym_carb_concentration[37]<-NA


ols_hostc<-lm(host_carb_concentration~treatment, data=all_data)
summary(ols_hostc)
summary(glht(ols_hostc, linfct = mcp(treatment = "Tukey")))
# Linear Hypotheses:
# Estimate Std. Error t value Pr(>|t|)
# fed_hot - fed_ctrl == 0     -0.0078735  0.0046888  -1.679    0.346
# unfed_ctrl - fed_ctrl == 0   0.0009532  0.0046888   0.203    0.997
# unfed_hot - fed_ctrl == 0   -0.0094871  0.0046888  -2.023    0.194
# unfed_ctrl - fed_hot == 0    0.0088267  0.0044892   1.966    0.216
# unfed_hot - fed_hot == 0    -0.0016136  0.0044892  -0.359    0.984
# unfed_hot - unfed_ctrl == 0 -0.0104403  0.0044892  -2.326    0.107
# (Adjusted p values reported -- single-step method)

marginal  = emmeans(ols_hostc, ~ treatment)
cld(marginal, Letters=letters)
# treatment  emmean      SE df lower.CL upper.CL .group
# unfed_hot  0.0194 0.00317 46   0.0130   0.0258  a    
# fed_hot    0.0210 0.00317 46   0.0146   0.0274  a    
# fed_ctrl   0.0289 0.00345 46   0.0219   0.0358  a    
# unfed_ctrl 0.0298 0.00317 46   0.0234   0.0362  a    


# Plot final sample data 
(g<-ggplot(all_data, aes(x=treatment, y=host_carb_concentration))+
  geom_jitter(colour="grey30",width = 0.2, height = 0)+
  geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA) +
  xlab("") + 
  ylab(expression(Host~Carbohydrates~(mg/cm^2)))+
  theme_bw()+
  scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                   labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
  scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                    values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                    labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
  guides(colour="none", fill="none")+
  scale_y_continuous(breaks=c(0.015,0.03,0.045), limits=c(0.01,0.049))+
  annotate("text", x=1.1, y=0.049, label="italic(Host~Carbs)", size=5, parse=T))

ols_symc<-lm(sym_carb_concentration~treatment, data=all_data)
summary(ols_symc)
summary(glht(ols_symc, linfct = mcp(treatment = "Tukey")))
# Linear Hypotheses:
# Estimate Std. Error t value Pr(>|t|)    
# fed_hot - fed_ctrl == 0     -0.0074686  0.0013475  -5.543  < 0.001 ***
# unfed_ctrl - fed_ctrl == 0  -0.0045498  0.0013475  -3.376  0.00801 ** 
# unfed_hot - fed_ctrl == 0   -0.0077333  0.0013730  -5.632  < 0.001 ***
# unfed_ctrl - fed_hot == 0    0.0029188  0.0012901   2.262  0.12239    
# unfed_hot - fed_hot == 0    -0.0002647  0.0013167  -0.201  0.99708    
# unfed_hot - unfed_ctrl == 0 -0.0031835  0.0013167  -2.418  0.08828 .  


marginal  = emmeans(ols_symc, ~ treatment)
cld(marginal, Letters=letters)
# treatment   emmean       SE df lower.CL upper.CL .group
# unfed_hot  0.00500 0.000950 45  0.00309  0.00691  a    
# fed_hot    0.00527 0.000912 45  0.00343  0.00710  a    
# unfed_ctrl 0.00818 0.000912 45  0.00635  0.01002  a    
# fed_ctrl   0.01273 0.000992 45  0.01074  0.01473   b   


(h<-ggplot(all_data, aes(x=treatment, y=sym_carb_concentration))+
  geom_jitter(colour="grey30",width = 0.2, height = 0)+
  geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA) +
  xlab("") + 
  ylab(expression(Sym~Carbohydrates~(mg/cm^2)))+
  theme_bw()+
  scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                   labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
  scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                    values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                    labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
  guides(colour="none", fill="none")+
  scale_y_continuous(breaks=c(0,0.005,0.01,0.015,0.02), limits=c(0,0.024))+
  annotate("text", x=1.1, y=0.024, label="italic(Sym~Carbs)", size=5, parse=T)+
    annotate("text", x=1, y=0.0215, label="A", size=4)+
    annotate("text", x=2, y=0.0215, label="B", size=4)+
    annotate("text", x=3, y=0.0215, label="B", size=4)+
    annotate("text", x=4, y=0.0215, label="B", size=4))

ggarrange(a,b,c,d,e,f,g,h, nrow=4, ncol=2, common.legend = T, legend = "top", labels = "AUTO")
ggsave("fig2.png", height=15, width=10, units="in", dpi=300)


## Physiology PCA ####
phys_data<-dplyr::select(all_data, coral_id, 
                         treatment,pam4, 
                         dil_corr_symdens_percm2,
                         sym_carb_concentration,
                         sym_prot_concentration,
                         host_carb_concentration, 
                         host_prot_concentration, 
                         chlA_total, chlC2_total)
#make some colnames nicer for plotting
colnames(phys_data)[4:10]<-c("sym_density", "sym_carbs", "sym_protein",
                             "host_carbs", "host_protein", "chlA", "chlC2")


#set missing values to zero for pca 
phys_data[is.na(phys_data)]<-0
phys_pca<-prcomp(phys_data[3:10], center = T, scale. = T)

adonis(phys_pca$x~treatment,data=phys_data, method = 'eu')
# Terms added sequentially (first to last)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  3    121.79  40.596  6.9109 0.31068  0.001 ***
# Residuals 46    270.21   5.874         0.68932           
# Total     49    392.00                 1.00000              
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

ggbiplot(phys_pca, scale = 1, groups = phys_data$treatment, ellipse = T)+
  scale_colour_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                    values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                    labels=c("Fed, Ambient","Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
  ylab("PC2 (16.8%)")+
  xlab("PCA1 (33.7%)")+
  theme_bw()+
  theme(legend.position = "bottom", legend.background = element_rect(fill="white", colour="grey60"))+
  labs(colour="Treatment")+
  annotate("text", x=1.5, y=-2.5, label=paste("italic(Adonis) * ' p'[treatment]<0.001"), parse=TRUE)

ggsave("figure3.png", height = 7, width=8, units="in", dpi=300)  

## Surface area supp ####
#these data are already in the master data sheet
ols_sa<-lm(surface_area_cm2~treatment, data=all_data)
summary(ols_sa)
summary(glht(ols_sa, linfct = mcp(treatment = "Tukey")))
#no differences in SA between treatments

ggplot(all_data, aes(x=treatment, y=surface_area_cm2))+
   geom_jitter(colour="grey30",width = 0.2, height = 0)+
   geom_boxplot(alpha=0.7, aes(fill=treatment), colour="grey10", outlier.colour = NA)+
   xlab("") + 
   ylab(expression(Live~tissue~surface~area~(cm^2)))+
   theme_bw()+
   scale_x_discrete(limits=c("fed_ctrl","unfed_ctrl","fed_hot", "unfed_hot"),
                    labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
   scale_fill_manual(limits=c("fed_ctrl","unfed_ctrl", "fed_hot", "unfed_hot"),
                     values=c("navyblue", "cyan3",  "darkgoldenrod2","firebrick"),
                     labels=c("Fed, Ambient","Unfed, Ambient", "Unfed, Heated", "Fed, Heated"))+
   guides(colour="none", fill="none")

ggsave("figureS2.png", height=4, width = 5, units = "in", dpi=300)

## Gene Expression ####
# See read_mapping_scripts for a description and scripts/commands used for the 
# processing raw reads that generated the input file used here. 

cts_raw<-read.table("input_files/gene_expression/ocu_het_counts.txt", header=TRUE, row.names = 1) 
# need to remove data for symbiont genes, find row when sym counts start
head(grep("Sym.*", row.names(cts_raw),ignore.case = T))
#[1] 23017 23018 23019 23020 23021 23022
# that's the row before symbiont contigs start
cts_raw<-cts_raw[1:23016,]

cts_sym<-cts_raw[23017:45019,]

coldata<-read.table("input_files/gene_expression/ocu_het_samp_data.txt", header=TRUE)
rownames(coldata)<-coldata$coral_id
#order coldata by coral_id to match with cts
coldata<-coldata[order(row.names(coldata)),]
#filter coldata to only have the samples that are in cts
coldata<-coldata[coldata$coral_id %in% colnames(cts_raw),]

#look at total counts by sample
totalCounts=colSums(cts_raw)
barplot(totalCounts, col="coral")
totalCounts
# B1     B2     B5     B8     C2     C4     E3     G1     G3     I5    J13    J15     J9     L4     M1     M5    N14     N4    O11 
# 180855 148620   1537 297235  47455   3872  11104   2121  30335  84049 811642 347231 396993 110910   1315 209896 445580 165288   3617 
# Q1     Q6     R2     R7 
# 197692 220261  20984 152505 
# Want to keep the samples with more than 30K reads to avoid biases DEG expression by low count samples
cts_raw<-cts_raw[,totalCounts>30000]
coldata<-coldata[totalCounts>30000,]

#make Table 1
table1<-data.frame(coral_id=colnames(cts_raw))
table1<-left_join(table1, dplyr::select(all_data, genotype, treatment, coral_id), by="coral_id")
table1$total_counts<-totalCounts

#write.table(table1, file="/Users/hannyrivera/Documents/BU/DaviesLab/_Data/Oculina_Exp/table1.txt")


## Check for outliers using PCA on rlog counts 
dds_nm<-DESeq2::DESeqDataSetFromMatrix(countData = cts_raw, colData = coldata, design=~1)  
rlog_nm=rlogTransformation(dds_nm, blind=TRUE) 
pca = prcomp(t(assay(rlog_nm)), center = TRUE, scale. = FALSE)

cts_pca_ind <- get_pca_ind(pca)
# 2. Coordinate of groups
cts_pca_ind_coord <- cts_pca_ind$coord %>% as_tibble() %>% dplyr::select(Dim.1, Dim.2)
cts_pca_ind_coord<-cbind(cts_pca_ind_coord, dplyr::select(coldata, coral_id, treatment))

ggplot(cts_pca_ind_coord, aes(Dim.1, Dim.2, colour=treatment, label=coral_id)) +
  geom_hline(yintercept =0)+geom_vline(xintercept =0)+
  geom_point()+
  geom_text(hjust = 0, nudge_y = 0.1, nudge_x = 0.2, check_overlap = T)+
  theme_bw()+
  scale_color_manual(limits=c("unfed_ctrl", "fed_ctrl", "unfed_hot", "fed_hot"),
                     values=c("navyblue", "cyan3", "firebrick", "darkgoldenrod2"),
                     labels=c("Unfed, Control", "Fed, Control", "Unfed, Heated", "Fed, Heated"))

# We don't seem to have any specific outlier samples after filtering out those with low reads



### Run DeSeq2 ####
dds<-DESeq2::DESeqDataSetFromMatrix(countData = cts_raw, colData = coldata, design=~treatment) # this is where you would need to edit the formula using design, 
# pre-filtering to remove low count genes 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
## The way the contrasts are written here will means 
## that a positive log fold change means it was up-regulated in the non control treatment 


### Unfed control temp vs fed control temp ####
ddsMF<-DESeq2::DESeq(dds)
resMF_ufct<-results(ddsMF, contrast=c("treatment", "unfed_ctrl",  "fed_ctrl"))
head(resMF_ufct)
summary(resMF_ufct)
# out of 17300 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 0, 0%
# outliers [1]       : 21, 0.12%
# low counts [2]     : 0, 0%
# (mean count < 0)

# generate the heats data for GOWMU
heats_ufct<-mutate(as.data.frame(resMF_ufct),negP=-log(pvalue))
heats_ufct<-mutate(heats_ufct, signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
heats_ufct$genes<-rownames(resMF_ufct)
heats_ufct_signedP<-dplyr::select(heats_ufct, genes, signedlogP)
write.csv(heats_ufct_signedP, file="input_files/gene_expression/heats_ufct_signedP.csv", row.names = F, quote = F)
heats_ufct_lfc<-dplyr::select(heats_ufct, genes,log2FoldChange)
write.csv(heats_ufct_lfc, file="input_files/gene_expression/heats_ufct_lfc.csv", row.names = F, quote = F)

valState_ufct=cbind(resMF_ufct$log2FoldChange, resMF_ufct$pvalue, resMF_ufct$padj)
colnames(valState_ufct)=c("log2change","pval.fc_fh", "padj.fc_fh")
rownames(valState_ufct)<-rownames(resMF_ufct)
head(resMF_ufct)

### Fed hot temp vs fed control temp ####
resMF_fdhot<-results(ddsMF, contrast=c("treatment", "fed_hot",  "fed_ctrl"))
head(resMF_fdhot)
summary(resMF_fdhot)
# out of 17300 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 319, 1.8%
# LFC < 0 (down)     : 356, 2.1%
# outliers [1]       : 21, 0.12%
# low counts [2]     : 10719, 62%
# (mean count < 4)

write.table(resMF_fdhot, file="input_files/gene_expression/deg_fc_fh.txt", quote=F, sep="\t")
valState_fdhot=cbind(resMF_fdhot$log2FoldChange, resMF_fdhot$pvalue, resMF_fdhot$padj)
colnames(valState_fdhot)=c("log2change","pval.fc_fh", "padj.fc_fh")
rownames(valState_fdhot)<-rownames(resMF_fdhot)
head(valState_fdhot)

# generate the heats data for GOWMU
heats_fdhot<-mutate(as.data.frame(resMF_fdhot),negP=-log(pvalue))
heats_fdhot<-mutate(heats_fdhot,signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
heats_fdhot$genes<-rownames(resMF_fdhot)
heats_fdhot_signedP<-dplyr::select(heats_fdhot, genes, signedlogP)
write.csv(heats_fdhot_signedP, file="input_files/gene_expression//heats_fdhot_signedP.csv", row.names = F, quote = F)
heats_fdhot_lfc<-dplyr::select(heats_fdhot, genes,log2FoldChange)
write.csv(heats_fdhot_lfc, file="input_files/gene_expression/heats_fdhot_lfc.csv", row.names = F, quote = F)


### Unfed hot temp vs fed control temp ####
resMF_ufhot<-results(ddsMF, contrast=c("treatment", "unfed_hot", "fed_ctrl"))
head(resMF_ufhot)
summary(resMF_ufhot)
# out of 17300 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 31, 0.18%
# LFC < 0 (down)     : 63, 0.36%
# outliers [1]       : 21, 0.12%
# low counts [2]     : 13400, 77%
# (mean count < 9)

# generate heats data for GO MWU
heats_ufhot<-mutate(as.data.frame(resMF_ufhot),negP=-log(pvalue))
heats_ufhot<-mutate(heats_ufhot,signedlogP=case_when(log2FoldChange<0 ~ negP*-1, log2FoldChange>0 ~ negP))
heats_ufhot$genes<-rownames(resMF_ufhot)
heats_ufhot_signedP<-dplyr::select(heats_ufhot, genes, signedlogP)
write.csv(heats_ufhot_signedP, file="input_files/gene_expression/heats_ufhot_signedP.csv", row.names = F, quote = F)
heats_ufhot_lfc<-dplyr::select(heats_ufhot, genes,log2FoldChange)
write.csv(heats_ufhot_lfc, file="input_files/gene_expression//heats_ufhot_lfc.csv", row.names = F, quote = F)


write.table(resMF_ufhot, file="input_files/gene_expression/deg_fc_uh.txt", quote=F, sep="\t")
valState_ufhot=cbind(resMF_ufhot$log2FoldChange, resMF_ufhot$pvalue, resMF_ufhot$padj)
head(valState_ufhot)
colnames(valState_ufhot)=c("log2change", "pval.fc_uh", "padj.fc_uh")
rownames(valState_ufhot)<-rownames(resMF_ufhot)
head(valState_ufhot)

### PCA ####
# rlog transform the data for better plotting 
rlogMF=rlogTransformation(ddsMF, blind=TRUE) 
rld=assay(rlogMF)
colnames(rld)<-paste(coldata$coral_id, coldata$treatment, sep="_")

pcaData<-DESeq2::plotPCA(rlogMF, intgroup=c("feeding", "temp"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca = prcomp(t(assay(rlogMF)), center = TRUE, scale. = FALSE)

adonis(pca$x~feeding,data=coldata, method = 'eu')
adonis(pca$x~temp,data=coldata, method = 'eu')
adonis(pca$x~treatment,data=coldata, method = 'eu')

(a<-ggplot(pcaData, aes(PC1, PC2)) +
  geom_hline(yintercept =0)+geom_vline(xintercept =0) +
  geom_point(aes(colour=group), size=3)+
  stat_ellipse(geom="polygon", aes(fill=group), alpha=0.3, level=0.95)+
  xlab(paste0("PC1 (",percentVar[1],"%)")) +
  ylab(paste0("PC2 (",percentVar[2],"%)")) + 
  theme_bw()+
  scale_fill_manual(limits=c("fed:ctrl","unfed:ctrl","fed:hot","unfed:hot"),
                        values=c("navyblue", "cyan3", "darkgoldenrod2", "firebrick"),
                        labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
  scale_colour_manual(limits=c("fed:ctrl","unfed:ctrl","fed:hot","unfed:hot"),
                     values=c("navyblue", "cyan3", "darkgoldenrod2", "firebrick"),
                     labels=c("Fed, Ambient", "Unfed, Ambient", "Fed, Heated", "Unfed, Heated"))+
  annotate("text", x=50, y=40, label=paste("italic(Adonis) * ' p'[Temp]<0.01"), parse=TRUE)+
  annotate("text", x=50, y=30, label=paste("italic(Adonis) * ' p'[Feed]>0.2"), parse=TRUE)+
  annotate("text", x=50, y=20, label=paste("italic(Adonis) * ' p'[Treatment]<0.01"), parse=TRUE)+
  labs(colour="Treatment")+guides(fill="none")+
  theme(panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour="black", fill=NA, size=1.5), 
        axis.text = element_text(face="bold", colour="black", size=8),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA),
        legend.background = element_blank(),
        legend.text = element_text(face="bold", colour="black", size=8),
        axis.title = element_text(face="bold", colour = "black", size=10),
        legend.title = element_text(face="bold", colour="black", size=10))+
  scale_y_continuous(limits=c(-50,50))+
  scale_x_continuous(limits=c(-50,70)))

ggsave("figures/PCA_allgenes.png", dpi=300, height=4.5, width=6, units="in")

### Venn Diagrams of differentially expressed genes ###

fdhot_deg<-dplyr::filter(as.data.frame(valState_fdhot), padj.fc_fh<0.1)
fdhot_deg_up<-dplyr::filter(as.data.frame(valState_fdhot), padj.fc_fh<0.1 & log2change>0)
fdhot_deg_down<-dplyr::filter(as.data.frame(valState_fdhot), padj.fc_fh<0.1 & log2change<0)

fdhot_deg$gene<-row.names(fdhot_deg)
fdhot_deg_up$gene<-row.names(fdhot_deg_up)
fdhot_deg_down$gene<-row.names(fdhot_deg_down)

ufhot_deg<-dplyr::filter(as.data.frame(valState_ufhot), padj.fc_uh<0.1)
ufhot_deg_up<-dplyr::filter(as.data.frame(valState_ufhot), padj.fc_uh<0.1 & log2change>0)
ufhot_deg_down<-dplyr::filter(as.data.frame(valState_ufhot), padj.fc_uh<0.1 & log2change<0)

ufhot_deg$gene<-row.names(ufhot_deg)
ufhot_deg_up$gene<-row.names(ufhot_deg_up)
ufhot_deg_down$gene<-row.names(ufhot_deg_down)

ufct_deg<-dplyr::filter(as.data.frame(valState_ufct), padj.fc_fh<0.1)
ufct_deg$gene<-row.names(ufct_deg)


x <- list(
  A = fdhot_deg$gene,
  B = ufhot_deg$gene, 
  C = ufct_deg$gene
)

names(x) <- c("Fed, Heated","Unfed, Heated", "Unfed, Ambient")

(b<-ggvenn(x, fill_color = c("darkgoldenrod2", "firebrick", "cyan3"), text_size = 5, set_name_size = 5, show_percentage = F ))
ggarrange(a,b, ncol=1, nrow=2, common.legend = T, legend = "top", widths = c(1,1.5))
ggsave("Tag_seq_counts/fig4.png", width = 6, height = 8, units="in", dpi=300)

x2 <- list(
  A = fdhot_deg_down$gene,
  B = ufhot_deg_down$gene,
  C = fdhot_deg_up$gene,
  D = ufhot_deg_up$gene
)

ggvenn(x2)



### Heatmaps ####

## contrast for fed, heated
#remove samples not relevant for fed, ambient (fc) vs fed, heated (fh) contrast
rld_fc_fh<-rld[,grep("_fed", colnames(rld))]

rldpvals_fdhot=cbind(rld_fc_fh,valState_fdhot)
head(rldpvals_fdhot)
dim(rldpvals_fdhot)
# [1] 17300  12
table(complete.cases(rldpvals_fdhot))
# FALSE  TRUE 
# 10740  6560 
# the lowcount and outlier genes are listed as false since their p adjusted is NA
rld_data_fdhot<-as.data.frame(rldpvals_fdhot)


#Read in gene annotations
gg=read.delim("input_files/gene_expression/coral_gene_names.tab",sep="\t", header=FALSE)
colnames(gg)<-c("contig","genename")

adjpval=0.10 # FDR cutoff
conds_fdhot=dplyr::filter(rld_data_fdhot, padj.fc_fh < adjpval)
conds_fdhot$contig<-rownames(conds_fdhot)
conds_fdhot<-left_join(conds_fdhot, gg, by="contig") # add in genenames from annotations file
#write.table(degs, file="deg_list_exp.txt", quote=FALSE, row.names = FALSE, sep = "\t")
#rename NAs to unknowns
conds_fdhot$genename[is.na(conds_fdhot$genename)]<-"unknown"
#get rid of species designation on the genename 
conds_fdhot$genename<-sub("OS=.*", "", conds_fdhot$genename)

# only keep the top 30 most differential expressed DEGs for heatmap figure, 
# will do top 30 up and top 30 down reg genes
conds_top_fdhot<-rbind(conds_fdhot[order(conds_fdhot$log2change),][1:25,],conds_fdhot[order(conds_fdhot$log2change),][651:675,])

# scale_color_manual(limits=c("unfed_ctrl", "fed_ctrl", "unfed_hot", "fed_hot"),
#                    values=c("navyblue", "cyan3", "firebrick", "darkgoldenrod2"),
#                    labels=c("Unfed, Control", "Fed, Control", "Unfed, Heated", "Fed, Heated"))

ccol<-rev(colorRampPalette(brewer.pal(n=11, name="BrBG"))(20))

pheatmap(conds_top_fdhot[,1:9], show_rownames=T, labels_row=conds_top_fdhot$genename,
         show_colnames = T, angle_col = "45", labels_col = colnames(conds_top_fdhot[1:9]),
         cellheight = 8, cellwidth = 25, fontsize_row=5, cluster_cols=T, scale='row',
         cluster_rows = T, color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 10, treeheight_row = 10, legend=T,
         filename="clean_repo/figures/heatmap_top60_fdhot.png", width=10, height=10)


# contrast for unfed, heated
# remove samples not relevant for unfed_heat (uh) vs fed, ambient (fc) contrast
rld_fc_uh<-rld[,c(grep("_fed_ctrl", colnames(rld)),grep("_unfed_hot", colnames(rld)))]

rldpvals_ufhot=cbind(rld_fc_uh,valState_ufhot)
head(rldpvals_ufhot)
dim(rldpvals_ufhot)
# [1] 17300  10
table(complete.cases(rldpvals_ufhot))
# FALSE  TRUE 
# 13421  3879  
# the lowcount and outlier genes are listed as false since their p adjusted is NA
rld_data_ufhot<-as.data.frame(rldpvals_ufhot)


#Read in gene annotations
adjpval=0.10 # FDR cutoff
conds_ufhot=dplyr::filter(rld_data_ufhot, padj.fc_uh < adjpval)
conds_ufhot$contig<-rownames(conds_ufhot)
conds_ufhot<-left_join(conds_ufhot, gg, by="contig") # add in genenames from annotations file
#write.table(degs, file="deg_list_exp.txt", quote=FALSE, row.names = FALSE, sep = "\t")
#rename NAs to unknowns
conds_ufhot$genename[is.na(conds_ufhot$genename)]<-"unknown"
#get rid of species designation on the genename 
conds_ufhot$genename<-sub("OS=.*", "", conds_ufhot$genename)

ccol<-rev(colorRampPalette(brewer.pal(n=11, name="BrBG"))(20))

pheatmap(conds_ufhot[,1:7], show_rownames=T, labels_row=conds_ufhot$genename,
         show_colnames = T, angle_col = "45", labels_col = colnames(conds_ufhot[1:9]),
         cellheight = 8, cellwidth = 25, fontsize_row=5, cluster_cols=T, scale='row',
         cluster_rows = T, color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 10, treeheight_row = 10, legend=T,
         filename="clean_repo/figures/heatmap_alldegs_ufhot.png", width=10, height=20)


#find overlap genes between unfed hot and fed hot relative to fed ambient
conds_overlapping<-conds_ufhot[conds_ufhot$contig %in% conds_fdhot$contig,]
conds_overlapping<-left_join(conds_overlapping, conds_fdhot[,c(2,3,6,8,13)], by="contig")
conds_overlapping<-conds_overlapping[,c(1:7,13:16,8:12)]

pheatmap(conds_overlapping[,1:11], show_rownames=T, labels_row=conds_overlapping$genename,
         show_colnames = T, angle_col = "45", labels_col = colnames(conds_overlapping[1:11]),
         cellheight = 8, cellwidth = 25, fontsize_row=5, cluster_cols=T, scale='row',
         cluster_rows = T, color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 10, treeheight_row = 10, legend=T,
         filename="clean_repo/figures/heatmap_alldegs_overlapping.png", width=10, height=15)

# look at non overlaps
`%nin%` = Negate(`%in%`)

conds_nooverlap<-conds_ufhot[conds_ufhot$contig %nin% conds_fdhot$contig,]

pheatmap(conds_nooverlap[,1:7], show_rownames=T, labels_row=conds_nooverlap$genename,
         show_colnames = T, angle_col = "45", labels_col = colnames(conds_nooverlap[1:11]),
         cellheight = 8, cellwidth = 25, fontsize_row=5, cluster_cols=T, scale='row',
         cluster_rows = T, color=ccol, cutree_rows = 2, cutree_cols = 2, 
         treeheight_col = 10, treeheight_row = 10, legend=T,
         filename="clean_repo/figures/heatmap_ughtdegs_nooverlap.png", width=10, height=10)

conds_nooverlap_fdhot<-conds_fdhot[conds_fdhot$contig %nin% conds_ufhot$contig,]



## Delta Ranks #### 
# load results BP 
ufht_goBP=read.table("input_files/gene_expression/MWU_BP_heats_ufhot_signedP.csv",header=T, sep = " ")
fdht_goBP=read.table("input_files/gene_expression/MWU_BP_heats_fdhot_signedP.csv",header=T, sep = " ")
ufct_goBP=read.table("input_files/gene_expression/MWU_BP_heats_ufct_signedP.csv",header=T, sep = " ")

# keep only significant GO terms (p<0.1)
ufht_goBP<-dplyr::filter(ufht_goBP, p.adj<0.1)
fdht_goBP<-dplyr::filter(fdht_goBP, p.adj<0.1)
ufct_goBP<-dplyr::filter(ufct_goBP, p.adj<0.1)

# Compare between fed hot and unfed hot 
goods_ufht_fdhot=intersect(ufht_goBP$term,fdht_goBP$term)

ufht_goBP_comp1=ufht_goBP[ufht_goBP$term %in% goods_ufht_fdhot,]
#get GO terms from unfed list that are from the shared list 
fdht_goBP_comp1=fdht_goBP[fdht_goBP$term %in% goods_ufht_fdhot,]
# get GO terms from the fed list that are in shared list 

# all overlapping GO terms
sigs_ufht_fdht=merge(ufht_goBP_comp1,fdht_goBP_comp1,by="term")
dim(sigs_ufht_fdht) # 110 GO terms

# make simple GO category labels for coloring
sigs_ufht_fdht$short<-NA
sigs_ufht_fdht$short[grep("potential", sigs_ufht_fdht$name.x, value=FALSE)]<-"action potential"
sigs_ufht_fdht$short[grep("catabolic", sigs_ufht_fdht$name.x, value=FALSE)]<-"catabolism"
sigs_ufht_fdht$short[grep("DNA", sigs_ufht_fdht$name.x, value=FALSE)]<-"DNA regulation and repair"
sigs_ufht_fdht$short[grep("transport", sigs_ufht_fdht$name.x, value=FALSE)]<-"ion and amino acid transport"
sigs_ufht_fdht$short[grep("homeosta*", sigs_ufht_fdht$name.x, value=FALSE)]<-"ion and cellular homeostasis"
sigs_ufht_fdht$short[grep("protein", sigs_ufht_fdht$name.x, value=FALSE)]<-"protein folding or localization"
sigs_ufht_fdht$short[grep("secretion", sigs_ufht_fdht$name.x, value=FALSE)]<-"regulation of secretion"
sigs_ufht_fdht$short[grep("stimulus", sigs_ufht_fdht$name.x, value=FALSE)]<-"response to stimulus"
sigs_ufht_fdht$short[grep("RNA", sigs_ufht_fdht$name.x, value=FALSE)]<-"RNA processing or splicing"
sigs_ufht_fdht$short[grep("signal", sigs_ufht_fdht$name.x, value=FALSE)]<-"signaling pathways"
sigs_ufht_fdht$short[is.na(sigs_ufht_fdht$short)]<-"other"

point_colors=c("#A6CEE3" ,"#1F78B4" ,"#B2DF8A" ,"#33A02C" ,"#FB9A99" ,"#E31A1C" ,"#FDBF6F" ,"#FF7F00", "#CAB2D6", "#6A3D9A","gray42")

ggplot(sigs_ufht_fdht, aes(x=delta.rank.x, y=delta.rank.y, color=short))+
  geom_point(size=2.5)+
  scale_color_manual(limits=c("action potential","catabolism","DNA regulation and repair","ion and amino acid transport",
                              "ion and cellular homeostasis","protein folding or localization","regulation of secretion",
                              "response to stimulus","RNA processing or splicing",
                              "signaling pathways","other"),
                     values = point_colors)+
  labs(color="General Biological Process")+
  geom_hline(yintercept =0)+
  geom_vline(xintercept =0)+
  labs(x="Unfed, heated GO Delta Rank",
       y="Fed, heated GO Delta Rank")+
  theme_bw()+
  theme(panel.grid.minor=element_blank(), 
        legend.position = "right",
        legend.background = element_rect(colour="grey"),
        legend.text = element_text(face="bold", colour="black", size=8),
        legend.title = element_text(face="bold", colour="black", size=10),
        legend.key.size = unit(10,"points"),
        axis.text = element_text(face="bold", colour="black", size=8),
        axis.title.x = ggtext::element_markdown(face="bold"),
        axis.title.y = ggtext::element_markdown(face="bold"))

ggsave("fig5.png", width=8, height=5, units="in", dpi=300)
tableS1<-dplyr::select(sigs_ufht_fdht, name.x, short)
colnames(tableS1)<-c("full GO term", "general GO term")
write.csv(tableS1,file="clean_repo/tables/tableS1.csv",quote=F, row.names = F )

# are there any GO terms that don't overlap between them? 
symdiff <- function( x, y) {setdiff( union(x, y), intersect(x, y))}

diffs=symdiff(ufht_goBP$term,fdht_goBP$term)
# 382 GO terms 
ufht_goBP_unique<-ufht_goBP[ufht_goBP$term %in% diffs,]
fdht_goBP_unique<-fdht_goBP[fdht_goBP$term %in% diffs,]

#plot the unique unfed hot GO terms 
ufht_goBP_unique$short<-NA
ufht_goBP_unique$short[grep("cell cycle", ufht_goBP_unique$name, value=FALSE)]<-"cell cycle"
ufht_goBP_unique$short[grep("mitotic", ufht_goBP_unique$name, value=FALSE)]<-"cell cycle"
ufht_goBP_unique$short[grep("replication", ufht_goBP_unique$name, value=FALSE)]<-"cell cycle"
ufht_goBP_unique$short[grep("spindle", ufht_goBP_unique$name, value=FALSE)]<-"cell cycle"
ufht_goBP_unique$short[grep("DNA", ufht_goBP_unique$name, value=FALSE)]<-"cell cycle"
ufht_goBP_unique$short[grep("chromosome", ufht_goBP_unique$name, value=FALSE)]<-"cell cycle"
ufht_goBP_unique$short[grep("chromatin", ufht_goBP_unique$name, value=FALSE)]<-"cell cycle"
ufht_goBP_unique$short[grep("telomer", ufht_goBP_unique$name, value=FALSE)]<-"telomere"
ufht_goBP_unique$short[grep("cellular response", ufht_goBP_unique$name, value=FALSE)]<-"cellular response"
ufht_goBP_unique$short[grep("homeostasis", ufht_goBP_unique$name, value=FALSE)]<-"homeostasis"
ufht_goBP_unique$short[grep("RNA", ufht_goBP_unique$name, value=FALSE)]<-"RNA processing"
ufht_goBP_unique$short[grep("lipid", ufht_goBP_unique$name, value=FALSE)]<-"metabolism"
ufht_goBP_unique$short[grep("transport", ufht_goBP_unique$name, value=FALSE)]<-"metabolism"
ufht_goBP_unique$short[grep("abolic", ufht_goBP_unique$name, value=FALSE)]<-"metabolism"
ufht_goBP_unique$short[grep("response", ufht_goBP_unique$name, value=FALSE)]<-"stress response"
ufht_goBP_unique$short[grep("axon", ufht_goBP_unique$name, value=FALSE)]<-"nervous system"
ufht_goBP_unique$short[grep("nerv", ufht_goBP_unique$name, value=FALSE)]<-"nervous system"
ufht_goBP_unique$short[grep("junction", ufht_goBP_unique$name, value=FALSE)]<-"nervous system"
ufht_goBP_unique$short[grep("neuron", ufht_goBP_unique$name, value=FALSE)]<-"nervous system"
ufht_goBP_unique$short[is.na(ufht_goBP_unique$short)]<-"other"
ufht_goBP_unique$short<-as.factor(ufht_goBP_unique$short)

tableS2<-dplyr::select(ufht_goBP_unique, name, short)
colnames(tableS2)<-c("full GO term", "general GO term")
write.csv(tableS2,file="clean_repo/tables/tableS2.csv",quote=F, row.names = F )




# shorten some of the name for visualization 
ufht_goBP_unique$name<-gsub("regulation","reg",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("positive","+",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("negative","-",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("maintenance","maint.",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("organization","org",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("process","proc",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("cellular","cell",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("microtubules to kinetochore","",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("from RNA polymerase II promoter","",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("involved in cell response to glucose stimulus","",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("via telomere lengthening","",ufht_goBP_unique$name)
ufht_goBP_unique$name<-gsub("mitotic sister","sister",ufht_goBP_unique$name)


empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(ufht_goBP_unique$short), ncol(ufht_goBP_unique)) )
colnames(to_add) <- colnames(ufht_goBP_unique)
to_add$short <- rep(levels(ufht_goBP_unique$short), each=empty_bar)
ufht_goBP_unique2 <- rbind(ufht_goBP_unique, to_add)
ufht_goBP_unique2 <- ufht_goBP_unique2 %>% arrange(short)
ufht_goBP_unique2$id <- seq(1, nrow(ufht_goBP_unique2))


label_data <- ufht_goBP_unique2
number_of_bar <- nrow(ufht_goBP_unique2)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)
label_data$lab_pos<-label_data$delta.rank
label_data$lab_pos[label_data$lab_pos<0]<-10


ggplot(ufht_goBP_unique2, aes(x=as.factor(id), y=delta.rank, fill=as.factor(short))) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-10000,10000)+
  theme_minimal() +
  coord_polar() + 
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +  
  labs(fill="General Biological Process")+
  geom_text(data=label_data, aes(x=id, y=lab_pos+20, label=name, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

ggsave("clean_repo/figures/fig6.png", height=15, width=10, units="in", dpi=300)
