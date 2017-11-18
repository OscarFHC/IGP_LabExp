library(MASS)
library(faraway)
library(lme4)
library(nlme)
library(languageR)
library(AICcmodavg)

library(ggplot2)
library(plyr)
library(magrittr)
library(reshape2)
dat=read.table(file="D:/Research/IGP_LabExp/Analysis/IGP_II_20170419.csv", sep=",",header=T,fill=T)
dat[,2:ncol(dat)] = dat[,2:ncol(dat)]+0.01

dat.l = dat %>%
  melt(id="hr", variable.name="trmt", value.name="den") %>%
  mutate(trmt.all = 0)

dat.l[which(grepl("B00", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B00"
dat.l[which(grepl("C00", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C00"
dat.l[which(grepl("B02", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B02"
dat.l[which(grepl("C02", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C02"
dat.l[which(grepl("B04", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B04"
dat.l[which(grepl("C04", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C04"
dat.l[which(grepl("B06", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B06"
dat.l[which(grepl("C06", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C06"
dat.l[which(grepl("B08", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B08"
dat.l[which(grepl("C08", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C08"
dat.l[which(grepl("B10", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B10"
dat.l[which(grepl("C10", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C10"

dat.summary = ddply(dat.l, c("trmt.all", "hr"), summarize,
                    rep = length(den),
                    mean = log(mean(den)),
                    sd = sd(log(den)),
                    se = sd/sqrt(rep)
                    )
B_sum = dat.summary[which(grepl("B", as.character(dat.summary[,"trmt.all"]))==TRUE),]
C_sum = dat.summary[which(grepl("C", as.character(dat.summary[,"trmt.all"]))==TRUE),]

B_sum[,"trmt.all"] = relevel(as.factor(B_sum[,"trmt.all"]), ref=c("B100"))
C_sum[,"trmt.all"] = relevel(as.factor(C_sum[,"trmt.all"]), ref=c("C100"))

ggplot(B_sum, aes(x=hr, y=mean, colour=trmt.all)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_line(aes(group=trmt.all)) +
  geom_point(size=4) + 
  xlab("hr") + 
  ylab("Log[density(ind./mL)]") + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
        legend.key=element_rect(color="white", fill="white"), 
        legend.title=element_text(size=14),
        legend.text=element_text(size=14)
        )+
  scale_colour_hue(l=45)

mod_B =  kruskal.test(mean~trmt.all, data=B_sum)
mod_B
dunn.test(B_sum[,"mean"], 
          B_sum[,"trmt.all"], method="bonferroni")

ggplot(C_sum, aes(x=hr, y=mean, colour=trmt.all)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_line(aes(group=trmt.all)) +
  geom_point(size=4) + 
  xlab("hr") + 
  ylab("Log[density(ind./mL)]") + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
        legend.key=element_rect(color="white", fill="white"), 
        legend.title=element_text(size=14),
        legend.text=element_text(size=14)
  )+
  scale_colour_hue(l=45)

mod_C =  kruskal.test(mean~trmt.all, data=C_sum)
mod_C
dunn.test(C_sum[,"mean"], 
          C_sum[,"trmt.all"], method="bonferroni")

# B_sum_box = B_sum[] %>%
#  dcast(trmt.all~hr, fun.aggregate=mean, value.var=c("mean"), subset=.(hr >= \250))


B_stat = B_sum %>%
  subset(hr>400 & hr<600) %>%
  ggplot(aes(trmt.all, mean)) + 
  geom_boxplot() + 
  xlab("treatment") + 
  ylab("Log[density(ind./mL)]") + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
        legend.key=element_rect(color="white", fill="white"), 
        legend.title=element_text(size=14),
        legend.text=element_text(size=14)
  ) 


C_sum[which(C_sum[,"hr"]>=250),] %>%
  ggplot(aes(trmt.all, mean)) + 
  geom_boxplot() + 
  xlab("treatment") + 
  ylab("Log[density(ind./mL)]") + 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
        legend.key=element_rect(color="white", fill="white"), 
        legend.title=element_text(size=14),
        legend.text=element_text(size=14)
  )  
 
