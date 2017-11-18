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
dat=read.table(file="D:/Research/IGP exp_in lab/Analysis/IGP_ParamEstimate_20170324.csv", sep=",",header=T,fill=T)
dat[,2:ncol(dat)] = dat[,2:ncol(dat)]+0.01

dat.l = dat %>%
  melt(id="hr", variable.name="trmt", value.name="den") %>%
  mutate(trmt.all = 0)

dat.l[which(grepl("B_m", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B"
dat.l[which(grepl("C_m", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C"
dat.l[which(grepl("B_250", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B_250"
dat.l[which(grepl("C_250", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C_250"
dat.l[which(grepl("B_10", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="B_10"
dat.l[which(grepl("C_10", as.character(dat.l[,"trmt"]))==TRUE),"trmt.all"]="C_10"

dat.summary = ddply(dat.l, c("trmt.all", "hr"), summarize,
                    N = length(den),
                    mean = log(mean(den)),
                    sd = sd(log(den)),
                    se = sd/sqrt(N)
                    )

ggplot(dat.summary, aes(x=hr, y=mean, colour=trmt.all)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1) +
  geom_line(aes(group=trmt.all)) +
  geom_point(size=4) + 
  xlab("hr") + 
  ylab("Log[density(ind./mL)]") + 
  scale_colour_manual("", 
                     values=c("#000033", "#000099", "#0066FF", "#003300", "#336600", "#66CC33")) + 
                     #labels=c("Blepharisma_W/ Colpidium", "Colpidium_W/ Blephsrisma", "Blepharisma", "Colpidium"))+ 
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
        legend.key=element_rect(color="white", fill="white"), 
        legend.title=element_text(size=14),
        legend.text=element_text(size=14)
        )
  

screen35=read.table(paste(paste(path, "35screen",sep="//"), ".csv", sep=""), sep=",",header=F,fill=T)
screen35=as.data.frame(screen35)
screen35=cbind(screen35, screen35[,1]*20/(screen35[,2]*20))
  colnames(screen35)=c("Colpidium", "Blepharisma", "treat", "BA", "ratio")
  
ggplot(screen35, aes(x = as.factor(treat), y = ratio, fill=BA)) +
  geom_boxplot() +
  labs(x="% volume screened", y="# of Colpidium per Blepharisma") + 
  scale_x_discrete(labels=c("20%","40%","60%", "80%"))+
  scale_fill_manual(values = c("green", "darkgreen"))+
  guides(fill=guide_legend(title=NULL, reverse = T))+
  theme_bw() +
  theme(panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.margin=unit(0, "cm"),
        legend.key=element_rect(color="white", fill="white"), 
        legend.title=element_text("", size=16),
        legend.text=element_text(size=16),
        axis.line=element_line(colour="black"),
        axis.title=element_text(size=18, margin=margin(0,20,0,0)),
        axis.text=element_text(size=18))

t.test(screen35[which(screen35$BA=="before"&screen35$treat=="0.2"),"ratio"], 
       screen35[which(screen35$BA=="after"&screen35$treat=="0.2"),"ratio"])
