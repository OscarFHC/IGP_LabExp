####### Loading required library #####

if (!require(ggplot2)) {
  install.packages("ggplot2", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(ggplot2)
}else{library(ggplot2)}

if (!require(cowplot)) {
  install.packages("cowplot", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(cowplot)
}else{library(cowplot)}

if (!require(tidyverse)) {
  install.packages("tidyverse", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(tidyverse)
}else{library(tidyverse)}

if (!require(viridis)) {
  install.packages("viridis", dependencies=TRUE, repos = 'http://cran.us.r-project.org')
  library(viridis)
}else{library(viridis)}

if (!require(RColorBrewer)) {
  install.packages("RColorBrewer", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
  library(RColorBrewer)
}else{library(RColorBrewer)}

if (!require(gridExtra)) {
  install.packages("gridExtra", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
  library(gridExtra)
}else{library(gridExtra)}

####### Loading required library #####

####### Loading protozoa dataset and organizing data for further plotting and analyses #####
dat <- read.table(file="https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/IGP_II_20170419_1.csv", sep = ",", header = TRUE, fill = TRUE)
dat[,2:ncol(dat)] <- dat[, 2:ncol(dat)] + 0.01

dat_l <- dat %>% 
  gather(key = trmt, value = den, -hr) %>% 
  mutate(trmt_all = 0)

dat_l[which(grepl("B00", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="B00"
dat_l[which(grepl("C00", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="C00"
dat_l[which(grepl("B02", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="B02"
dat_l[which(grepl("C02", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="C02"
dat_l[which(grepl("B04", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="B04"
dat_l[which(grepl("C04", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="C04"
dat_l[which(grepl("B06", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="B06"
dat_l[which(grepl("C06", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="C06"
dat_l[which(grepl("B08", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="B08"
dat_l[which(grepl("C08", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="C08"
dat_l[which(grepl("B10", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="B10"
dat_l[which(grepl("C10", as.character(dat_l[,"trmt"]))==TRUE), "trmt_all"]="C10"

NonParam_log <- function(den){ # This function log transforms the data when calculatingnon-parametric SE
  set.seed(1032)
  means = replicate(5000, mean(log(sample(den, size=length(den), replace=TRUE))))
  sd(means)
}

NonParam <- function(den){ # This function does not log-transform data before doing permutation
  set.seed(1032)
  means = replicate(5000, mean(sample(den, size=length(den), replace=TRUE)))
  sd(means)
}

dat_sum <- dat_l %>% 
  group_by(trmt_all, hr) %>% 
  summarize(
    rep = length(den),
    avg = mean(den),
    avg_log = mean(log(den)),
    sd = sd(den),
    sd_log = sd(log(den)),
    se = sd/sqrt(rep),
    se_log = sd_log/sqrt(rep),
    se_permu = NonParam(den),
    se_permu_log = NonParam_log(den)
  ) %>%
  ungroup()

tp_start <- 2
tp_end <- length(unique(dat_sum[["hr"]]))-4
timepoint <- unique(dat_sum[["hr"]])[tp_start:tp_end]

tf <- matrix(0, 12, length(timepoint))
for (i in 1:length(timepoint)){
  mean_t <- dat_sum %>%
    subset(hr %in% unique(dat_sum[["hr"]])[(i+1):(i+5)]) %>%
    group_by(trmt_all) %>%
    summarize(
      mean = mean(avg)
    )
  mean_tpre <- dat_sum%>%
    subset(hr %in% unique(dat_sum[["hr"]])[i:(i+4)]) %>%
    group_by(trmt_all) %>%
    summarize(
      mean = mean(avg)
    )
  tf[,i] <- mean_t$"mean" / mean_tpre$"mean"
}

NetDenChange <- colSums(abs(tf-1))
best <- which(NetDenChange==min(NetDenChange)) + 1

B_sum_00 <- dat_sum %>% 
  filter(hr %in% timepoint[best:(best+4)]) %>%
  filter(trmt_all == "B00")

C_sum_00 <- dat_sum %>% 
  filter(hr %in% timepoint[best:(best+4)]) %>%
  filter(trmt_all == "C00")

B_cal <- dat_sum %>%
  filter(hr %in% timepoint[best:(best+4)] & substr(trmt_all, 1, 1) == "B") %>%
  group_by(trmt_all) %>%
  summarize(
    rep = length(avg),
    avg_trmt = mean(avg/mean(B_sum_00[["avg"]])),
    sd_trmt = sd(avg/mean(B_sum_00[["avg"]])),
    se_trmt = sd_trmt/sqrt(rep),
    se_permu = NonParam(avg/mean(B_sum_00[["avg"]]))
  ) %>%
  mutate(alpha = seq(0, 1, by = 0.2))

C_cal <- dat_sum %>%
  filter(hr %in% timepoint[best:(best+4)] & substr(trmt_all, 1, 1) == "C") %>%
  group_by(trmt_all) %>%
  summarize(
    rep = length(avg),
    avg_trmt = mean(avg/mean(C_sum_00[["avg"]])),
    sd_trmt = sd(avg/mean(C_sum_00[["avg"]])),
    se_trmt = sd_trmt/sqrt(rep),
    se_permu = NonParam(avg/mean(C_sum_00[["avg"]]))
  ) %>%
  mutate(alpha = seq(0, 1, by = 0.2))

dat_cal <- rbind(B_cal, C_cal)
####### Loading protozoa dataset and organizing data for further plotting and analyses #####

####### Loading bacteria dataset and organizing data for further plotting and analyses #####
Bac_raw = read.table(file="https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/IGP_II_20170419_bac2.csv",
                     fill=TRUE, header=TRUE, sep=",")

Bac = Bac_raw %>%
  mutate(time = rep(rep(c("T0", "Tend"), each=nrow(Bac_raw)/8), each=4)) %>%
  subset(Gate.Name %in% c("culture population")) %>%
  mutate(trmt = rep(paste0(
    c(rep(c("00_", "02_", "04_", "06_"), each=3), 
      rep(("08_"), each=4),
      rep(c("10_", "ctrl_"), each=5)),
    c(rep(seq(from=1, to=3), 4),
      rep(seq(from=1, to=4), 1),
      rep(seq(from=1, to=5), 2))
  ), 2)
  ) %>%
  mutate(trmt_all = rep(
    c(rep(c("Bac00", "Bac02", "Bac04", "Bac06"), each=3), 
      rep(c("Bac08"), each=4),
      rep(c("Bac10", "Bacctrl"), each=5)),2)
  ) %>%
  subset(select=c("time", "trmt_all", "trmt", "Gate.Name", "Concentration")) 
colnames(Bac) = c("time", "trmt_all", "trmt", "pop", "den")

Cons_ctrl = mean(Bac[which(Bac[,"trmt_all"] == "Bacctrl" & Bac[,"time"] == "Tend"), "den"])
Cons_00 = mean(Bac[which(Bac[,"trmt_all"] == "Bac00" & Bac[,"time"] == "Tend"), "den"])
Eqm_00 = mean(Bac[which(Bac[,"trmt_all"] == "Bac00" & Bac[,"time"] == "Tend"), "den"])

Bac_pop <- Bac %>%
  subset(select = -pop) %>%
  spread(key = time, value = den) %>%
  mutate(Consumption = Cons_ctrl - Tend) %>%
  mutate(Change_Cons = Consumption / Cons_00) %>%
  mutate(Change_Eqm = Tend/Eqm_00)

Bac_cal <- Bac_pop %>%
  group_by(trmt_all) %>%
  summarize(
    rep = length(Change_Eqm),
    avg_trmt = mean(Change_Eqm),
    sd_trmt = sd(Change_Eqm),
    se_trmt = sd_trmt/sqrt(rep),
    se_permu = NonParam(Change_Eqm)
  ) %>%
  mutate(alpha = c(seq(0, 1, by = 0.2), 0))

dat_cal <- rbind(dat_cal, Bac_cal)

####### Loading bacteria dataset and organizing data for further plotting and analyses #####

###########################################################################################
##### Finding the best combination of efficiency (e2) and resource partitioning (p) #######
###########################################################################################

##### Estimating Attack rate (A) and handling time (h) from experiment ####################
BC_Type2=read.table(file="https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/BonC_Type2_20171112.csv", 
                    sep = ",", header = T, fill = T)

BC_nls = BC_Type2 %>% 
  group_by(Trmt, rep) %>%
  summarize(
    InitRatio = mean(T0),
    End_W = mean(T4_W),
    End_WO = mean(T4_WO),
    Cons = (End_WO-End_W)
  )

Type2_mod = nls(Cons ~ (s2 * InitRatio)/(1 + (h2 * s2 * InitRatio)), 
                start = list(s2 = 0.4, h2 = 0.34), data = BC_nls)

Type1_mod = lm(Cons ~ InitRatio, data = BC_nls)
summary(Type2_mod)
anova(Type1_mod, Type2_mod)

SST <- sum((mean(BC_nls$Cons) - BC_nls$Cons)^2) 
SSR <- sum((predict(Type2_mod) - BC_nls$Cons)^2)
R2 <- 1-(SSR/SST)
##### Estimating Attack rate (A) and handling time (h) from experiment ####################

##### Loading other parameter value estimation ############################################
mod.lst <- expand.grid(e2 = seq(0, 0.5, 0.01), partition = seq(0, 0.15, 0.01))

Bac_pop_mod <- Bac_pop[-grep("ctrl", Bac_pop[,"trmt_all"]),] %>%
  mutate(trmt = as.numeric(c(rep(c("0", "0.2", "0.4", "0.6"), each=3), rep(("0.8"), each=4), rep(c("1.0"), each=5))),
         trmt2 = trmt^2)

for (i in 1:nrow(mod.lst)){
  temp <- read.table(file = paste0("https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/Type2ModResults/Type2_e2_", 
                                   mod.lst[i, "e2"]*10, "_p_", mod.lst[i,"partition"]*100, "..csv"), sep=",", header=FALSE)
    colnames(temp) <- c("alpha", "R1", "R2", "Cons", "Pred")
  temp = temp %>%
    mutate(Res = (R1+R2) / (temp[1,"R1"] + temp[1,"R2"]),
           Cons = Cons / (temp[1,"Cons"]),
           Pred = Pred / (temp[1,"Pred"])) %>%
    subset(alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))

  mod.lst[i, "RSS_Pred"] <- sum((B_cal[["avg_trmt"]] - temp[,"Pred"])^2)/
                            sum((B_cal[["avg_trmt"]] - mean(temp[,"Pred"]))^2)
  mod.lst[i, "RSS_Prey"] <- sum((C_cal[["avg_trmt"]] - temp[,"Cons"])^2)/
                            sum((C_cal[["avg_trmt"]] - mean(temp[,"Cons"]))^2)
  
  tempdat_Bac <- data.frame(alpha = c(rep(c(0, 0.2, 0.4, 0.6), each = 3), rep(0.8, 4), rep(1.0, 5))) %>%
    inner_join(temp, by = "alpha")
  mod.lst[i, "RSS_Bac"] <- sum((Bac_pop_mod[["Change_Eqm"]] - tempdat_Bac[,"Res"])^2)/
    sum((Bac_pop_mod[["Change_Eqm"]] - mean(Bac_pop_mod[["Change_Eqm"]]))^2)
  # 
  # tempdat_BC <- temp %>%
  #   subset(alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))
  # 
  # mod.lst[i, "RSS_Pred"] <- sum((B_pop_mod[["Change_Eqm"]] - rep(tempdat_BC[,"Pred"], each = 5))^2)/
  #                           sum((B_pop_mod[["Change_Eqm"]] - mean(B_pop_mod[["Change_Eqm"]]))^2)
  # mod.lst[i, "RSS_Prey"] <- sum((C_pop_mod[,"Change_Eqm"] - rep(tempdat_BC[,"Cons"], each = 5))^2)/
  #                           sum((C_pop_mod[,"Change_Eqm"] - mean(C_pop_mod[["Change_Eqm"]]))^2)
}

mod.lst <- mod.lst %>%
  mutate(RSS_tot = RSS_Bac + RSS_Prey + RSS_Pred) 

e2 <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["e2"]]
partition <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["partition"]]
##### Loading other parameter value estimation ############################################

###########################################################################################
##### Finding the best combination of efficiency (e2) and resource partitioning (p) #######
###########################################################################################

###########################################################################################
##### Loading the model results from Type II model ########################################
###########################################################################################
dat_mod_II = read.table(file = 
                        paste0("https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/Type2ModResults/Type2_e2_",
                               e2*10, "_p_", partition*100, "..csv"),
                        sep = ",", header = FALSE, fill = TRUE)
colnames(dat_mod_II) <- c("alpha", "R1", "R2", "Cons", "Pred")
dat_mod_II = dat_mod_II %>%
  mutate(Res = (R1+R2) / (dat_mod_II[1,"R1"] + dat_mod_II[1,"R2"]),
         Cons = Cons / (dat_mod_II[1,"Cons"]),
         Pred = Pred / (dat_mod_II[1,"Pred"])) 

# Type I R2 for IG predator
B_II_R2 <- cbind(as.data.frame(B_cal[,1:6]), subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
B_II_Cor <- cor(B_II_R2$avg_trmt, B_II_R2$Pred)
SST_II_pred = sum((B_II_R2[,"avg_trmt"] - mean(B_II_R2[,"avg_trmt"]))^2)
SSR_II_pred = sum((B_II_R2[,"avg_trmt"] - B_II_R2[,"Pred"])^2)
R2_II_pred = ifelse((1-(SSR_II_pred/SST_II_pred)) > 0, 1-(SSR_II_pred/SST_II_pred), 0)

# Type I R2 for IG prey
C_II_R2 <- cbind(as.data.frame(C_cal[,1:6]), subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
C_II_Cor <- cor(C_II_R2$avg_trmt, C_II_R2$Cons)
SST_II_prey = sum((C_II_R2[,"avg_trmt"] - mean(C_II_R2[,"avg_trmt"]))^2)
SSR_II_prey = sum((C_II_R2[,"avg_trmt"] - C_II_R2[,"Cons"])^2)
R2_II_prey = ifelse((1-(SSR_II_prey/SST_II_prey)) > 0, (1-(SSR_II_prey/SST_II_prey)), 0)

# Type I R2 for Bacteria
Bac_II_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
Bac_II_Cor <- cor(Bac_II_R2$avg_trmt, Bac_II_R2$Res)
SST_II_bac = sum((Bac_II_R2[,"avg_trmt"] - mean(Bac_II_R2[,"avg_trmt"]))^2)
SSR_II_bac = sum((Bac_II_R2[,"avg_trmt"] - Bac_II_R2[,"Res"])^2)
R2_II_bac = ifelse((1-(SSR_II_bac/SST_II_bac)) > 0, (1-(SSR_II_bac/SST_II_bac)), 0)
###########################################################################################
##### Loading the model results from Type II model ########################################
###########################################################################################

###########################################################################################
##### Loading the model results from Type I model #########################################
###########################################################################################
dat_mod_I = read.table(file="https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/Mod_Type1.csv", 
                       sep=",", header=TRUE, fill=TRUE)
dat_mod_I = dat_mod_I %>%
  mutate(Res = (R1+R2) / (dat_mod_I[1,"R1"] + dat_mod_I[1,"R2"]),
         Cons = C / (dat_mod_I[1,"C"]),
         Pred = P / (dat_mod_I[1,"P"]))

# Type I R2 for IG predator
B_I_R2 <- cbind(as.data.frame(B_cal[1:6,1:6]), subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
B_I_Cor <- cor(B_I_R2$avg_trmt, B_I_R2$Pred)
SST_I_pred = sum((B_I_R2[,"avg_trmt"] - mean(B_I_R2[,"avg_trmt"]))^2)
SSR_I_pred = sum((B_I_R2[,"avg_trmt"] - B_I_R2[,"Pred"])^2)
R2_I_pred = ifelse((1-(SSR_I_pred/SST_I_pred)) > 0, 1-(SSR_I_pred/SST_I_pred), 0)

# Type I R2 for IG prey
C_I_R2 <- cbind(as.data.frame(C_cal[1:6,1:6]), subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
C_I_Cor <- cor(C_I_R2$avg_trmt, C_I_R2$Cons)
SST_I_prey = sum((C_I_R2[,"avg_trmt"] - mean(C_I_R2[,"avg_trmt"]))^2)
SSR_I_prey = sum((C_I_R2[,"avg_trmt"] - C_I_R2[,"Cons"])^2)
R2_I_prey = ifelse((1-(SSR_I_prey/SST_I_prey)) > 0, (1-(SSR_I_prey/SST_I_prey)), 0)

# Type I R2 for Bacteria
Bac_I_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
Bac_I_Cor <- cor(Bac_I_R2$avg_trmt, Bac_I_R2$Res)
SST_I_bac = sum((Bac_I_R2[,"avg_trmt"] - mean(Bac_I_R2[,"avg_trmt"]))^2)
SSR_I_bac = sum((Bac_I_R2[,"avg_trmt"] - Bac_I_R2[,"Res"])^2)
R2_I_bac = ifelse((1-(SSR_I_bac/SST_I_bac)) > 0, (1-(SSR_I_bac/SST_I_bac)), 0)

###########################################################################################
##### Loading the model results from Type I model #########################################
###########################################################################################

############################################################################################
##### for internal hump (Fig 1) ############################################################
############################################################################################
if (!require(vegan)) {
  install.packages("vegan", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
  library(vegan)
}else{library(vegan)}

dat_emp_I = read.table(file = 
                         paste0("https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/Empirical_Type1.csv"),
                       sep = ",", header = FALSE, fill = TRUE)
colnames(dat_emp_I) <- c("alpha", "R1", "R2", "Cons", "Pred")
dat_emp_I = dat_emp_I %>%
  mutate(Res = (R1+R2) / (dat_emp_I[1,"R1"] + dat_emp_I[1,"R2"]),
         Cons = Cons / (dat_emp_I[1,"Cons"]),
         Pred = Pred / (dat_emp_I[1,"Pred"])) 

Bac_pop_mod <- Bac_pop[-grep("ctrl", Bac_pop[,"trmt_all"]),] %>%
  mutate(trmt = as.numeric(c(rep(c("0", "0.2", "0.4", "0.6"), each=3), rep(("0.8"), each=4), rep(c("1.0"), each=5))),
         trmt2 = trmt^2)
mod.t1 <- lm(Change_Eqm ~ trmt, data = Bac_pop_mod)
mod.t2 <- lm(Change_Eqm ~ trmt + trmt2, data = Bac_pop_mod)

anova(mod.t1, mod.t2)
anova(mod.t1, mod.t2, test = "Chisq")

summary(mod.t2)

# II <- c(rep(subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))[1:4,"Res"], each=3), 
#         rep(subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))[5,"Res"], each=4),
#         rep(subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))[6,"Res"], each=5))
# I <- c(rep(subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))[1:4,"Res"], each=3), 
#        rep(subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))[5,"Res"], each=4),
#        rep(subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0))[6,"Res"], each=5))
# SSR_II <- sum((Bac_pop_mod[,"Change_Eqm"] - II)^2)
# SSR_I <- sum((Bac_pop_mod[,"Change_Eqm"] - I)^2)
# SST <- sum((Bac_pop_mod[,"Change_Eqm"] - mean(Bac_pop_mod[,"Change_Eqm"]))^2)
# 1 - SSR_II/SST
##### Plotting bacteria density at the steady state along with a quadratic function                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
Bacplot_mod <- ggplot() +
  geom_point(data = Bac_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = Bac_pop_mod, aes(x = trmt, y = Change_Eqm), method = 'lm', formula = y ~ poly(x, 2), se = TRUE, size = 1, color = "black") +
  #geom_smooth(data = Bac_II_R2, aes(x = alpha, y = Res), method = 'lm', formula = y ~ poly(x, 2), se = TRUE, size = 1, color = "black") +
  geom_line(data = dat_mod_I, aes(x = alpha, y = Res), size = 1, linetype = "longdash") +
  geom_smooth(data = dat_emp_I[dat_emp_I$Res<2, ], aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, linetype = 3, color = "black") +
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", expression(bold("80%*")), expression(bold("100%*")))) + 
  scale_y_continuous(limits = c(0.5, 2), expand = c(0, 0), breaks = seq(0.5, 2, 0.5)) + 
  labs(x = "", y = "")+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0, r = 0.4, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20))
Bacplot_mod_up <- ggplot() +
  geom_smooth(data = dat_emp_I[dat_emp_I$Res>2, ], aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, linetype = 3, color = "black") +
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", expression(bold("80%*")), expression(bold("100%*")))) + 
  scale_y_continuous(limits = c(2, 6), expand = c(0, 0), breaks = seq(3.5, 6.0, 1.5)) + 
  labs(x = "", y = "")+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 0.4, b = -0.35, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text = element_text(size = 20), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

Fig1 <- plot_grid(Bacplot_mod_up, Bacplot_mod, nrow = 2, rel_heights = c(1, 2)) + 
  draw_label("Availability of IG prey (IGP strength)", 
             x = 0.54, y = 0.02, size = 20) +
  draw_label(expression("error bars represent standard error of the mean"), 
             x = 0.87, y = 0.01, size = 12) + 
  draw_label(expression(atop("Bacteria density at the steady state")), 
             angle = 90, x = 0.023, y = 0.5, size = 20)

# ggsave(filename = "D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/Fig2_Bac_IGP.eps",
#        plot = Fig1, width = 36, height = 24, units = c("cm"), dpi = 600)
############################################################################################
##### for internal hump (Fig. 1) ###########################################################
############################################################################################

############################################################################################
####### Fig. S1 protozoa population dynamics ###############################################
############################################################################################
NetDenChange <- colSums(abs(tf-1))
i <- which(NetDenChange==min(NetDenChange)) + 1

B_plot <- dat_sum %>%
  subset(substr(trmt_all, 1, 1) == "B") %>% 
  ggplot(aes(x = hr, y = avg_log, colour = trmt_all)) + 
  geom_errorbar(aes(ymin = avg_log - se_permu_log, ymax = avg_log + se_permu_log), width = 0.1) +
  geom_line() +
  geom_vline(xintercept = c(timepoint[i], timepoint[i+4]), linetype = 2) +
  geom_point(size = 4) +
  scale_x_continuous(expand = c(0, 0), limits = c(-10, 620)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-6, 6.5)) +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
                      guide = guide_legend(byrow=TRUE, title = expression("IGP \nstrength"))) + 
  labs(x = "", y = "", title = "") + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.3, r = 0.4, b = 0.1, l = 1.4, "cm"),
        legend.justification = c("right", "top"),
        legend.margin = margin(t = 0, r = 0.4, b = 0, l = 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size=14))

C_plot <- dat_sum %>%
  subset(substr(trmt_all, 1, 1) == "C") %>% 
  ggplot(aes(x = hr, y = avg_log, colour = trmt_all)) + 
  geom_errorbar(aes(ymin = avg_log - se_permu_log, ymax = avg_log + se_permu_log), width = 0.1) +
  geom_line() +
  geom_vline(xintercept = c(timepoint[i], timepoint[i+4]), linetype = 2) +
  geom_point(size = 4) +
  scale_x_continuous(expand = c(0, 0), limits = c(-10, 620)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-6, 6.5)) +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
                      guide = guide_legend(byrow = TRUE, title = expression("IGP \nstrength"))) + 
  labs(x = "", y = "", title = "") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0, r = 0.4, b = 0.6, l = 1.4, "cm"),
        legend.justification = c("right", "top"),
        legend.margin = margin(t = 0, r = 0.4, b = 0, l = 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 14))

FigS1 <- plot_grid(B_plot + theme(legend.position="none"), 
                   C_plot + theme(legend.position="none"), 
                   labels = c("a.", "b."), ncol = 1, align = 'v')
shared_legend <- get_legend(B_plot + theme(legend.position = c(1.0, 0.58)))
FigS1 <- plot_grid(FigS1, shared_legend, nrow = 1, rel_widths = c(3, 0.25))
FigS1 <- ggdraw(FigS1) + 
  draw_label("Hours", 
             x = 0.52, y = 0.02, size = 18) + 
  draw_label("Protozoa density (Log[0.01 + ind./mL])", angle = 90, 
             x = 0.035, y = 0.5, size = 18) + 
  draw_label(expression("IG prey (" * italic("Colpidium") * ")"), 
             x = 0.52, y = 0.5, size = 14) + 
  draw_label(expression("IG predator (" * italic("Blespharisma") * ")"), 
             x = 0.52, y = 0.98, size = 14) + 
  draw_label(expression("error bars represent standard error of the mean"), 
             x = 0.85, y = 0.01, size = 12)

# ggsave(filename = "D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/FigS1_PopDyna.pdf", 
#        plot = FigS1, width = 35, height = 24, units = c("cm"), dpi = 600)
############################################################################################
####### Fig. S1 protozoa population dynamics ###############################################
############################################################################################

############################################################################################
####### Fig. 2 Showing Type II IGP #########################################################
############################################################################################
##### Estimating Attack rate (A) and handling time (h) from experiment ####################
BC_Type2=read.table(file="https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/BonC_Type2_20171112.csv", 
                    sep = ",", header = T, fill = T)

BC_nls = BC_Type2 %>% 
  group_by(Trmt, rep) %>%
  summarize(
    InitRatio = mean(T0),
    End_W = mean(T4_W),
    End_WO = mean(T4_WO),
    Cons = (End_WO-End_W)
  )

Type2_mod = nls(Cons ~ (s2 * InitRatio)/(1 + (h2 * s2 * InitRatio)), 
                start = list(s2 = 0.4, h2 = 0.34), data = BC_nls)
##### Estimating Attack rate (A) and handling time (h) from experiment ####################

##### Loading other parameter value estimation ############################################
mod.lst <- expand.grid(e2 = seq(0, 0.5, 0.01), partition = seq(0, 0.15, 0.01))

for (i in 1:nrow(mod.lst)){
  tempdat <- read.table(file = paste0("https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/Type2ModResults/Type2_e2_", 
                                      mod.lst[i, "e2"]*10, "_p_", mod.lst[i,"partition"]*100, "..csv"), sep=",", header=FALSE)
  colnames(tempdat) <- c("alpha", "P1", "P2", "C", "B")
  
  PredSSR <- cbind(B_cal, subset(tempdat, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
  PredSSR[,"B"] <- PredSSR[,"B"]/PredSSR[1,"B"]
  mod.lst[i, "RSS_Pred"] <- sum((PredSSR[,"avg_trmt"] - PredSSR[,"B"])^2)/
    sum((PredSSR[,"avg_trmt"] - mean(PredSSR[,"avg_trmt"]))^2)
  
  PreySSR <- cbind(C_cal, subset(tempdat, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
  PreySSR[,"C"] <- PreySSR[,"C"]/PreySSR[1,"C"]
  mod.lst[i, "RSS_Prey"] <- sum((PreySSR[,"avg_trmt"] - PreySSR[,"C"])^2)/
    sum((PreySSR[,"avg_trmt"] - mean(PreySSR[,"avg_trmt"]))^2)
  
  ResSSR <- cbind(Bac_cal[-nrow(Bac_cal),], subset(tempdat, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
  ResSSR[,"Bac"] <- (ResSSR[,"P1"] + ResSSR[,"P2"]) / (ResSSR[1,"P1"] + ResSSR[1,"P2"])
  mod.lst[i, "RSS_Bac"] <- sum((ResSSR[,"avg_trmt"] - ResSSR[,"Bac"])^2)/
    sum((ResSSR[,"avg_trmt"] - mean(ResSSR[,"avg_trmt"]))^2)
}

mod.lst <- mod.lst %>%
  mutate(RSS_tot = RSS_Bac + RSS_Prey + RSS_Pred) 

e2 <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["e2"]]
partition <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["partition"]]
##### Loading other parameter value estimation ############################################

##### Plotting A and h estimation ##################################################
TyprII_Param_hA <- BC_nls %>%
  ggplot(aes(x = InitRatio)) +
  geom_jitter(aes(y = Cons), size = 5, height = 0.05) + 
  geom_line(aes(y = predict(Type2_mod)), linetype = 2, size = 2)+
  labs(x = "IG prey density per IG predator (ind./ml)", 
       y = expression(atop("Number of IG prey consumed", "(ind./ " ~ mL%.%day ~ ")")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0, r = 0.4, b = 0.6, l = 1.5, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.x = element_text(size = 20, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)))

##### Plotting A and h estimation ##################################################

##### Plotting other parameter estimation ##########################################
TypeII_Param <- mod.lst %>%
  subset(RSS_Pred<1 & RSS_Prey<1 & RSS_Bac<1) %>%
  ggplot(aes(x = e2, y = 1-partition, fill = RSS_tot)) +
  geom_tile(color = "white", size = 0.5) +
  coord_fixed(ratio = 2, expand = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_viridis() + 
  labs(x = "Assimilation efficiency of IG prey", 
       y = expression(atop("Degree of", "resource partitioning")), 
       fill = expression("Normalized \ntotal residual \nsum of squared")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0, r = 0.4, b = 0.6, l = 1.4, "cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.justification = c("right", "center"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.x = element_text(size = 20, margin = margin(t = 12, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 4, b = 0, l = 0)))
##### Plotting other parameter estimation ##########################################

##### Combining two plots ##########################################################
Fig2_Param <- plot_grid(TyprII_Param_hA, TypeII_Param, 
                        labels = c("a.", "b."), ncol = 1)

ggsave(filename="D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/Fig3_TypeII_Param.pdf",
       plot = Fig2_Param, width = 36, height = 24, units = c("cm"), dpi = 600)
##### Combining two plots ##########################################################
##########################################################################################
##### Fig. 2 Showing Type II IGP #########################################################
##########################################################################################

###########################################################################################
##### Finding the best combination of Accack rate (A) and handling time (h) ###############
###########################################################################################
e2 <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["e2"]]
partition <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["partition"]]

dat_mod_II = read.table(file = 
                        paste0("https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/Type2ModResults/Type2_e2_",
                               e2*10, "_p_", partition*100, "..csv"),
                        sep = ",", header = FALSE, fill = TRUE)
colnames(dat_mod_II) <- c("alpha", "R1", "R2", "Cons", "Pred")
dat_mod_II = dat_mod_II %>%
  mutate(Res = (R1+R2) / (dat_mod_II[1,"R1"] + dat_mod_II[1,"R2"]),
         Cons = Cons / (dat_mod_II[1,"Cons"]),
         Pred = Pred / (dat_mod_II[1,"Pred"]))

# Type I R2 for IG predator
B_II_R2 <- cbind(as.data.frame(B_cal[,1:6]), subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_pred = sum((B_II_R2[,"avg_trmt"] - mean(B_II_R2[,"avg_trmt"]))^2)
SSR_II_pred = sum((B_II_R2[,"avg_trmt"] - B_II_R2[,"Pred"])^2)
R2_II_pred = ifelse((1-(SSR_II_pred/SST_II_pred)) > 0, 1-(SSR_II_pred/SST_II_pred), 0)

# Type I R2 for IG prey
C_II_R2 <- cbind(as.data.frame(C_cal[,1:6]), subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_prey = sum((C_II_R2[,"avg_trmt"] - mean(C_II_R2[,"avg_trmt"]))^2)
SSR_II_prey = sum((C_II_R2[,"avg_trmt"] - C_II_R2[,"Cons"])^2)
R2_II_prey = ifelse((1-(SSR_II_prey/SST_II_prey)) > 0, (1-(SSR_II_prey/SST_II_prey)), 0)

# Type I R2 for Bacteria
Bac_II_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(dat_mod_II, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_bac = sum((Bac_II_R2[,"avg_trmt"] - mean(Bac_II_R2[,"avg_trmt"]))^2)
SSR_II_bac = sum((Bac_II_R2[,"avg_trmt"] - Bac_II_R2[,"Res"])^2)
R2_II_bac = ifelse((1-(SSR_II_bac/SST_II_bac)) > 0, (1-(SSR_II_bac/SST_II_bac)), 0)

###########################################################################################
##### Finding the best combination of Accack rate (A) and handling time (h) ###############
###########################################################################################

###########################################################################################
##### Loading the model results from Type I model #########################################
###########################################################################################
dat_mod_I = read.table(file="https://raw.githubusercontent.com/OscarFHC/IGP_LabExp/master/data/Mod_Type1.csv", 
                       sep=",", header=TRUE, fill=TRUE)
dat_mod_I = dat_mod_I %>%
  mutate(Res = (R1+R2) / (dat_mod_I[1,"R1"] + dat_mod_I[1,"R2"]),
         Cons = C / (dat_mod_I[1,"C"]),
         Pred = P / (dat_mod_I[1,"P"]))

# Type I R2 for IG predator
B_I_R2 <- cbind(as.data.frame(B_cal[1:6,1:6]), subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_I_pred = sum((B_I_R2[,"avg_trmt"] - mean(B_I_R2[,"avg_trmt"]))^2)
SSR_I_pred = sum((B_I_R2[,"avg_trmt"] - B_I_R2[,"Pred"])^2)
R2_I_pred = ifelse((1-(SSR_I_pred/SST_I_pred)) > 0, 1-(SSR_I_pred/SST_I_pred), 0)

# Type I R2 for IG prey
C_I_R2 <- cbind(as.data.frame(C_cal[1:6,1:6]), subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_I_prey = sum((C_I_R2[,"avg_trmt"] - mean(C_I_R2[,"avg_trmt"]))^2)
SSR_I_prey = sum((C_I_R2[,"avg_trmt"] - C_I_R2[,"Cons"])^2)
R2_I_prey = ifelse((1-(SSR_I_prey/SST_I_prey)) > 0, (1-(SSR_I_prey/SST_I_prey)), 0)

# Type I R2 for Bacteria
Bac_I_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(dat_mod_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_I_bac = sum((Bac_I_R2[,"avg_trmt"] - mean(Bac_I_R2[,"avg_trmt"]))^2)
SSR_I_bac = sum((Bac_I_R2[,"avg_trmt"] - Bac_I_R2[,"Res"])^2)
R2_I_bac = ifelse((1-(SSR_I_bac/SST_I_bac)) > 0, (1-(SSR_I_bac/SST_I_bac)), 0)
###########################################################################################
##### Loading the model results from Type I model #########################################
###########################################################################################

###########################################################################################
##### Fig 3. Overlay experimental data with model prediction ##############################
###########################################################################################
Bplot_mod <- ggplot() +
  geom_point(data = B_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = B_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_mod_II, aes(x = alpha, y = Pred), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Pred), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 1.5), expand = c(0, 0.1), breaks = seq(0.4, 2, 0.4)) + 
  labs(x = "",
       y = expression(atop("IG predator density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 0.4, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bplot_mod <- ggdraw(Bplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

Cplot_mod <- ggplot() +
  geom_point(data = C_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = C_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_mod_II, aes(x = alpha, y = Cons), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Cons), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 1.5), expand = c(0, 0.1), breaks = seq(0.4, 2, 0.4)) + 
  labs(x = "",
       y = expression(atop("IG prey density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 0.4, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Cplot_mod <- ggdraw(Cplot_mod) + 
  draw_label("(ratio to 0% availability )", angle = 90, x = 0.085, y = 0.58, size = 14)

Bacplot_mod <- ggplot() +
  geom_point(data = Bac_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_mod_II, aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Res), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 1.5), expand = c(0, 0.1), breaks = seq(0.4, 2, 0.4)) + 
  labs(x = "",
       y = expression(atop("Bacteria density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 0.4, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bacplot_mod <- ggdraw(Bacplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

Fig3 <- plot_grid(Bacplot_mod, Cplot_mod, Bplot_mod,
                  labels = c("a.", "b.", "c."), ncol = 1, align = 'v')

Fig3 <- ggdraw(Fig3) + 
  draw_label("Availability of IG prey (IGP strength)", 
             x = 0.52, y = 0.02, size = 20) +
  draw_label(expression("error bars represent standard error of the mean"), 
             x = 0.86, y = 0.01, size = 12)

ggsave(filename="D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/Fig3_ModOverlay.pdf", 
       plot = Fig3, width = 36, height = 24, units = c("cm"), dpi = 600)
###########################################################################################
##### Fig 3. Overlay experimental data with model prediction ##############################
###########################################################################################




###########################################################################################
##### Fig. S2 model population dynamics ###################################################
###########################################################################################
alpha <- c("a00", "a02", "a04", "a06", "a08", "a10")

TS <- read.table(file = paste0("D:/Research/IGP_LabExp_Public/check/TS_", alpha[1], ".csv"), 
                 sep = ",", header = FALSE, fill = TRUE) %>%
      rename(TimeSteps = V1, R1 = V2, R2 = V3, C = V4, B = V5) %>%
      mutate(alpha = alpha[1])

for (i in 2:length(alpha)){
  TS <- rbind(TS,
              read.table(file = paste0("D:/Research/IGP_LabExp_Public/check/TS_", alpha[i], ".csv"), 
                         sep = ",", header = FALSE, fill = TRUE) %>%
                rename(TimeSteps = V1, R1 = V2, R2 = V3, C = V4, B = V5) %>%
                mutate(alpha = alpha[i])
              )
}
#colnames(TS) <- c("R1", "R2", "C", "B", "TimeSteps", "trmt")
TS <- TS %>% mutate(C = log(C), 
                    B= log(B),
                    int = ifelse(TimeSteps<2500, "First 1000 time steps", "Final 1000 time steps"),
                    int_f = factor(int, levels = c("First 1000 time steps", "Final 1000 time steps")))
head(TS)
str(TS)
B_TS <- ggplot(data = TS, aes(x = TimeSteps, y = B, colour = alpha)) + 
  geom_line(size = 1) + 
  facet_grid(.~int_f, scales = "free_x") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-6, 6.5)) +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
                      guide = guide_legend(byrow=TRUE, title = expression("IGP \nstrength"))) + 
  labs(x = "", y = "", title = "") + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.3, r = 0.4, b = 0.1, l = 1.4, "cm"),
        legend.justification = c("right", "top"),
        legend.margin = margin(t = 0, r = 0.4, b = 0, l = 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(size=14, angle = 18),
        axis.text.y = element_text(size=14),
        panel.spacing = unit(2.1, "lines"),
        strip.text.x = element_text(size = 14))

C_TS <- ggplot(data = TS, aes(x = TimeSteps, y = C + 5.2, colour = alpha)) + 
  geom_line(size = 1) + 
  facet_grid(.~int_f, scales = "free_x") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(-6, 6.5)) +
  scale_color_viridis(discrete = TRUE,
                      labels = c("0%", "20%", "40%", "60%", "80%", "100%"),
                      guide = guide_legend(byrow=TRUE, title = expression("IGP \nstrength"))) + 
  labs(x = "", y = "", title = "") + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.3, r = 0.4, b = 0.1, l = 1.4, "cm"),
        legend.justification = c("right", "top"),
        legend.margin = margin(t = 0, r = 0.4, b = 0, l = 0, "cm"),
        legend.key = element_rect(color = "white", fill = "white"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(size=14, angle = 18),
        axis.text.y = element_text(size=14),
        panel.spacing = unit(2.1, "lines"),
        strip.text.x = element_text(size = 14))

FigS2 <- plot_grid(B_TS + theme(legend.position="none"), 
                   C_TS + theme(legend.position="none"), 
                   labels = c("a.", "b."), ncol = 1, align = 'v')
shared_legend <- get_legend(C_TS + theme(legend.position = c(1.1, 0.58)))
FigS2 <- plot_grid(FigS2, shared_legend, nrow = 1, rel_widths = c(3, 0.25))
FigS2 <- ggdraw(FigS2) + 
  draw_label("Time Steps", 
             x = 0.52, y = 0.02, size = 18) + 
  draw_label("Log(Simulated density)", angle = 90, 
             x = 0.035, y = 0.5, size = 18) + 
  draw_label(expression("IG prey (" * italic("Colpidium") * ")"), 
             x = 0.52, y = 0.5, size = 14) + 
  draw_label(expression("IG predator (" * italic("Blespharisma") * ")"), 
             x = 0.52, y = 0.98, size = 14) 

ggsave(filename = "D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/FigS2_SimTS.jpeg", 
       plot = FigS2, width = 36, height = 24, units = c("cm"), dpi = 600)
###########################################################################################
##### Fig. S2 model population dynamics ###################################################
###########################################################################################

###########################################################################################
##### Fig. S3 empirical Type I ############################################################
###########################################################################################
dat_emp_I = read.table(file = 
                       paste0("https://raw.githubusercontent.com/OscarFHC/IGP_LabExp_Public/master/Empirical_Type1.csv"),
                       sep = ",", header = FALSE, fill = TRUE)
colnames(dat_emp_I) <- c("alpha", "R1", "R2", "Cons", "Pred")
dat_emp_I = dat_emp_I %>%
  mutate(Res = (R1+R2) / (dat_emp_I[1,"R1"] + dat_emp_I[1,"R2"]),
         Cons = Cons / (dat_emp_I[1,"Cons"]),
         Pred = Pred / (dat_emp_I[1,"Pred"])) 

# Type I R2 for IG predator
B_emp_I_R2 <- cbind(as.data.frame(B_cal[,1:6]), subset(dat_emp_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
B_I_Cor <- cor(B_emp_I_R2$avg_trmt, B_emp_I_R2$Pred)
SST_emp_I_pred = sum((B_emp_I_R2[,"avg_trmt"] - mean(B_emp_I_R2[,"avg_trmt"]))^2)
SSR_emp_I_pred = sum((B_emp_I_R2[,"avg_trmt"] - B_emp_I_R2[,"Pred"])^2)
R2_emp_I_pred = ifelse((1-(SSR_emp_I_pred/SST_emp_I_pred)) > 0, 1-(SSR_emp_I_pred/SST_emp_I_pred), 0)

# Type I R2 for IG prey
C_emp_I_R2 <- cbind(as.data.frame(C_cal[,1:6]), subset(dat_emp_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
C_emp_I_Cor <- cor(C_emp_I_R2$avg_trmt, C_emp_I_R2$Cons)
SST_emp_I_prey = sum((C_emp_I_R2[,"avg_trmt"] - mean(C_emp_I_R2[,"avg_trmt"]))^2)
SSR_emp_I_prey = sum((C_emp_I_R2[,"avg_trmt"] - C_emp_I_R2[,"Cons"])^2)
R2_emp_I_prey = ifelse((1-(SSR_emp_I_prey/SST_emp_I_prey)) > 0, (1-(SSR_emp_I_prey/SST_emp_I_prey)), 0)

# Type I R2 for Bacteria
Bac_emp_I_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(dat_emp_I, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
Bac_emp_I_Cor <- cor(Bac_emp_I_R2$avg_trmt, Bac_emp_I_R2$Res)
SST_emp_I_bac = sum((Bac_emp_I_R2[,"avg_trmt"] - mean(Bac_emp_I_R2[,"avg_trmt"]))^2)
SSR_emp_I_bac = sum((Bac_emp_I_R2[,"avg_trmt"] - Bac_emp_I_R2[,"Res"])^2)
R2_emp_I_bac = ifelse((1-(SSR_emp_I_bac/SST_emp_I_bac)) > 0, (1-(SSR_emp_I_bac/SST_emp_I_bac)), 0)

Bplot_mod <- ggplot() +
  geom_point(data = B_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = B_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_emp_I, aes(x = alpha, y = Pred), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Pred), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 1.5), expand = c(0, 0.1), breaks = seq(0.4, 2, 0.4)) + 
  labs(x = "",
       y = expression(atop("IG predator density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 0.4, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bplot_mod <- ggdraw(Bplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

Cplot_mod <- ggplot() +
  geom_point(data = C_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = C_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_emp_I, aes(x = alpha, y = Cons), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Cons), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(-0.1, 2), expand = c(0, 0.1), breaks = seq(0, 2, 0.4)) + 
  labs(x = "",
       y = expression(atop("IG prey density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 0.4, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Cplot_mod <- ggdraw(Cplot_mod) + 
  draw_label("(ratio to 0% availability )", angle = 90, x = 0.085, y = 0.58, size = 14)

Bacplot_mod <- ggplot() +
  geom_point(data = Bac_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_emp_I, aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Res), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 6), expand = c(0, 0.1), breaks = seq(0, 6, 1.5)) + 
  labs(x = "",
       y = expression(atop("Bacteria density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 0.4, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bacplot_mod <- ggdraw(Bacplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

FigS3 <- plot_grid(Bacplot_mod, Cplot_mod, Bplot_mod,
                   labels = c("a.", "b.", "c."), ncol = 1, align = 'v')

FigS3 <- ggdraw(FigS3) + 
  draw_label("Availability of IG prey (IGP strength)", 
             x = 0.52, y = 0.02, size = 20) +
  draw_label(expression("error bars represent standard error of the mean"), 
             x = 0.86, y = 0.01, size = 12)

ggsave(filename="D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/FigS3_Empirical_TypeI.tiff", 
       plot = FigS3, width = 36, height = 24, units = c("cm"), dpi = 600)

###########################################################################################
##### Fig. S3 empirical Type I ############################################################
###########################################################################################

###########################################################################################
##### Fig. S4 Sensitivity on s1 and e1 ####################################################
###########################################################################################
##### e1 15% lower ####################################################
s1_lo15 = read.table(file =  paste0("D:/Research/IGP_LabExp_Public/Type2_s1_lo15.csv"),
                     sep = ",", header = FALSE, fill = TRUE)
colnames(s1_lo15) <- c("alpha", "R1", "R2", "Cons", "Pred")
s1_lo15 = s1_lo15 %>%
  mutate(Res = (R1+R2) / (s1_lo15[1,"R1"] + s1_lo15[1,"R2"]),
         Cons = Cons / (s1_lo15[1,"Cons"]),
         Pred = Pred / (s1_lo15[1,"Pred"]))
B_II_R2 <- cbind(as.data.frame(B_cal[,1:6]), subset(s1_lo15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_pred = sum((B_II_R2[,"avg_trmt"] - mean(B_II_R2[,"avg_trmt"]))^2)
SSR_II_pred = sum((B_II_R2[,"avg_trmt"] - B_II_R2[,"Pred"])^2)
R2_II_pred = ifelse((1-(SSR_II_pred/SST_II_pred)) > 0, 1-(SSR_II_pred/SST_II_pred), 0)

# Type I R2 for IG prey
C_II_R2 <- cbind(as.data.frame(C_cal[,1:6]), subset(s1_lo15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_prey = sum((C_II_R2[,"avg_trmt"] - mean(C_II_R2[,"avg_trmt"]))^2)
SSR_II_prey = sum((C_II_R2[,"avg_trmt"] - C_II_R2[,"Cons"])^2)
R2_II_prey = ifelse((1-(SSR_II_prey/SST_II_prey)) > 0, (1-(SSR_II_prey/SST_II_prey)), 0)

# Type I R2 for Bacteria
Bac_II_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(s1_lo15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_bac = sum((Bac_II_R2[,"avg_trmt"] - mean(Bac_II_R2[,"avg_trmt"]))^2)
SSR_II_bac = sum((Bac_II_R2[,"avg_trmt"] - Bac_II_R2[,"Res"])^2)
R2_II_bac = ifelse((1-(SSR_II_bac/SST_II_bac)) > 0, (1-(SSR_II_bac/SST_II_bac)), 0)

Bplot_mod <- ggplot() +
  geom_point(data = B_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = B_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = s1_lo15, aes(x = alpha, y = Pred), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = e1_lo15, aes(x = alpha, y = Pred), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG predator density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bplot_mod <- ggdraw(Bplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

Cplot_mod <- ggplot() +
  geom_point(data = C_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = C_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = s1_lo15, aes(x = alpha, y = Cons), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Cons), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG prey density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Cplot_mod <- ggdraw(Cplot_mod) + 
  draw_label("(ratio to 0% availability )", angle = 90, x = 0.175, y = 0.58, size = 14)

Bacplot_mod <- ggplot() +
  geom_point(data = Bac_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = s1_lo15, aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Res), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("Bacteria density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bacplot_mod <- ggdraw(Bacplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

FigS4_1 <- plot_grid(Bacplot_mod, Cplot_mod, Bplot_mod,
                     labels = c("", "", ""), ncol = 1, align = 'v')

FigS4_1 <- ggdraw(FigS4_1) + 
  draw_label("Availability of IG prey (IGP strength)", 
             x = 0.6, y = 0.02, size = 20) 
##### s1 15% lower ####################################################

##### s1 15% higher ####################################################
s1_up15 = read.table(file =  paste0("D:/Research/IGP_LabExp_Public/Type2_s1_up15.csv"),
                     sep = ",", header = FALSE, fill = TRUE)
colnames(s1_up15) <- c("alpha", "R1", "R2", "Cons", "Pred")
s1_up15 = s1_up15 %>%
  mutate(Res = (R1+R2) / (s1_up15[1,"R1"] + s1_up15[1,"R2"]),
         Cons = Cons / (s1_up15[1,"Cons"]),
         Pred = Pred / (s1_up15[1,"Pred"]))
B_II_R2 <- cbind(as.data.frame(B_cal[,1:6]), subset(s1_up15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_pred = sum((B_II_R2[,"avg_trmt"] - mean(B_II_R2[,"avg_trmt"]))^2)
SSR_II_pred = sum((B_II_R2[,"avg_trmt"] - B_II_R2[,"Pred"])^2)
R2_II_pred = ifelse((1-(SSR_II_pred/SST_II_pred)) > 0, 1-(SSR_II_pred/SST_II_pred), 0)

# Type I R2 for IG prey
C_II_R2 <- cbind(as.data.frame(C_cal[,1:6]), subset(s1_up15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_prey = sum((C_II_R2[,"avg_trmt"] - mean(C_II_R2[,"avg_trmt"]))^2)
SSR_II_prey = sum((C_II_R2[,"avg_trmt"] - C_II_R2[,"Cons"])^2)
R2_II_prey = ifelse((1-(SSR_II_prey/SST_II_prey)) > 0, (1-(SSR_II_prey/SST_II_prey)), 0)

# Type I R2 for Bacteria
Bac_II_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(s1_up15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_bac = sum((Bac_II_R2[,"avg_trmt"] - mean(Bac_II_R2[,"avg_trmt"]))^2)
SSR_II_bac = sum((Bac_II_R2[,"avg_trmt"] - Bac_II_R2[,"Res"])^2)
R2_II_bac = ifelse((1-(SSR_II_bac/SST_II_bac)) > 0, (1-(SSR_II_bac/SST_II_bac)), 0)

Bplot_mod <- ggplot() +
  geom_point(data = B_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = B_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = s1_up15, aes(x = alpha, y = Pred), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Pred), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG predator density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bplot_mod <- ggdraw(Bplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

Cplot_mod <- ggplot() +
  geom_point(data = C_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = C_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = s1_up15, aes(x = alpha, y = Cons), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Cons), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG prey density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Cplot_mod <- ggdraw(Cplot_mod) + 
  draw_label("(ratio to 0% availability )", angle = 90, x = 0.175, y = 0.58, size = 14)

Bacplot_mod <- ggplot() +
  geom_point(data = Bac_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = s1_up15, aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Res), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("Bacteria density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bacplot_mod <- ggdraw(Bacplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

FigS4_2 <- plot_grid(Bacplot_mod, Cplot_mod, Bplot_mod,
                     labels = c("", "", ""), ncol = 1, align = 'v')

FigS4_2 <- ggdraw(FigS4_2) + 
  draw_label("Availability of IG prey (IGP strength)", 
             x = 0.6, y = 0.02, size = 20) 
##### s1 15% higher ####################################################

##### e1 15% lower ####################################################
e1_lo15 = read.table(file =  paste0("D:/Research/IGP_LabExp_Public/Type2_e1_lo15.csv"),
                     sep = ",", header = FALSE, fill = TRUE)
colnames(e1_lo15) <- c("alpha", "R1", "R2", "Cons", "Pred")
e1_lo15 = e1_lo15 %>%
  mutate(Res = (R1+R2) / (e1_lo15[1,"R1"] + e1_lo15[1,"R2"]),
         Cons = Cons / (e1_lo15[1,"Cons"]),
         Pred = Pred / (e1_lo15[1,"Pred"]))
B_II_R2 <- cbind(as.data.frame(B_cal[,1:6]), subset(e1_lo15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_pred = sum((B_II_R2[,"avg_trmt"] - mean(B_II_R2[,"avg_trmt"]))^2)
SSR_II_pred = sum((B_II_R2[,"avg_trmt"] - B_II_R2[,"Pred"])^2)
R2_II_pred = ifelse((1-(SSR_II_pred/SST_II_pred)) > 0, 1-(SSR_II_pred/SST_II_pred), 0)

# Type I R2 for IG prey
C_II_R2 <- cbind(as.data.frame(C_cal[,1:6]), subset(e1_lo15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_prey = sum((C_II_R2[,"avg_trmt"] - mean(C_II_R2[,"avg_trmt"]))^2)
SSR_II_prey = sum((C_II_R2[,"avg_trmt"] - C_II_R2[,"Cons"])^2)
R2_II_prey = ifelse((1-(SSR_II_prey/SST_II_prey)) > 0, (1-(SSR_II_prey/SST_II_prey)), 0)

# Type I R2 for Bacteria
Bac_II_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(e1_lo15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_bac = sum((Bac_II_R2[,"avg_trmt"] - mean(Bac_II_R2[,"avg_trmt"]))^2)
SSR_II_bac = sum((Bac_II_R2[,"avg_trmt"] - Bac_II_R2[,"Res"])^2)
R2_II_bac = ifelse((1-(SSR_II_bac/SST_II_bac)) > 0, (1-(SSR_II_bac/SST_II_bac)), 0)

Bplot_mod <- ggplot() +
  geom_point(data = B_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = B_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = e1_lo15, aes(x = alpha, y = Pred), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = e1_lo15, aes(x = alpha, y = Pred), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG predator density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bplot_mod <- ggdraw(Bplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

Cplot_mod <- ggplot() +
  geom_point(data = C_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = C_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = e1_lo15, aes(x = alpha, y = Cons), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Cons), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG prey density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Cplot_mod <- ggdraw(Cplot_mod) + 
  draw_label("(ratio to 0% availability )", angle = 90, x = 0.175, y = 0.58, size = 14)

Bacplot_mod <- ggplot() +
  geom_point(data = Bac_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = e1_lo15, aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Res), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("Bacteria density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bacplot_mod <- ggdraw(Bacplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

FigS4_3 <- plot_grid(Bacplot_mod, Cplot_mod, Bplot_mod,
                     labels = c("", "", ""), ncol = 1, align = 'v')

FigS4_3 <- ggdraw(FigS4_3) + 
  draw_label("Availability of IG prey (IGP strength)", 
             x = 0.6, y = 0.02, size = 20) 
##### e1 15% lower ####################################################

##### e1 15% higher ####################################################
e1_up15 = read.table(file =  paste0("D:/Research/IGP_LabExp_Public/Type2_e1_up15.csv"),
                     sep = ",", header = FALSE, fill = TRUE)
colnames(e1_up15) <- c("alpha", "R1", "R2", "Cons", "Pred")
e1_up15 = e1_up15 %>%
  mutate(Res = (R1+R2) / (e1_up15[1,"R1"] + e1_up15[1,"R2"]),
         Cons = Cons / (e1_up15[1,"Cons"]),
         Pred = Pred / (e1_up15[1,"Pred"]))
B_II_R2 <- cbind(as.data.frame(B_cal[,1:6]), subset(e1_up15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_pred = sum((B_II_R2[,"avg_trmt"] - mean(B_II_R2[,"avg_trmt"]))^2)
SSR_II_pred = sum((B_II_R2[,"avg_trmt"] - B_II_R2[,"Pred"])^2)
R2_II_pred = ifelse((1-(SSR_II_pred/SST_II_pred)) > 0, 1-(SSR_II_pred/SST_II_pred), 0)

# Type I R2 for IG prey
C_II_R2 <- cbind(as.data.frame(C_cal[,1:6]), subset(e1_up15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_prey = sum((C_II_R2[,"avg_trmt"] - mean(C_II_R2[,"avg_trmt"]))^2)
SSR_II_prey = sum((C_II_R2[,"avg_trmt"] - C_II_R2[,"Cons"])^2)
R2_II_prey = ifelse((1-(SSR_II_prey/SST_II_prey)) > 0, (1-(SSR_II_prey/SST_II_prey)), 0)

# Type I R2 for Bacteria
Bac_II_R2 <- cbind(as.data.frame(Bac_cal[1:6,1:6]), subset(e1_up15, alpha %in% c(0, 0.2, 0.4, 0.6, 0.8, 1.0)))
SST_II_bac = sum((Bac_II_R2[,"avg_trmt"] - mean(Bac_II_R2[,"avg_trmt"]))^2)
SSR_II_bac = sum((Bac_II_R2[,"avg_trmt"] - Bac_II_R2[,"Res"])^2)
R2_II_bac = ifelse((1-(SSR_II_bac/SST_II_bac)) > 0, (1-(SSR_II_bac/SST_II_bac)), 0)

Bplot_mod <- ggplot() +
  geom_point(data = B_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = B_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = e1_up15, aes(x = alpha, y = Pred), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Pred), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG predator density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bplot_mod <- ggdraw(Bplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

Cplot_mod <- ggplot() +
  geom_point(data = C_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = C_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = e1_up15, aes(x = alpha, y = Cons), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Cons), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("IG prey density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Cplot_mod <- ggdraw(Cplot_mod) + 
  draw_label("(ratio to 0% availability )", angle = 90, x = 0.175, y = 0.58, size = 14)

Bacplot_mod <- ggplot() +
  geom_point(data = Bac_II_R2, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = e1_up15, aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = e1_up15, aes(x = alpha, y = Res), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0.1), breaks = seq(0, 2.5, 0.5)) + 
  labs(x = "",
       y = expression(atop("Bacteria density", "at the steady state")))+ 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(t = 0.5, r = 1.0, b = 0.6, l = 1.0, "cm"),
        axis.line = element_line(colour = "black"), 
        axis.text = element_text(size = 20), 
        axis.title.y = element_text(size = 20, margin = margin(t = 0, r = 22, b = 0, l = 0)))
Bacplot_mod <- ggdraw(Bacplot_mod) + 
  draw_label("", angle = 90, x = 0.085, y = 0.55, size = 14)

FigS4_4 <- plot_grid(Bacplot_mod, Cplot_mod, Bplot_mod,
                     labels = c("", "", ""), ncol = 1, align = 'v')

FigS4_4 <- ggdraw(FigS4_4) + 
  draw_label("Availability of IG prey (IGP strength)", 
             x = 0.6, y = 0.02, size = 20) 
##### e1 15% higher ####################################################

FigS4 <- plot_grid(FigS4_1, FigS4_2, FigS4_3, FigS4_4,
                   labels = c("a.", "b.", "c.", "d."), ncol = 2, align = 'v')
# ggsave(filename="D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/FigS4_Sensitivity.eps", 
#        plot = FigS4, width = 36, height = 42, units = c("cm"), dpi = 600)
###########################################################################################
##### Fig. S4 Sensitivity on s1 and e1 ####################################################
###########################################################################################

##### Calculating the Blephasrima's consumption rate differences between 100% vs 20% IGP treatments
Bac02 <- mean(c(Bac_pop[which(Bac_pop$trmt_all == "Bac02"),]$Tend, Bac_pop[which(Bac_pop$trmt_all == "Bac10"),]$Tend))/100
Bac00 <- mean(c(Bac_pop[which(Bac_pop$trmt_all == "Bac02"),]$Tend, Bac_pop[which(Bac_pop$trmt_all == "Bac10"),]$Tend))/100
Ble00 <- B_sum_00$avg#*1.085409
Ble02 <- B_sum_00$avg#*1.071174
Col00 <- C_sum_00$avg#*0.616693
Col02 <- mean(C_sum_00$avg)*1.074068

all <- Ble00 * ( (1.25 * Bac00)/(1 + 1.25 * 0.8 * Bac00 + 0.3858 * 0.35959 * Col00) )

all02 <- Ble02 * ( (1.25 * Bac02*0.2)/(1 + 1.25 * 0.8 * Bac02*0.2 + 0.3858 * 0.35959 * Col00*0.2) )

t.test(all, all02)
mean(all)
sd(all)/sqrt(5)

mean(all02)
sd(all02)/sqrt(5)
