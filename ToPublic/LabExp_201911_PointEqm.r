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

##### Loading other parameter value estimation ############################################

### refit 
mod.lst <- expand.grid(e2 = seq(0, 1, 0.1), partition = seq(0, 1, 0.1))

Bac_pop_mod <- Bac_pop[-grep("ctrl", Bac_pop[,"trmt_all"]),] %>%
  mutate(trmt = as.numeric(c(rep(c("0", "0.2", "0.4", "0.6"), each=3), rep(("0.8"), each=4), rep(c("1.0"), each=5))),
         trmt2 = trmt^2)

for (i in 1:nrow(mod.lst)){
  temp <- read.table(file = paste0("D:/Research/IGP_LabExp/Analysis/ParamRefit/r015/Type2_e2_", 
                                   mod.lst[i, "e2"]*10, "._p_", mod.lst[i,"partition"]*100, "..csv"), sep=",", header=FALSE)
  colnames(temp) <- c("alpha", "R1", "R2", "Cons", "Pred")
  temp = temp %>%
    mutate(Res = (R1+R2) / (temp[1,"R1"] + temp[1,"R2"]),
           Cons = Cons / (temp[1,"Cons"]),
           Pred = Pred / (temp[1,"Pred"])) 
  
  mod.lst[i, "RSS_Pred"] <- sum((B_cal[["avg_trmt"]] - temp[,"Pred"])^2)/
    sum((B_cal[["avg_trmt"]] - mean(temp[,"Pred"]))^2)
  mod.lst[i, "RSS_Prey"] <- sum((C_cal[["avg_trmt"]] - temp[,"Cons"])^2)/
    sum((C_cal[["avg_trmt"]] - mean(temp[,"Cons"]))^2)
  
  tempdat_Bac <- data.frame(alpha = c(rep(c(0, 0.2, 0.4, 0.6), each = 3), rep(0.8, 4), rep(1.0, 5))) %>%
    mutate(Res = c(rep(temp[c(1:4), "Res"], 3), rep(temp[5, "Res"], 4), rep(temp[6, "Res"], 5)))
  mod.lst[i, "RSS_Bac"] <- sum((Bac_pop_mod[["Change_Eqm"]] - tempdat_Bac[,"Res"])^2)/
    sum((Bac_pop_mod[["Change_Eqm"]] - mean(Bac_pop_mod[["Change_Eqm"]]))^2)
  }

mod.lst <- mod.lst %>%
  mutate(RSS_tot = RSS_Bac)# + RSS_Prey + RSS_Pred) 

e2 <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["e2"]]
partition <- subset(mod.lst, RSS_tot == min(mod.lst[,"RSS_tot"]))[["partition"]]


##### Loading other parameter value estimation ############################################


###########################################################################################
##### Revised Type II #####################################################################
###########################################################################################
# Revised Type II
dat_mod_II_R = read.table(file = "D:/Research/IGP_LabExp_Public/Type2ModResults_R/Type2_e1040_r015.csv",
                          sep = ",", header = FALSE, fill = TRUE)
colnames(dat_mod_II_R) <- c("alpha", "R1", "R2", "Cons", "Pred")
dat_mod_II_R = dat_mod_II_R %>%
  mutate(Res = (R1+R2) / (dat_mod_II_R[1,"R1"] + dat_mod_II_R[1,"R2"]),
         Cons = Cons / (dat_mod_II_R[1,"Cons"]),
         Pred = Pred / (dat_mod_II_R[1,"Pred"])) 

B_II_R2_R <- cbind(as.data.frame(B_cal[,1:6]), dat_mod_II_R)
SST_II_pred_R = sum((B_II_R2_R[,"avg_trmt"] - mean(B_II_R2_R[,"avg_trmt"]))^2)
SSR_II_pred_R = sum((B_II_R2_R[,"avg_trmt"] - B_II_R2_R[,"Pred"])^2)
R2_II_pred_R = ifelse((1-(SSR_II_pred_R/SST_II_pred_R)) > 0, 1-(SSR_II_pred_R/SST_II_pred_R), 0)

# Type I R2 for IG prey
C_II_R2_R <- cbind(as.data.frame(C_cal[,1:6]), dat_mod_II_R)
SST_II_prey_R = sum((C_II_R2_R[,"avg_trmt"] - mean(C_II_R2_R[,"avg_trmt"]))^2)
SSR_II_prey_R = sum((C_II_R2_R[,"avg_trmt"] - C_II_R2_R[,"Cons"])^2)
R2_II_prey_R = ifelse((1-(SSR_II_prey_R/SST_II_prey_R)) > 0, (1-(SSR_II_prey_R/SST_II_prey_R)), 0)

# Type I R2 for Bacteria
Bac_II_R2_R <- cbind(as.data.frame(Bac_cal[1:6,1:6]), dat_mod_II_R)
SST_II_bac_R = sum((Bac_II_R2_R[,"avg_trmt"] - mean(Bac_II_R2_R[,"avg_trmt"]))^2)
SSR_II_bac_R = sum((Bac_II_R2_R[,"avg_trmt"] - Bac_II_R2_R[,"Res"])^2)
R2_II_bac_R = ifelse((1-(SSR_II_bac_R/SST_II_bac_R)) > 0, (1-(SSR_II_bac_R/SST_II_bac_R)), 0)

###########################################################################################
##### Revised Type II #####################################################################
###########################################################################################

###########################################################################################
##### Fig. S2 model population dynamics ###################################################
###########################################################################################
alpha <- c("a00", "a02", "a04", "a06", "a08", "a10")

TS <- read.table(file = paste0("D:/Research/IGP_LabExp_Public/TScheck/TS_e1040_r015_", alpha[1], ".csv"), 
                 sep = ",", header = FALSE, fill = TRUE) %>%
  rename(TimeSteps = V1, R1 = V2, R2 = V3, C = V4, B = V5) %>%
  mutate(alpha = alpha[1])

for (i in 2:length(alpha)){
  TS <- rbind(TS,
              read.table(file = paste0("D:/Research/IGP_LabExp_Public/TScheck/TS_e1040_r015_", alpha[i], ".csv"), 
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
# head(TS)
# str(TS)
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

ggsave(filename = "D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/ParamRefit/SimTS_e1040_r015.jpeg", 
       plot = FigS2, width = 36, height = 24, units = c("cm"), dpi = 600)
###########################################################################################
##### Fig. S2 model population dynamics ###################################################
###########################################################################################


###########################################################################################
##### Overlay experimental data with model prediction #####################################
###########################################################################################
Bplot_mod <- ggplot() +
  geom_point(data = B_II_R2_R, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = B_II_R2_R, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_mod_II_R, aes(x = alpha, y = Pred), method = 'loess', se = FALSE, size = 1, color = "black") +
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 2), expand = c(0, 0.1), breaks = seq(0.4, 2, 0.4)) + 
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
  geom_point(data = C_II_R2_R, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = C_II_R2_R, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_mod_II_R, aes(x = alpha, y = Cons), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Cons), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 2), expand = c(0, 0.1), breaks = seq(0.4, 2, 0.4)) + 
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
  geom_point(data = Bac_II_R2_R, aes(x = alpha, y = avg_trmt), size = 2) +  
  geom_errorbar(data = Bac_II_R2_R, aes(x = alpha, ymin = avg_trmt - se_permu, ymax = avg_trmt + se_permu), width = 0.02) + 
  geom_smooth(data = dat_mod_II_R, aes(x = alpha, y = Res), method = 'loess', se = FALSE, size = 1, color = "black") +
  #geom_line(data = dat_mod_I, aes(x = alpha, y = Res), size = 1, linetype = 2) + 
  scale_x_continuous(limits = c(-0.02, 1.02), expand = c(0, 0), breaks = seq(0, 1, 0.2),
                     labels = c("0%", "20%", "40%", "60%", "80%", "100%")) + 
  scale_y_continuous(limits = c(0.4, 4), expand = c(0, 0.1), breaks = seq(0.4, 4, 1.2)) + 
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

ggsave(filename="D:/Manuscript/IGP_DivEffects_MS_Figs/LabExp_MS_Figs/Ver2/ParamRefit/Overlay_e1040_r015.jpeg", 
       plot = Fig3, width = 36, height = 24, units = c("cm"), dpi = 600)
###########################################################################################
##### Fig 3. Overlay experimental data with model prediction ##############################
###########################################################################################


