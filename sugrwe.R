################################################################

# Set working directory
setwd("D:/OneDrive/OneDrive - Aarhus universitet/Desktop/Experiment/CEH/stat") # Desktop at home
setwd("D:/AU/OneDrive - Aarhus Universitet/Desktop/Experiment/CEH/stat")# Huawei Laptop
setwd("C:/Users/au700351/OneDrive - Aarhus Universitet/Desktop/Experiment/CEH/stat") # AU laptop

setwd("E:/AU/Desktop/Experiment/CEH/stat") # AU laptop


# Load library
library(readr)
library(tidyverse)
library(dplyr)
library(plotrix)
library(ggplot2)
library(scales)
library(rgl)
library(patchwork)
require(dplyr)
library(extrafont)
library(emmeans)
library(robustbase)


font_import(paths = "font/")
# 
cu.list <- c(4.2,20.9,40.4)
swp.list <- c(0.5, 0.7, 1.0, 1.6, 1.8, 2.0)
temp.list <- c(4.5, 9.7, 20.1, 24.1, 27.7)

# Input origin data
dat.org <- read_csv('sugrwe.csv') %>% as_tibble() %>% mutate(., vict = number-survivor)

dat.desc <- dat.org %>% 
  group_by(sqrcu, temp, lgswp) %>% 
  summarise(., survave = mean(survival, na.rm = T),
            survse = std.error(survival, na.rm = T),
            grwave = mean(growth, na.rm = T),
            grwse = std.error(growth, na.rm = T),
            repave = mean(repro, na.rm = T),
            repse = std.error(repro, na.rm = T))


# Build model for visualisation
draw.line <- function() {
lgswp.seq <- seq(swp.list[1], swp.list[6], 0.01)
line.p <- vector("list", length = length(temp.list)*length(cu.list))
for (t in 1:length(temp.list)) {
  if (t == 1) i <- 1
  for (cu in 1:length(cu.list)) {
    dat.p1 <- subset(dat.org, sqrcu == cu.list[cu] & temp == temp.list[t])
    sur.mdl <- glm(cbind(survivor, vict) ~ lgswp, data = dat.p1, family = binomial(link = "logit"))
    grw.mdl <- glm(growth ~ lgswp, data = dat.p1)
    repr.mdl <- glm(repro ~ lgswp, data = dat.p1, family = poisson(link = "log"))
    
    line.p[[i]] <- tibble(lgswp = lgswp.seq, temp = temp.list[t], sqrcu = cu.list[cu], 
                      surv = predict(sur.mdl, newdata = data.frame(lgswp = lgswp.seq), type = "response"),
                      growth = predict(grw.mdl, newdata = data.frame(lgswp = lgswp.seq), type = "response"),
                      repro = predict(repr.mdl, newdata = data.frame(lgswp = lgswp.seq), type = "response"))
    i <- i+1
  }
}

line.t <- do.call(bind_rows, line.p)
return(line.t)
}

line.t <- draw.line()

# Plot by group ------
p.sur <- ggplot()+
  geom_point(aes(x = lgswp, y = survave*1, group = interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, size = 3, position = position_dodge(width = 0.1))+
  geom_errorbar(aes(x = lgswp, ymin = (survave-survse)*1, ymax = (survave+survse)*1, group =  interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, linewidth = 1, width = 0.2, position = position_dodge(width = 0.1))+
  geom_line(aes(x = lgswp, y = surv, color = factor(sqrcu)), data = line.t )+
  scale_color_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  scale_fill_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  # scale_y_continuous(transform = 'log', breaks = c(1, 10, 100, 500, 1000, 3000))+
  labs(x = 'Soil Water potential /-kPa', y = 'Survival /%', color = 'Cu level')+
  # ylim(c(0,200))+
  facet_wrap(~temp, nrow = 1)+
  theme_bw(base_family = "Aptos", base_size = 18)

# Plot by group ------
p.grw <- ggplot()+
  geom_point(aes(x = lgswp, y = grwave*1, group = interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, size = 3, position = position_dodge(width = 0.1))+
  geom_errorbar(aes(x = lgswp, ymin = (grwave-grwse)*1, ymax = (grwave+grwse)*1, group =  interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, linewidth = 1, width = 0.2, position = position_dodge(width = 0.1))+
  geom_line(aes(x = lgswp, y = growth, color = factor(sqrcu)), data = line.t )+
  scale_color_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  scale_fill_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  # scale_y_continuous(transform = 'log', breaks = c(1, 10, 100, 500, 1000, 3000))+
  labs(x = 'Soil Water potential /-kPa', y = 'Average Weight /%', color = 'Cu level')+
  # ylim(c(0,200))+
  facet_wrap(~temp, nrow = 1)+
  theme_bw(base_family = "Aptos", base_size = 18)

p.repr <- ggplot()+
  geom_point(aes(x = lgswp, y = repave*1, group = interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, size = 3, position = position_dodge(width = 0.1))+
  geom_errorbar(aes(x = lgswp, ymin = (repave-repse)*1, ymax = (repave+repse)*1, group =  interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, linewidth = 1, width = 0.2, position = position_dodge(width = 0.1))+
  geom_line(aes(x = lgswp, y = repro, color = factor(sqrcu)), data = line.t )+
  scale_color_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  scale_fill_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  # scale_y_continuous(transform = 'log', breaks = c(1, 10, 100, 500, 1000, 3000))+
  labs(x = 'Soil Water potential /-kPa', y = 'Num. of Offsprings', color = 'Cu level')+
  # ylim(c(0,200))+
  facet_wrap(~temp, nrow = 1)+
  theme_bw(base_family = "Aptos", base_size = 18)


# Statistics

# Surival affected by moisture and Cu effect under selected temperature
sur.under.temp.f <- function() {
for (t in 1:length(temp.list)) {
  if (t == 1) {i <- 1;stat.temp <- vector("list", length = length(temp.list))}
    lgswp.mdl <- glm(cbind(survivor, vict) ~ lgswp, data = subset(dat.org, temp == temp.list[t]), family = binomial(link = "logit"))
    sqrcu.mdl <- glm(cbind(survivor, vict) ~ sqrcu, data = subset(dat.org, temp == temp.list[t]), family = binomial(link = "logit"))
    int.mdl <- glm(cbind(survivor, vict) ~ lgswp*sqrcu, data = subset(dat.org, temp == temp.list[t]), family = binomial(link = "logit"))

    stat.temp[[i]] <- tibble(temp = temp.list[t], 
                          swp.coef = round(summary(lgswp.mdl)$coefficients[2,1], digits = 4),
                          swp.lw = round(confint(lgswp.mdl)[2,1], digits = 4),
                          swp.up = round(confint(lgswp.mdl)[2,2], digits = 4),
                          swp.p = round(summary(lgswp.mdl)$coefficients[2,4], digits = 4),
                          cu.coef = round(summary(sqrcu.mdl)$coefficients[2,1], digits = 4),
                          cu.lw = round(confint(sqrcu.mdl)[2,1], digits = 4),
                          cu.up = round(confint(sqrcu.mdl)[2,2], digits = 4),
                          cu.p = round(summary(sqrcu.mdl)$coefficients[2,4], digits = 4),
                          int.coef = round(summary(int.mdl)$coefficients[4,1], digits = 4),
                          int.lw = round(confint(int.mdl)[4,1], digits = 4),
                          int.up = round(confint(int.mdl)[4,2], digits = 4),
                          int.p = round(summary(int.mdl)$coefficients[4,4], digits = 4)) %>% 
      mutate(., swp.ps = case_when(
        swp.p < 0.0001 ~ '****',
        between(swp.p, 0.0001, 0.001) ~ '***',
        between(swp.p, 0.001, 0.01) ~ '**',
        between(swp.p, 0.01, 0.05) ~ '*',
        swp.p > 0.05 ~ 'n.s.'),
        cu.ps = case_when(
          cu.p < 0.0001 ~ '****',
          between(cu.p, 0.0001, 0.001) ~ '***',
          between(cu.p, 0.001, 0.01) ~ '**',
          between(cu.p, 0.01, 0.05) ~ '*',
          cu.p > 0.05 ~ 'n.s.'),
        int.ps = case_when(
          int.p < 0.0001 ~ '****',
          between(int.p, 0.0001, 0.001) ~ '***',
          between(int.p, 0.001, 0.01) ~ '**',
          between(int.p, 0.01, 0.05) ~ '*',
          int.p > 0.05 ~ 'n.s.')) %>% 
      relocate(swp.ps, .after = swp.p) %>% 
      relocate(cu.ps, .after = cu.p)
    
    i <- i+1
    if (t == length(temp.list)) stat.temp.table <- do.call(bind_rows, stat.temp)
}
 return(stat.temp.table)
}

sur.under.temp <- sur.under.temp.f()
write.csv(sur.under.temp, file = 'stat_sur_table1.csv')

# Survival affected by temp and lgswp under selected Cu levels
sur.under.Cu.f <- function(){
for (t in 1:length(cu.list)) {
  if (t == 1) {i <- 1; stat.cu1 <- vector("list", length = length(cu.list)); stat.cu2 <- vector("list", length = length(cu.list))}
  
  lgswp.mdl1 <- glm(cbind(survivor, vict) ~ lgswp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22), family = binomial(link = "logit"))
  temp.mdl1 <- glm(cbind(survivor, vict) ~ temp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22), family = binomial(link = "logit"))
  int.mdl1 <- glm(cbind(survivor, vict) ~ lgswp*temp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22), family = binomial(link = "logit"))
  
  lgswp.mdl2 <- glm(cbind(survivor, vict) ~ lgswp, data = subset(dat.org,  sqrcu == cu.list[t] & temp > 19), family = binomial(link = "logit"))
  temp.mdl2 <- glm(cbind(survivor, vict) ~ temp, data = subset(dat.org, sqrcu == cu.list[t] & temp > 19), family = binomial(link = "logit"))
  int.mdl2 <- glm(cbind(survivor, vict) ~ lgswp*temp, data = subset(dat.org, sqrcu == cu.list[t] & temp > 19), family = binomial(link = "logit"))
  
  stat.cu1[[i]] <- tibble(temp = "<20",
                          sqrcu = cu.list[t], 
                           swp.coef = round(summary(lgswp.mdl1)$coefficients[2,1], digits = 2),
                           swp.lw = round(confint(lgswp.mdl1)[2,1], digits = 2),
                           swp.up = round(confint(lgswp.mdl1)[2,2], digits = 2),
                           swp.p = round(summary(lgswp.mdl1)$coefficients[2,4], digits = 4),
                           temp.coef = round(summary(temp.mdl1)$coefficients[2,1], digits = 2),
                           temp.lw = round(confint(temp.mdl1)[2,1], digits = 2),
                           temp.up = round(confint(temp.mdl1)[2,2], digits = 2),
                           temp.p = round(summary(temp.mdl1)$coefficients[2,4], digits = 4),
                           int.coef = round(summary(int.mdl1)$coefficients[4,1], digits = 2),
                           int.lw = round(confint(int.mdl1)[4,1], digits = 2),
                           int.up = round(confint(int.mdl1)[4,2], digits = 2),
                           int.p = round(summary(int.mdl1)$coefficients[4,4], digits = 4)) %>% 
    mutate(., swp.ps = case_when(
      swp.p < 0.0001 ~ '****',
      between(swp.p, 0.0001, 0.001) ~ '***',
      between(swp.p, 0.001, 0.01) ~ '**',
      between(swp.p, 0.01, 0.05) ~ '*',
      swp.p > 0.05 ~ 'n.s.'),
      temp.ps = case_when(
        temp.p < 0.0001 ~ '****',
        between(temp.p, 0.0001, 0.001) ~ '***',
        between(temp.p, 0.001, 0.01) ~ '**',
        between(temp.p, 0.01, 0.05) ~ '*',
        temp.p > 0.05 ~ 'n.s.'),
      int.ps = case_when(
        int.p < 0.0001 ~ '****',
        between(int.p, 0.0001, 0.001) ~ '***',
        between(int.p, 0.001, 0.01) ~ '**',
        between(int.p, 0.01, 0.05) ~ '*',
        int.p > 0.05 ~ 'n.s.')) %>% 
    relocate(swp.ps, .after = swp.p) %>% 
    relocate(temp.ps, .after = temp.p)
  
  stat.cu2[[i]] <- tibble(temp = ">20",
                          sqrcu = cu.list[t],
                          swp.coef = round(summary(lgswp.mdl2)$coefficients[2,1], digits = 2),
                          swp.lw = round(confint(lgswp.mdl2)[2,1], digits = 2),
                          swp.up = round(confint(lgswp.mdl2)[2,2], digits = 2),
                          swp.p = round(summary(lgswp.mdl2)$coefficients[2,4], digits = 4),
                          temp.coef = round(summary(temp.mdl2)$coefficients[2,1], digits = 2),
                          temp.lw = round(confint(temp.mdl2)[2,1], digits = 2),
                          temp.up = round(confint(temp.mdl2)[2,2], digits = 2),
                          temp.p = round(summary(temp.mdl2)$coefficients[2,4], digits = 4),
                          int.coef = round(summary(int.mdl2)$coefficients[4,1], digits = 2),
                          int.lw = round(confint(int.mdl2)[4,1], digits = 2),
                          int.up = round(confint(int.mdl2)[4,2], digits = 2),
                          int.p = round(summary(int.mdl2)$coefficients[4,4], digits = 4)) %>% 
    mutate(., swp.ps = case_when(
      swp.p < 0.0001 ~ '****',
      between(swp.p, 0.0001, 0.001) ~ '***',
      between(swp.p, 0.001, 0.01) ~ '**',
      between(swp.p, 0.01, 0.05) ~ '*',
      swp.p > 0.05 ~ 'n.s.'),
      temp.ps = case_when(
        temp.p < 0.0001 ~ '****',
        between(temp.p, 0.0001, 0.001) ~ '***',
        between(temp.p, 0.001, 0.01) ~ '**',
        between(temp.p, 0.01, 0.05) ~ '*',
        temp.p > 0.05 ~ 'n.s.'),
      int.ps = case_when(
        int.p < 0.0001 ~ '****',
        between(int.p, 0.0001, 0.001) ~ '***',
        between(int.p, 0.001, 0.01) ~ '**',
        between(int.p, 0.01, 0.05) ~ '*',
        int.p > 0.05 ~ 'n.s.')) %>% 
    relocate(swp.ps, .after = swp.p) %>% 
    relocate(temp.ps, .after = temp.p)
  
  i <- i+1
  if (t == length(cu.list)) {
    stat.cu.table1 <- do.call(bind_rows, stat.cu1)
    stat.cu.table2 <- do.call(bind_rows, stat.cu2)
    stat.cu.table <- bind_rows(stat.cu.table1, stat.cu.table2)}
}
  return(stat.cu.table)
}

sur.under.cu <- sur.under.Cu.f()
write.csv(sur.under.cu, file = 'stat_sur_table2.csv')

# Growth affected by moisture and Cu effect under selected temperature
grw.under.temp.f <- function() {
  for (t in 1:length(temp.list)) {
    if (t == 1) {i <- 1;stat.temp <- vector("list", length = length(temp.list))}
    lgswp.mdl <- glm(growth ~ lgswp, data = subset(dat.org, temp == temp.list[t]))
    sqrcu.mdl <- glm(growth ~ sqrcu, data = subset(dat.org, temp == temp.list[t]))
    int.mdl <- glm(growth ~ lgswp*sqrcu, data = subset(dat.org, temp == temp.list[t]))
    
    stat.temp[[i]] <- tibble(temp = temp.list[t], 
                             swp.coef = round(summary(lgswp.mdl)$coefficients[2,1], digits = 6),
                             swp.lw = round(confint(lgswp.mdl)[2,1], digits = 6),
                             swp.up = round(confint(lgswp.mdl)[2,2], digits = 6),
                             swp.p = round(summary(lgswp.mdl)$coefficients[2,4], digits = 4),
                             cu.coef = round(summary(sqrcu.mdl)$coefficients[2,1], digits = 6),
                             cu.lw = round(confint(sqrcu.mdl)[2,1], digits = 6),
                             cu.up = round(confint(sqrcu.mdl)[2,2], digits = 6),
                             cu.p = round(summary(sqrcu.mdl)$coefficients[2,4], digits = 4),
                             int.coef = round(summary(int.mdl)$coefficients[4,1], digits = 6),
                             int.lw = round(confint(int.mdl)[4,1], digits = 6),
                             int.up = round(confint(int.mdl)[4,2], digits = 6),
                             int.p = round(summary(int.mdl)$coefficients[4,4], digits = 4)) %>% 
      mutate(., swp.ps = case_when(
        swp.p < 0.0001 ~ '****',
        between(swp.p, 0.0001, 0.001) ~ '***',
        between(swp.p, 0.001, 0.01) ~ '**',
        between(swp.p, 0.01, 0.05) ~ '*',
        swp.p > 0.05 ~ 'n.s.'),
        cu.ps = case_when(
          cu.p < 0.0001 ~ '****',
          between(cu.p, 0.0001, 0.001) ~ '***',
          between(cu.p, 0.001, 0.01) ~ '**',
          between(cu.p, 0.01, 0.05) ~ '*',
          cu.p > 0.05 ~ 'n.s.'),
        int.ps = case_when(
          int.p < 0.0001 ~ '****',
          between(int.p, 0.0001, 0.001) ~ '***',
          between(int.p, 0.001, 0.01) ~ '**',
          between(int.p, 0.01, 0.05) ~ '*',
          int.p > 0.05 ~ 'n.s.')) %>% 
      relocate(swp.ps, .after = swp.p) %>% 
      relocate(cu.ps, .after = cu.p)
    
    i <- i+1
    if (t == length(temp.list)) stat.temp.table <- do.call(bind_rows, stat.temp)
  }
  return(stat.temp.table)
}

grw.under.temp <- grw.under.temp.f()
write.csv(grw.under.temp, file = 'stat_grw_table1.csv')

# Growth affected by temp and lgswp under selected Cu levels
grw.under.Cu.f <- function(){
  for (t in 1:length(cu.list)) {
    if (t == 1) {i <- 1; stat.cu1 <- vector("list", length = length(cu.list)); stat.cu2 <- vector("list", length = length(cu.list))}
    
    lgswp.mdl1 <- glm(growth ~ lgswp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22))
    temp.mdl1 <- glm(growth ~ temp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22))
    int.mdl1 <- glm(growth ~ lgswp*temp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22))
    
    lgswp.mdl2 <- glm(growth ~ lgswp, data = subset(dat.org,  sqrcu == cu.list[t] & temp > 19))
    temp.mdl2 <- glm(growth ~ temp, data = subset(dat.org, sqrcu == cu.list[t] & temp > 19))
    int.mdl2 <- glm(growth ~ lgswp*temp, data = subset(dat.org, sqrcu == cu.list[t] & temp > 19))
    
    stat.cu1[[i]] <- tibble(temp = "<20",
                            sqrcu = cu.list[t], 
                            swp.coef = round(summary(lgswp.mdl1)$coefficients[2,1], digits = 4),
                            swp.lw = round(confint(lgswp.mdl1)[2,1], digits = 4),
                            swp.up = round(confint(lgswp.mdl1)[2,2], digits = 4),
                            swp.p = round(summary(lgswp.mdl1)$coefficients[2,4], digits = 4),
                            temp.coef = round(summary(temp.mdl1)$coefficients[2,1], digits = 4),
                            temp.lw = round(confint(temp.mdl1)[2,1], digits = 4),
                            temp.up = round(confint(temp.mdl1)[2,2], digits = 4),
                            temp.p = round(summary(temp.mdl1)$coefficients[2,4], digits = 4),
                            int.coef = round(summary(int.mdl1)$coefficients[4,1], digits = 6),
                            int.lw = round(confint(int.mdl1)[4,1], digits = 6),
                            int.up = round(confint(int.mdl1)[4,2], digits = 6),
                            int.p = round(summary(int.mdl1)$coefficients[4,4], digits = 4)) %>% 
      mutate(., swp.ps = case_when(
        swp.p < 0.0001 ~ '****',
        between(swp.p, 0.0001, 0.001) ~ '***',
        between(swp.p, 0.001, 0.01) ~ '**',
        between(swp.p, 0.01, 0.05) ~ '*',
        swp.p > 0.05 ~ 'n.s.'),
        temp.ps = case_when(
          temp.p < 0.0001 ~ '****',
          between(temp.p, 0.0001, 0.001) ~ '***',
          between(temp.p, 0.001, 0.01) ~ '**',
          between(temp.p, 0.01, 0.05) ~ '*',
          temp.p > 0.05 ~ 'n.s.'),
        int.ps = case_when(
          int.p < 0.0001 ~ '****',
          between(int.p, 0.0001, 0.001) ~ '***',
          between(int.p, 0.001, 0.01) ~ '**',
          between(int.p, 0.01, 0.05) ~ '*',
          int.p > 0.05 ~ 'n.s.')) %>% 
      relocate(swp.ps, .after = swp.p) %>% 
      relocate(temp.ps, .after = temp.p)
    
    stat.cu2[[i]] <- tibble(temp = ">20",
                            sqrcu = cu.list[t],
                            swp.coef = round(summary(lgswp.mdl2)$coefficients[2,1], digits = 4),
                            swp.lw = round(confint(lgswp.mdl2)[2,1], digits = 4),
                            swp.up = round(confint(lgswp.mdl2)[2,2], digits = 4),
                            swp.p = round(summary(lgswp.mdl2)$coefficients[2,4], digits = 4),
                            temp.coef = round(summary(temp.mdl2)$coefficients[2,1], digits = 4),
                            temp.lw = round(confint(temp.mdl2)[2,1], digits = 4),
                            temp.up = round(confint(temp.mdl2)[2,2], digits = 4),
                            temp.p = round(summary(temp.mdl2)$coefficients[2,4], digits = 4),
                            int.coef = round(summary(int.mdl2)$coefficients[4,1], digits = 6),
                            int.lw = round(confint(int.mdl2)[4,1], digits = 6),
                            int.up = round(confint(int.mdl2)[4,2], digits = 6),
                            int.p = round(summary(int.mdl2)$coefficients[4,4], digits = 4)) %>% 
      mutate(., swp.ps = case_when(
        swp.p < 0.0001 ~ '****',
        between(swp.p, 0.0001, 0.001) ~ '***',
        between(swp.p, 0.001, 0.01) ~ '**',
        between(swp.p, 0.01, 0.05) ~ '*',
        swp.p > 0.05 ~ 'n.s.'),
        temp.ps = case_when(
          temp.p < 0.0001 ~ '****',
          between(temp.p, 0.0001, 0.001) ~ '***',
          between(temp.p, 0.001, 0.01) ~ '**',
          between(temp.p, 0.01, 0.05) ~ '*',
          temp.p > 0.05 ~ 'n.s.'),
        int.ps = case_when(
          int.p < 0.0001 ~ '****',
          between(int.p, 0.0001, 0.001) ~ '***',
          between(int.p, 0.001, 0.01) ~ '**',
          between(int.p, 0.01, 0.05) ~ '*',
          int.p > 0.05 ~ 'n.s.')) %>% 
      relocate(swp.ps, .after = swp.p) %>% 
      relocate(temp.ps, .after = temp.p)
    
    i <- i+1
    if (t == length(cu.list)) {
      stat.cu.table1 <- do.call(bind_rows, stat.cu1)
      stat.cu.table2 <- do.call(bind_rows, stat.cu2)
      stat.cu.table <- bind_rows(stat.cu.table1, stat.cu.table2)}
  }
  return(stat.cu.table)
}

grw.under.cu <- grw.under.Cu.f()
write.csv(grw.under.cu, file = 'stat_grw_table2.csv')

# Reproduction affected by moisture and Cu effect under selected temperature
repr.under.temp.f <- function() {
  for (t in 1:length(temp.list)) {
    if (t == 1) {i <- 1;stat.temp <- vector("list", length = length(temp.list))}
    lgswp.mdl <- glm(repro ~ lgswp, data = subset(dat.org, temp == temp.list[t]), family = poisson(link = "log"))
    sqrcu.mdl <- glm(repro ~ sqrcu, data = subset(dat.org, temp == temp.list[t]), family = poisson(link = "log"))
    int.mdl <- glm(repro ~ lgswp*sqrcu, data = subset(dat.org, temp == temp.list[t]), family = poisson(link = "log"))
    
    stat.temp[[i]] <- tibble(temp = temp.list[t], 
                             swp.coef = round(summary(lgswp.mdl)$coefficients[2,1], digits = 2),
                             swp.lw = round(confint(lgswp.mdl)[2,1], digits = 2),
                             swp.up = round(confint(lgswp.mdl)[2,2], digits = 2),
                             swp.p = round(summary(lgswp.mdl)$coefficients[2,4], digits = 4),
                             cu.coef = round(summary(sqrcu.mdl)$coefficients[2,1], digits = 2),
                             cu.lw = round(confint(sqrcu.mdl)[2,1], digits = 2),
                             cu.up = round(confint(sqrcu.mdl)[2,2], digits = 2),
                             cu.p = round(summary(sqrcu.mdl)$coefficients[2,4], digits = 4),
                             int.coef = round(summary(int.mdl)$coefficients[4,1], digits = 2),
                             int.lw = round(confint(int.mdl)[4,1], digits = 2),
                             int.up = round(confint(int.mdl)[4,2], digits = 2),
                             int.p = round(summary(int.mdl)$coefficients[4,4], digits = 4)) %>% 
      mutate(., swp.ps = case_when(
        swp.p < 0.0001 ~ '****',
        between(swp.p, 0.0001, 0.001) ~ '***',
        between(swp.p, 0.001, 0.01) ~ '**',
        between(swp.p, 0.01, 0.05) ~ '*',
        swp.p > 0.05 ~ 'n.s.'),
        cu.ps = case_when(
          cu.p < 0.0001 ~ '****',
          between(cu.p, 0.0001, 0.001) ~ '***',
          between(cu.p, 0.001, 0.01) ~ '**',
          between(cu.p, 0.01, 0.05) ~ '*',
          cu.p > 0.05 ~ 'n.s.'),
        int.ps = case_when(
          int.p < 0.0001 ~ '****',
          between(int.p, 0.0001, 0.001) ~ '***',
          between(int.p, 0.001, 0.01) ~ '**',
          between(int.p, 0.01, 0.05) ~ '*',
          int.p > 0.05 ~ 'n.s.')) %>% 
      relocate(swp.ps, .after = swp.p) %>% 
      relocate(cu.ps, .after = cu.p)
    
    i <- i+1
    if (t == length(temp.list)) stat.temp.table <- do.call(bind_rows, stat.temp)
  }
  return(stat.temp.table)
}

repr.under.temp <- repr.under.temp.f()
write.csv(repr.under.temp, file = 'stat_repr_table1.csv')


# Reproduction affected by temp and lgswp under selected Cu levels
repr.under.Cu.f <- function(){
  for (t in 1:length(cu.list)) {
    if (t == 1) {i <- 1; stat.cu1 <- vector("list", length = length(cu.list)); stat.cu2 <- vector("list", length = length(cu.list))}
    
    lgswp.mdl1 <- glm(repro ~ lgswp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22), family = poisson(link = "log"))
    temp.mdl1 <- glm(repro ~ temp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22), family = poisson(link = "log"))
    int.mdl1 <- glm(repro ~ lgswp*temp, data = subset(dat.org, sqrcu == cu.list[t] & temp < 22), family = poisson(link = "log"))
    
    lgswp.mdl2 <- glm(repro ~ lgswp, data = subset(dat.org,  sqrcu == cu.list[t] & temp > 19), family = poisson(link = "log"))
    temp.mdl2 <- glm(repro ~ temp, data = subset(dat.org, sqrcu == cu.list[t] & temp > 19), family = poisson(link = "log"))
    int.mdl2 <- glm(repro ~ lgswp*temp, data = subset(dat.org, sqrcu == cu.list[t] & temp > 19), family = poisson(link = "log"))
    
    stat.cu1[[i]] <- tibble(temp = "<20",
                            sqrcu = cu.list[t], 
                            swp.coef = round(summary(lgswp.mdl1)$coefficients[2,1], digits = 2),
                            swp.lw = round(confint(lgswp.mdl1)[2,1], digits = 2),
                            swp.up = round(confint(lgswp.mdl1)[2,2], digits = 2),
                            swp.p = round(summary(lgswp.mdl1)$coefficients[2,4], digits = 4),
                            temp.coef = round(summary(temp.mdl1)$coefficients[2,1], digits = 2),
                            temp.lw = round(confint(temp.mdl1)[2,1], digits = 2),
                            temp.up = round(confint(temp.mdl1)[2,2], digits = 2),
                            temp.p = round(summary(temp.mdl1)$coefficients[2,4], digits = 4),
                            int.coef = round(summary(int.mdl1)$coefficients[4,1], digits = 2),
                            int.lw = round(confint(int.mdl1)[4,1], digits = 2),
                            int.up = round(confint(int.mdl1)[4,2], digits = 2),
                            int.p = round(summary(int.mdl1)$coefficients[4,4], digits = 4)) %>% 
      mutate(., swp.ps = case_when(
        swp.p < 0.0001 ~ '****',
        between(swp.p, 0.0001, 0.001) ~ '***',
        between(swp.p, 0.001, 0.01) ~ '**',
        between(swp.p, 0.01, 0.05) ~ '*',
        swp.p > 0.05 ~ 'n.s.'),
        temp.ps = case_when(
          temp.p < 0.0001 ~ '****',
          between(temp.p, 0.0001, 0.001) ~ '***',
          between(temp.p, 0.001, 0.01) ~ '**',
          between(temp.p, 0.01, 0.05) ~ '*',
          temp.p > 0.05 ~ 'n.s.'),
        int.ps = case_when(
          int.p < 0.0001 ~ '****',
          between(int.p, 0.0001, 0.001) ~ '***',
          between(int.p, 0.001, 0.01) ~ '**',
          between(int.p, 0.01, 0.05) ~ '*',
          int.p > 0.05 ~ 'n.s.')) %>% 
      relocate(swp.ps, .after = swp.p) %>% 
      relocate(temp.ps, .after = temp.p)
    
    stat.cu2[[i]] <- tibble(temp = ">20",
                            sqrcu = cu.list[t],
                            swp.coef = round(summary(lgswp.mdl2)$coefficients[2,1], digits = 2),
                            swp.lw = round(confint(lgswp.mdl2)[2,1], digits = 2),
                            swp.up = round(confint(lgswp.mdl2)[2,2], digits = 2),
                            swp.p = round(summary(lgswp.mdl2)$coefficients[2,4], digits = 4),
                            temp.coef = round(summary(temp.mdl2)$coefficients[2,1], digits = 2),
                            temp.lw = round(confint(temp.mdl2)[2,1], digits = 2),
                            temp.up = round(confint(temp.mdl2)[2,2], digits = 2),
                            temp.p = round(summary(temp.mdl2)$coefficients[2,4], digits = 4),
                            int.coef = round(summary(int.mdl2)$coefficients[4,1], digits = 2),
                            int.lw = round(confint(int.mdl2)[4,1], digits = 2),
                            int.up = round(confint(int.mdl2)[4,2], digits = 2),
                            int.p = round(summary(int.mdl2)$coefficients[4,4], digits = 4)) %>% 
      mutate(., swp.ps = case_when(
        swp.p < 0.0001 ~ '****',
        between(swp.p, 0.0001, 0.001) ~ '***',
        between(swp.p, 0.001, 0.01) ~ '**',
        between(swp.p, 0.01, 0.05) ~ '*',
        swp.p > 0.05 ~ 'n.s.'),
        temp.ps = case_when(
          temp.p < 0.0001 ~ '****',
          between(temp.p, 0.0001, 0.001) ~ '***',
          between(temp.p, 0.001, 0.01) ~ '**',
          between(temp.p, 0.01, 0.05) ~ '*',
          temp.p > 0.05 ~ 'n.s.'),
        int.ps = case_when(
          int.p < 0.0001 ~ '****',
          between(int.p, 0.0001, 0.001) ~ '***',
          between(int.p, 0.001, 0.01) ~ '**',
          between(int.p, 0.01, 0.05) ~ '*',
          int.p > 0.05 ~ 'n.s.')) %>% 
      relocate(swp.ps, .after = swp.p) %>% 
      relocate(temp.ps, .after = temp.p)
    
    i <- i+1
    if (t == length(cu.list)) {
      stat.cu.table1 <- do.call(bind_rows, stat.cu1)
      stat.cu.table2 <- do.call(bind_rows, stat.cu2)
      stat.cu.table <- bind_rows(stat.cu.table1, stat.cu.table2)}
  }
  return(stat.cu.table)
}

repr.under.cu <- repr.under.Cu.f()
write.csv(repr.under.cu, file = 'stat_repr_table2.csv')
