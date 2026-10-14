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
dat.org <- read_csv('cuaspb.csv') %>% as_tibble() %>% 
  mutate(., 
   sqrcu = case_when(
    cul == 'L1' ~ cu.list[1],
    cul == 'L2' ~ cu.list[2],
    cul == 'L3' ~ cu.list[3]
  ),
  tem = case_when(
    temp == '5' ~ temp.list[1],
    temp == '10' ~ temp.list[2],
    temp == '20' ~ temp.list[3],
    temp == '24' ~ temp.list[4],
    temp == '28' ~ temp.list[5]
  ),
  lgswp = case_when(
    swp == '3' ~ swp.list[1],
    swp == '5' ~ swp.list[2],
    swp == '10' ~ swp.list[3],
    swp == '40' ~ swp.list[4],
    swp == '70' ~ swp.list[5],
    swp == '100' ~ swp.list[6]
  )) %>% dplyr::select(tem, sqrcu, lgswp, weight, cu, pb, repli) %>% 
  filter(., between(cu, 0, 2000) & between(pb, 0, 2000))



dat.desc <- dat.org %>% 
  group_by(sqrcu, tem, lgswp) %>% 
  summarise(., cuave = mean(cu, na.rm = T),
            cuse = std.error(cu, na.rm = T),
            pbave = mean(pb, na.rm = T),
            pbse = std.error(pb, na.rm = T))

# Build model for visualisation
draw.line <- function() {
  lgswp.seq <- seq(swp.list[1], swp.list[6], 0.01)
  line.p <- vector("list", length = length(temp.list)*length(cu.list))
  for (t in 1:length(temp.list)) {
    if (t == 1) i <- 1
    for (c in 1:length(cu.list)) {
      dat.p1 <- subset(dat.org, sqrcu == cu.list[c] & tem == temp.list[t])
      cu.mdl <- lmrob(cu ~ lgswp, data = dat.p1)
      pb.mdl <- lmrob(pb ~ lgswp, data = dat.p1)

      line.p[[i]] <- tibble(lgswp = lgswp.seq, tem = temp.list[t], sqrcu = cu.list[c], 
                            cu = predict(cu.mdl, newdata = data.frame(lgswp = lgswp.seq), type = "response"),
                            pb = predict(pb.mdl, newdata = data.frame(lgswp = lgswp.seq), type = "response"))
                            
      i <- i+1
    }
  }
  
  line.t <- do.call(bind_rows, line.p)
  return(line.t)
}

line.t <- draw.line()

# Plot by group ------
p.cu <- ggplot()+
  geom_point(aes(x = lgswp, y = cuave, group = interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, size = 3, position = position_dodge(width = 0.1))+
  geom_errorbar(aes(lgswp, ymin = cuave-cuse, ymax = cuave+cuse, group =  interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, linewidth = 1, width = 0.2, position = position_dodge(width = 0.1))+
  geom_line(aes(x = lgswp, y = cu, color = factor(sqrcu)), data = line.t)+
  scale_color_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  scale_y_continuous(transform = 'log', breaks = c(1, 10, 100, 500, 1000, 3000))+
  labs(x = 'Soil Water potential /-kPa', y = 'Internal Cu conc. /mg kg-1', color = 'Cu level')+
  #ylim(c(0,200))+
  facet_wrap(~tem, nrow = 1)+
  theme_bw(base_family = "Aptos", base_size = 18)


p.pb <- ggplot()+
  geom_point(aes(x = lgswp, y = pbave, group = interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, size = 3, position = position_dodge(width = 0.1))+
  geom_errorbar(aes(lgswp, y = pbave, ymin = pbave-pbse, ymax = pbave+pbse, group =  interaction(lgswp, sqrcu), color = factor(sqrcu)),data = dat.desc, linewidth = 1, width = 0.2, position = position_dodge(width = 0.1))+
  geom_line(aes(x = lgswp, y = pb, color = factor(sqrcu)), data = line.t)+
  scale_color_manual(
    breaks = c('4.2','20.9','40.4'),
    values = c("#44ce9a", "#ffb870","#fb7070"),
    labels = c("Low","Medium", "High"))+
  scale_y_continuous(transform = 'log', breaks = c(1, 10, 100, 500, 1000, 2000))+
  labs(x = 'Soil Water potential /-kPa', y = 'Internal pb conc. /mg kg-1', color = 'Cu level')+
  facet_wrap(~tem, nrow = 1)+
  theme_bw(base_family = "Aptos", base_size = 18)


# Stat maybe
dat.lm <- dat.org %>% dplyr::select(., sqrcu, tem, lgswp, cu, pb) 


res_m0 <- lmrob(cu ~ tem+lgswp, data = subset(dat.lm, between(tem,20,30) & sqrcu == cu.list[1])); summary(res_m0);cbind(coef = coef(res_m0), ci = confint(res_m0))


res_m1 <- lmrob(cu ~ lgswp*sqrcu, data = subset(dat.lm, tem == temp.list[1])); summary(res_m1);cbind(coef = coef(res_m1), ci = confint(res_m1))


res_m2 <- lmrob(pb ~ tem*lgswp, data = subset(dat.lm, between(tem,0,22) & sqrcu == cu.list[3])); summary(res_m2);cbind(coef = coef(res_m2), ci = confint(res_m2))


res_m3 <- lmrob(pb ~ lgswp*sqrcu, data = subset(dat.lm, tem == temp.list[5])); summary(res_m3);cbind(coef = coef(res_m3), ci = confint(res_m3))

# m0 <- rlm(cu ~ lgswp+tem, data = subset(dat.lm, between(tem,0,21) & sqrcu == cu.list[1]));summary(m0)
# res_m0 <- residuals(m0)
# qqnorm(res_m0)
# qqline(res_m0, col = "red")
# 
# boot_model <- function(data, indices) {
#   d <- data[indices, ] # resample the data
#   coef(rlm(cu ~ lgswp+tem, data = d))
# }
# 
# # Bootstrapping the model
# results <- boot(subset(dat.lm, between(tem,0,21) & sqrcu == cu.list[1]), boot_model, R = 1000)
# boot.ci(results, type = "bca", index = 4)

m1 <- glm(cu ~ lgswp, data = subset(dat.lm, tem == temp.list[1] & sqrcu == cu.list[3]), family = gaussian());summary(m1)

m3 <- glm(cu ~ tem, data = subset(dat.lm, between(tem,20, 32) & sqrcu == cu.list[1] & lgswp == swp.list[1]), family = gaussian());summary(m3)
em3 <- emmeans(m3, ~ tem);pairs(em3, adjust = "tukey")

m2 <- lm(cu ~ as, data = dat.lm);summary(m2)

m4 <- glm(cu ~ sqrcu, data = subset(dat.lm, tem == temp.list[1] & lgswp == swp.list[3]));summary(m4)
em4 <- emmeans(m4, ~ sqrcu);contrast(em4, adjust = "bonferroni")


attach(subset(dat.lm, lgswp == swp.list[1] & sqrcu == cu.list[2] & between(tem,0, 21)))
pairwise.t.test(cu, tem, p.adjust.method = "bonf")
pairwise.t.test(pb, tem, p.adjust.method = "bonf")
pairwise.wilcox.test(cu, tem, p.adjust.method = "bonf")$p.value
pairwise.wilcox.test(pb, tem, p.adjust.method = "bonf")
detach()

kruskal.test(cu~tem, data = subset(dat.lm, lgswp == swp.list[1] & sqrcu == cu.list[1] & between(tem,0, 22)))
kwAllpairsNemenyiTest(cu~tem, data = subset(dat.lm, lgswp == swp.list[1] & sqrcu == cu.list[1] & between(tem,0, 22)))

