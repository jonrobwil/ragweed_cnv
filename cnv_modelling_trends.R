#############################################################################################################################################################################
################################################################   MERGED CNV TEMPORAL ANALYSIS   ###########################################################################
#############################################################################################################################################################################


library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(matrixStats)
library(lme4)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(circlize)
library(car)
library(emmeans)
library(ggeffects)
library(ppclust)
library(ggeffects)
library(pander)
library(kableExtra)
library(hash)


setwd("/home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30/figs/genotyping/modelling")

neutral <- as.matrix(read.table("neutral_hm.cov"))
samples <- read.table("samples_neut.list")
neutral.pcs <- eigen(neutral)
neutral.pcs <- as.data.frame(neutral.pcs$vectors[,1:2])
neutral.pcs <- cbind(samples, neutral.pcs)
colnames(neutral.pcs) <- c("V1", "neut_PC1", "neut_PC2")
rownames(neutral.pcs) = neutral.pcs[,1]
neutral.pcs[,1] = NULL

combined.pops=read.csv("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/combined_pops.csv")
rownames(combined.pops) = combined.pops[,1]
combined.pops[,1]=NULL

combined.pops[which(combined.pops$range == "eu" | combined.pops$range == "Europe"),][2] = "eu"
combined.pops[which(combined.pops$range == "East" | combined.pops$range == "middle" | combined.pops$range == "na" | combined.pops$range == "South" | combined.pops$range == "South-west" | combined.pops$range == "west" | combined.pops$range == "West"),][2] = "na"







cnv.name = hash()

cnv.name[["45530001"]] = "cnv-chr2a"
cnv.name[["30540001"]] = "cnv-chr4a"
cnv.name[["53800001"]] = "cnv-chr4b"
cnv.name[["37400001"]] = "cnv-chr9a"
cnv.name[["11570001"]] = "cnv-chr12a"
cnv.name[["18280001"]] = "cnv-chr14a"
cnv.name[["24890001"]] = "cnv-chr17b"
cnv.name[["21650001"]] = "cnv-chr18a"
cnv.name[["31240001"]] = "cnv-chr18b"
cnv.name[["12550001"]] = "cnv-chr5a"
cnv.name[["30730001"]] = "cnv-chr17c"
cnv.name[["10700001"]] = "cnv-chr8a"
cnv.name[["19190001"]] = "cnv-chr10a"
cnv.name[["21790001"]] = "cnv-chr13a"
cnv.name[["20190001"]] = "cnv-chr17a"
cnv.name[["29680001"]] = "cnv-chr6a"
cnv.name[["20640001"]] = "cnv-chr3a"








### cnv-chr4a: year:lat interaction ###


chr = "h1s4"
chr.start = 30540001
chr.end = 42390000


cnv.gt.mod = read.table("h1s4:30540001-42390000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s4:30540001-42390000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:lat + pc1, family=binomial)
Anova(glm, type = 3)
an_4a = Anova(glm, type = 3)
an_4a


#Response: cbind(h1count, h2count)
#         LR Chisq Df Pr(>Chisq)   
#year       5.0406  1    0.02476 * 
#range      3.9998  1    0.04551 * 
#lat        4.6797  1    0.03052 * 
#pc1        8.8819  1    0.00288 **
#year:lat   4.8990  1    0.02687 * 

# driven by southern american frequency shift


reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr4a"

write.table(reduced, "reduced_cont_cnv-chr4a.txt", quote=FALSE, sep="\t")




#get estimates of latitude in each range at specific timepoints
hist<-quantile(cnv.gt.group$year[cnv.gt.group$year< 2000])
early_hist<-hist[2]
med_hist<-hist[3]
late_hist<-hist[4]
med_modern<-median(cnv.gt.group$year[cnv.gt.group$year > 2000])

lat<-emtrends(glm, pairwise ~ year, var="lat", at=list(year=c(early_hist,med_hist,late_hist, med_modern)))
lat
# year lat.trend     SE  df asymp.LCL asymp.UCL
# 1892    0.0273 0.0466 Inf   -0.0640    0.1186
# 1906    0.0132 0.0422 Inf   -0.0696    0.0959
# 1930   -0.0111 0.0363 Inf   -0.0821    0.0600
# 2014   -0.0959 0.0397 Inf   -0.1738   -0.0181

lat<-as.data.frame(lat$emtrends)
write.table(lat, "cnv-chr4a-slopes-lat-emt.txt", quote=FALSE, sep="\t")


### plotting


cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )

cnv4a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)








### cnv-chr4b: range:lat interaction ###


chr = "h1s4"
chr.start = 53800001
chr.end = 54500000



cnv.gt.mod = read.table("h1s4:53800001-54500000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s4:53800001-54500000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + lat:range + year:range + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + lat:range + pc1, family=binomial)
an_4b = Anova(glm, type = 3)
an_4b


#Response: cbind(h1count, h2count)
#          LR Chisq Df Pr(>Chisq)    
#year         0.631  1   0.427055    
#range        4.040  1   0.044436 *  
#lat         37.290  1  1.018e-09 ***
#pc1         79.482  1  < 2.2e-16 ***
#range:lat    7.806  1   0.005206 ** 
#---


# stronger lat in NA


reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr4b"

write.table(reduced, "reduced_cont_cnv-chr4b.txt", quote=FALSE, sep="\t")





# emtrends 

means1_an_4b = (emtrends(glm, "range", var = "lat", adjust="fdr"))
means1_an_4b
rangec = as.data.frame(means1_an_4b)

write.table(rangec, "cnv-chr4b-range-emt.txt", quote=FALSE, sep="\t")

# range lat.trend     SE  df asymp.LCL asymp.UCL
# eu        0.117 0.0449 Inf    0.0165     0.218
# na        0.343 0.0705 Inf    0.1848     0.501
#
#Confidence level used: 0.95 
#Conf-level adjustment: bonferroni method for 2 estimates



cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv4b_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)








### cnv-chr5a - year:range + year:lat interactions ###




chr = "h1s5"
chr.start = 12550001
chr.end = 13470000


cnv.gt.mod = read.table("h1s5:12550001-13470000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s5:12550001-13470000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:lat + year:range + pc1, family=binomial)
Anova(glm, type = 3)

an_5a = Anova(glm, type = 3)
an_5a



#Response: cbind(h1count, h2count)
#           LR Chisq Df Pr(>Chisq)   
#year         3.9917  1   0.045726 * 
#range        6.6830  1   0.009734 **
#lat          4.4550  1   0.034800 * 
#pc1          2.9177  1   0.087613 . 
#year:lat     4.1577  1   0.041445 * 
#year:range   6.5605  1   0.010427 * 




reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr5a"

write.table(reduced, "reduced_cont_cnv-chr5a.txt", quote=FALSE, sep="\t")



means1_an_5a = (emtrends(glm, "range", var = "year", adjust="fdr"))
means1_an_5a
# range year.trend      SE  df asymp.LCL asymp.UCL
# eu       0.00763 0.00389 Inf   -0.0011  0.016365
# na      -0.00900 0.00420 Inf   -0.0184  0.000417

rangec = as.data.frame(means1_an_5a)

write.table(rangec, "cnv-chr5a-range-emt.txt", quote=FALSE, sep="\t")



#get estimates of latitude in each range at specific timepoints
hist<-quantile(cnv.gt.group$year[cnv.gt.group$year< 2000])
early_hist<-hist[2]
med_hist<-hist[3]
late_hist<-hist[4]
med_modern<-median(cnv.gt.group$year[cnv.gt.group$year > 2000])

lat<-emtrends(glm, pairwise ~ year, var="lat", at=list(year=c(early_hist,med_hist,late_hist, med_modern)))
lat
# year lat.trend     SE  df asymp.LCL asymp.UCL
# 1892    0.1538 0.0498 Inf    0.0562    0.2513
# 1906    0.1384 0.0438 Inf    0.0526    0.2242
# 1930    0.1120 0.0350 Inf    0.0434    0.1805
# 2014    0.0196 0.0384 Inf   -0.0556    0.0948

lat<-as.data.frame(lat$emtrends)
write.table(lat, "cnv-chr5a-lat-emt.txt", quote=FALSE, sep="\t")



cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv5a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)






### cnv-chr8a - range:lat and year:range ###


chr = "h1s8"
chr.start = 10700001
chr.end = 12580000




cnv.gt.mod = read.table("h1s8:10700001-12580000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n*2)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n*2) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s8:10700001-12580000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n*2)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n*2) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)


glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + range:lat + year:range + pc1, family=binomial)
an_8a = Anova(glm, type = 3)
an_8a


#Response: cbind(h1count, h2count)
#           LR Chisq Df Pr(>Chisq)   
#year         1.0137  1   0.314011   
#range        3.1130  1   0.077672 . 
#lat          0.4618  1   0.496770   
#pc1          3.8335  1   0.050238 . 
#range:lat    5.0557  1   0.024544 * 
#year:range   6.8422  1   0.008903 **



reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr8a"

write.table(reduced, "reduced_cont_cnv-chr8a.txt", quote=FALSE, sep="\t")





# emtrends 

means1_an_8a = (emtrends(glm, "range", var = "lat", adjust="fdr"))
means1_an_8a

# range lat.trend     SE  df asymp.LCL asymp.UCL
# eu      -0.0395 0.0391 Inf  -0.12710    0.0481
# na       0.0726 0.0312 Inf   0.00265    0.1426

rangec = as.data.frame(means1_an_8a)

write.table(rangec, "cnv-chr8a-latrange-emt.txt", quote=FALSE, sep="\t")


means1_an_8a = (emtrends(glm, "range", var = "year", adjust="fdr"))
means1_an_8a

# range year.trend      SE  df asymp.LCL asymp.UCL
# eu       0.00681 0.00246 Inf   0.00129   0.01233
# na      -0.00304 0.00284 Inf  -0.00940   0.00333

rangec = as.data.frame(means1_an_8a)

write.table(rangec, "cnv-chr8a-yearrange-emt.txt", quote=FALSE, sep="\t")



cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv8a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)







### cnv-chr9a no interactions ###


chr = "h1s9"
chr.start = 37400001
chr.end = 38820000


cnv.gt.mod = read.table("h1s9:37400001-38820000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n*2)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n*2) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s9:37400001-38820000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n*2)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n*2) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + lat:year + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + pc1, family=binomial)
an_9a = Anova(glm, type = 3)
an_9a



#Analysis of Deviance Table (Type III tests)
#
#Response: cbind(h1count, h2count)
#      LR Chisq Df Pr(>Chisq)  
#year    4.1191  1    0.04240 *
#range   5.1728  1    0.02294 *
#lat     3.2198  1    0.07275 .
#pc1     1.8278  1    0.17639  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr9a"

write.table(reduced, "reduced_cont_cnv-chr9a.txt", quote=FALSE, sep="\t")





cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )



cnv9a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)





### cnv-chr10a - range*lat*year interactions ###



chr = "h1s10"
chr.start = 19190001
chr.end = 19660000



cnv.gt.mod = read.table("h1s10:19190001-19660000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s10:19190001-19660000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

an_10a = Anova(glm, type = 3)
an_10a

reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr10a"

write.table(reduced, "reduced_cont_cnv-chr10a.txt", quote=FALSE, sep="\t")


#Response: cbind(h1count, h2count)
#               LR Chisq Df Pr(>Chisq)    
#year              3.550  1   0.059529 .  
#range             7.352  1   0.006697 ** 
#lat               2.458  1   0.116897    
#pc1              99.764  1  < 2.2e-16 ***
#year:range        7.430  1   0.006414 ** 
#year:lat          2.386  1   0.122417    
#range:lat         7.045  1   0.007949 ** 
#year:range:lat    7.151  1   0.007493 ** 



me = ggemmeans(glm, terms = c("lat", "year", "range"))
me = as.data.frame(me)

#get estimates of latitude in each range at specific timepoints
hist<-quantile(cnv.gt.group$year[cnv.gt.group$year< 2000])
early_hist<-hist[2]
med_hist<-hist[3]
late_hist<-hist[4]
med_modern<-median(cnv.gt.group$year[cnv.gt.group$year > 2000])

lat<-emtrends(glm, pairwise ~ year, var="lat", at=list(range="na", year=c(early_hist,med_hist,late_hist, med_modern)))
lat<-as.data.frame(lat$emtrends)
lat$range<-"na"

lat2<-emtrends(glm, pairwise ~ year, var="lat", at=list(range="eu", year=c(early_hist,med_hist,late_hist, med_modern)))
lat2<-as.data.frame(lat2$emtrends)
lat2$range<-"eu"
lat3<-rbind(lat,lat2)

#  year    lat.trend         SE  df     asymp.LCL  asymp.UCL range
#1 1892 -0.008260829 0.07224355 Inf -0.1498555888 0.13333393    na
#2 1906  0.000148943 0.06346489 Inf -0.1242399535 0.12453784    na
#3 1930  0.014565696 0.04982788 Inf -0.0830951560 0.11222655    na
#4 2014  0.065024331 0.04497151 Inf -0.0231182005 0.15316686    na
#5 1892  0.151535445 0.04792285 Inf  0.0576083833 0.24546251    eu
#6 1906  0.121804300 0.04220501 Inf  0.0390840007 0.20452460    eu
#7 1930  0.070836622 0.03635260 Inf -0.0004131691 0.14208641    eu
#8 2014 -0.107550249 0.06266512 Inf -0.2303716218 0.01527112    eu


write.table(lat3, "cnv-chr10a-rangelat-emt.txt", quote=FALSE, sep="\t")


#get estimate of time in each  range at specific timepoints
summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])
mina<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[1]
lowa<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[2]
mida<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[3]
higha<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[5]
maxa<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[6]
emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(lowa, mida, higha)))
emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(mina, lowa, mida, higha, maxa)))

#  lat year.trend      SE  df asymp.LCL asymp.UCL
# 25.1   -0.00503 0.01203 Inf  -0.02860   0.01854
# 38.5    0.00300 0.00354 Inf  -0.00394   0.00993
# 41.1    0.00457 0.00308 Inf  -0.00146   0.01060
# 44.1    0.00638 0.00390 Inf  -0.00126   0.01402
# 49.8    0.00980 0.00729 Inf  -0.00448   0.02409

year1<-emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(mina, lowa, mida, higha, maxa)))
year1<-as.data.frame(year1$emtrends)
year1$range<-"na"


mine<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[1]
lowe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[2]
mide<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[3]
highe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[5]
maxe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[6]
emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(lowe, mide, highe)))
emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(mine, lowe, mide, highe, maxe)))

#  lat year.trend      SE  df asymp.LCL asymp.UCL
# 43.1   0.010964 0.00391 Inf  0.003302   0.01863
# 45.7   0.005470 0.00262 Inf  0.000329   0.01061
# 48.5  -0.000626 0.00220 Inf -0.004935   0.00368
# 52.2  -0.008299 0.00362 Inf -0.015404  -0.00119
# 63.3  -0.031969 0.01073 Inf -0.053005  -0.01093




year2<-emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(mine, lowe, mide, highe, maxe)))
year2<-as.data.frame(year2$emtrends)
year2$range<-"eu"

year<-rbind(year1, year2)
write.table(year, "cnv-chr10a-rangeyear-emt.txt", quote=FALSE, sep="\t")





cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )

cnv10a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)









### cnv-chr14a - year*range*lat interactions ###




chr = "h1s14"
chr.start = 18280001
chr.end = 18750000



cnv.gt.mod = read.table("h1s14:18280001-18750000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s14:18280001-18750000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
an_14a = Anova(glm, type = 3)
an_14a

#Analysis of Deviance Table (Type III tests)
#
#Response: cbind(h1count, h2count)
#               LR Chisq Df Pr(>Chisq)  
#year             5.3510  1    0.02071 *
#range            5.5352  1    0.01864 *
#lat              4.7772  1    0.02884 *
#pc1              5.0555  1    0.02455 *
#year:range       5.6976  1    0.01699 *
#year:lat         5.1444  1    0.02332 *
#range:lat        5.4662  1    0.01939 *
#year:range:lat   5.5768  1    0.01820 *
#---


reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr14a"

write.table(reduced, "reduced_cont_cnv-chr14a.txt", quote=FALSE, sep="\t")




me = ggemmeans(glm, terms = c("lat", "year", "range"))
me = as.data.frame(me)

#get estimates of latitude in each range at specific timepoints
hist<-quantile(cnv.gt.group$year[cnv.gt.group$year< 2000])
early_hist<-hist[2]
med_hist<-hist[3]
late_hist<-hist[4]
med_modern<-median(cnv.gt.group$year[cnv.gt.group$year > 2000])

lat<-emtrends(glm, pairwise ~ year, var="lat", at=list(range="na", year=c(early_hist,med_hist,late_hist, med_modern)))
lat<-as.data.frame(lat$emtrends)
lat$range<-"na"

lat2<-emtrends(glm, pairwise ~ year, var="lat", at=list(range="eu", year=c(early_hist,med_hist,late_hist, med_modern)))
lat2<-as.data.frame(lat2$emtrends)
lat2$range<-"eu"
lat3<-rbind(lat,lat2)
lat3

#  year   lat.trend         SE  df   asymp.LCL asymp.UCL range
#1 1892  0.08005183 0.07311655 Inf -0.06325397 0.2233576    na
#2 1906  0.07951228 0.06440575 Inf -0.04672066 0.2057452    na
#3 1930  0.07858736 0.05091066 Inf -0.02119571 0.1783704    na
#4 2014  0.07535011 0.04574806 Inf -0.01431444 0.1650146    na
#5 1892 -0.06246054 0.09597295 Inf -0.25056407 0.1256430    eu
#6 1906 -0.01532830 0.08429628 Inf -0.18054599 0.1498894    eu
#7 1930  0.06546982 0.07089900 Inf -0.07348967 0.2044293    eu
#8 2014  0.34826325 0.11384993 Inf  0.12512149 0.5714050    eu

write.table(lat3, "cnv-chr14a-lattrend-emt.txt", quote=FALSE, sep="\t")


#get estimate of time in each  range at specific timepoints
summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])
mina<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[1]
lowa<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[2]
mida<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[3]
higha<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[5]
maxa<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[6]
emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(lowa, mida, higha)))
emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(mina, lowa, mida, higha, maxa)))

#  lat year.trend      SE  df asymp.LCL asymp.UCL
# 25.1   5.39e-04 0.01229 Inf  -0.02356   0.02463
# 38.5   2.36e-05 0.00404 Inf  -0.00790   0.00794
# 41.1  -7.74e-05 0.00357 Inf  -0.00707   0.00692
# 44.1  -1.93e-04 0.00423 Inf  -0.00848   0.00809
# 49.8  -4.13e-04 0.00739 Inf  -0.01490   0.01407

year1<-emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(mina, lowa, mida, higha, maxa)))
year1<-as.data.frame(year1$emtrends)
year1$range<-"na"


mine<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[1]
lowe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[2]
mide<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[3]
highe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[5]
maxe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[6]
emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(lowe, mide, highe)))
emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(mine, lowe, mide, highe, maxe)))

#  lat year.trend      SE  df asymp.LCL asymp.UCL
# 43.1  -0.017952 0.00948 Inf  -0.03653  0.000625
# 45.7  -0.009243 0.00668 Inf  -0.02233  0.003843
# 48.5   0.000421 0.00454 Inf  -0.00848  0.009322
# 52.2   0.012585 0.00546 Inf   0.00189  0.023277
# 63.3   0.050109 0.01827 Inf   0.01430  0.085918

year2<-emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(mine, lowe, mide, highe, maxe)))
year2<-as.data.frame(year2$emtrends)
year2$range<-"eu"

year<-rbind(year1, year2)
write.table(year, "cnv-chr14a-yeartrend-emt.txt", quote=FALSE, sep="\t")



cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv14a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)







### cnv-chr17a - range*lat*year interactions ###





chr = "h1s17"
chr.start = 20190001
chr.end = 20740000


cnv.gt.mod = read.table("h1s17:20190001-20740000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s17:20190001-20740000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

an_17a = Anova(glm, type = 3)
an_17a


#Response: cbind(h1count, h2count)
#               LR Chisq Df Pr(>Chisq)    
#year             0.0026  1  0.9590881    
#range            6.6995  1  0.0096438 ** 
#lat              0.3477  1  0.5553995    
#pc1             14.2199  1  0.0001626 ***
#year:range       6.1297  1  0.0132927 *  
#year:lat         0.1970  1  0.6571746    
#range:lat        8.0968  1  0.0044345 ** 
#year:range:lat   7.4237  1  0.0064370 ** 


reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr17a"

write.table(reduced, "reduced_cont_cnv-chr17a.txt", quote=FALSE, sep="\t")





me = ggemmeans(glm, terms = c("lat", "year", "range"))
me = as.data.frame(me)

#get estimates of latitude in each range at specific timepoints
hist<-quantile(cnv.gt.group$year[cnv.gt.group$year< 2000])
early_hist<-hist[2]
med_hist<-hist[3]
late_hist<-hist[4]
med_modern<-median(cnv.gt.group$year[cnv.gt.group$year > 2000])

lat<-emtrends(glm, pairwise ~ year, var="lat", at=list(range="na", year=c(early_hist,med_hist,late_hist, med_modern)))
lat<-as.data.frame(lat$emtrends)
lat$range<-"na"

lat2<-emtrends(glm, pairwise ~ year, var="lat", at=list(range="eu", year=c(early_hist,med_hist,late_hist, med_modern)))
lat2<-as.data.frame(lat2$emtrends)
lat2$range<-"eu"
lat3<-rbind(lat,lat2)

#  year    lat.trend         SE  df   asymp.LCL   asymp.UCL range
#1 1892 -0.341810991 0.07835058 Inf -0.49537531 -0.18824668    na
#2 1906 -0.321510273 0.06925895 Inf -0.45725531 -0.18576523    na
#3 1930 -0.286709043 0.05443914 Inf -0.39340780 -0.18001028    na
#4 2014 -0.164904738 0.03459752 Inf -0.23271464 -0.09709484    na
#5 1892  0.041108144 0.04258742 Inf -0.04236167  0.12457796    eu
#6 1906  0.026714743 0.03747777 Inf -0.04674033  0.10016982    eu
#7 1930  0.002040342 0.03272382 Inf -0.06209717  0.06617786    eu
#8 2014 -0.084320061 0.05966003 Inf -0.20125157  0.03261145    eu


write.table(lat3, "cnv-chr17a-lat-year-emt.txt", quote=FALSE, sep="\t")


#get estimate of time in each  range at specific timepoints
summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])
mina<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[1]
lowa<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[2]
mida<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[3]
higha<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[5]
maxa<-summary(cnv.gt.group$lat[cnv.gt.group$range=="na"])[6]
emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(lowa, mida, higha)))
emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(mina, lowa, mida, higha, maxa)))

#  lat year.trend      SE  df asymp.LCL asymp.UCL
# 25.1   -0.01442 0.01090 Inf  -0.03578   0.00695
# 38.5    0.00496 0.00322 Inf  -0.00135   0.01126
# 41.1    0.00876 0.00325 Inf   0.00239   0.01513
# 44.1    0.01313 0.00441 Inf   0.00448   0.02178
# 49.8    0.02139 0.00786 Inf   0.00598   0.03680

year1<-emtrends(glm, pairwise ~ lat, var="year", at=list(range="na", lat=c(mina, lowa, mida, higha, maxa)))
year1<-as.data.frame(year1$emtrends)
year1$range<-"na"


mine<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[1]
lowe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[2]
mide<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[3]
highe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[5]
maxe<-summary(cnv.gt.group$lat[cnv.gt.group$range=="eu"])[6]
emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(lowe, mide, highe)))
emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(mine, lowe, mide, highe, maxe)))

#   lat year.trend      SE  df asymp.LCL asymp.UCL
#  43.1    0.00442 0.00377 Inf  -0.00297   0.01181
#  45.7    0.00176 0.00259 Inf  -0.00331   0.00683
#  48.5   -0.00119 0.00215 Inf  -0.00541   0.00303
#  52.2   -0.00491 0.00339 Inf  -0.01155   0.00174
#  63.3   -0.01636 0.00995 Inf  -0.03586   0.00313

year2<-emtrends(glm, pairwise ~ lat, var="year", at=list(range="eu", lat=c(mine, lowe, mide, highe, maxe)))
year2<-as.data.frame(year2$emtrends)
year2$range<-"eu"

year<-rbind(year1, year2)
write.table(year, "cnv-chr17a-year-range-emt.txt", quote=FALSE, sep="\t")





cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )
cnv17a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)






### cnv-chr17b - range:lat interactions ###



chr = "h1s17"
chr.start = 24890001
chr.end = 26890000


cnv.gt.mod = read.table("h1s17:24890001-26890000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s17:24890001-26890000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + range:lat + pc1, family=binomial)
an_17b = Anova(glm, type = 3)
an_17b


#Response: cbind(h1count, h2count)
#          LR Chisq Df Pr(>Chisq)    
#year        0.0165  1   0.897661    
#range      10.1860  1   0.001415 ** 
#lat         6.7214  1   0.009527 ** 
#pc1        19.9724  1  7.857e-06 ***
#range:lat   9.8237  1   0.001723 **


reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr17b"

write.table(reduced, "reduced_cont_cnv-chr17b.txt", quote=FALSE, sep="\t")




means1_an_17b = (emtrends(glm, "range", var = "lat", adjust="fdr"))
means1_an_17b
# range lat.trend     SE  df asymp.LCL asymp.UCL
# eu      -0.0143 0.0312 Inf   -0.0843    0.0557
# na       0.1178 0.0275 Inf    0.0562    0.1793

rangec = as.data.frame(means1_an_17b)

write.table(rangec, "cnv-chr17b-range-emt.txt", quote=FALSE, sep="\t")


cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )

cnv17b_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)






### cnv-chr17c - range:lat interactions ###


chr = "h1s17"
chr.start = 30730001
chr.end = 31210000





cnv.gt.mod = read.table("h1s17:30730001-31210000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n*2)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n*2) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s17:30730001-31210000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n*2)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n*2) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + range:lat + year:range + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + range:lat + pc1, family=binomial)
an_17c = Anova(glm, type = 3)
an_17c

#Response: cbind(h1count, h2count)
#          LR Chisq Df Pr(>Chisq)    
#year        7.7484  1  0.0053759 ** 
#range       6.5533  1  0.0104688 *  
#lat         5.2877  1  0.0214767 *  
#pc1        12.7703  1  0.0003522 ***
#range:lat   5.6110  1  0.0178478 * 




reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr17c"

write.table(reduced, "reduced_cont_cnv-chr17c.txt", quote=FALSE, sep="\t")




# emtrends 

means1_an_17c = (emtrends(glm, "range", var = "lat", adjust="fdr"))
means1_an_17c

# range lat.trend     SE  df asymp.LCL asymp.UCL
# eu      0.10544 0.0324 Inf    0.0327     0.178
# na     -0.00462 0.0311 Inf   -0.0743     0.065

rangec = as.data.frame(means1_an_17c)

write.table(rangec, "cnv-chr17c-range-emt.txt", quote=FALSE, sep="\t")




cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv17c_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)





### cnv-chr18a - year:range + year:lat interactions ###


chr = "h1s18"
chr.start = 21650001
chr.end = 22250000


cnv.gt.mod = read.table("h1s18:21650001-22250000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s18:21650001-22250000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + pc1, family=binomial)
Anova(glm, type = 3)

an_18a = Anova(glm, type = 3)
an_18a


#Response: cbind(h1count, h2count)
#           LR Chisq Df Pr(>Chisq)  
#year         5.4899  1    0.01913 *
#range        5.8531  1    0.01555 *
#lat          5.6295  1    0.01766 *
#pc1          4.8289  1    0.02799 *
#year:range   5.7652  1    0.01635 *
#year:lat     5.4212  1    0.01989 *



reduced<-as.data.frame(t(Anova(glm, type = 3)))
reduced$HB<-"cnv-chr18a"

write.table(reduced, "reduced_cont_cnv-chr18a.txt", quote=FALSE, sep="\t")




means1_an_18a = (emtrends(glm, "range", var = "year", adjust="fdr"))
means1_an_18a
# range year.trend      SE  df asymp.LCL asymp.UCL
# eu       0.00756 0.00359 Inf -0.000482   0.01561
# na      -0.00708 0.00397 Inf -0.015970   0.00181

rangec = as.data.frame(means1_an_18a)

write.table(rangec, "cnv-chr18a-range-emt.txt", quote=FALSE, sep="\t")



#get estimates of latitude in each range at specific timepoints
hist<-quantile(cnv.gt.group$year[cnv.gt.group$year< 2000])
early_hist<-hist[2]
med_hist<-hist[3]
late_hist<-hist[4]
med_modern<-median(cnv.gt.group$year[cnv.gt.group$year > 2000])

lat<-emtrends(glm, pairwise ~ year, var="lat", at=list(year=c(early_hist,med_hist,late_hist, med_modern)))
lat
# year lat.trend     SE  df asymp.LCL asymp.UCL
# 1892    0.1254 0.0458 Inf    0.0357    0.2151
# 1906    0.1097 0.0403 Inf    0.0308    0.1887
# 1930    0.0829 0.0321 Inf    0.0201    0.1457
# 2014   -0.0110 0.0333 Inf   -0.0763    0.0542

lat<-as.data.frame(lat$emtrends)
write.table(lat, "cnv-chr18a-lat-emt.txt", quote=FALSE, sep="\t")



cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )

cnv18a_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)






png("cnv_predicted_sig.png", height=2000, width=1500)
plot_grid(cnv4a_model, cnv4b_model, cnv5a_model, cnv8a_model, cnv9a_model, cnv10a_model, cnv14a_model, cnv17a_model, cnv17b_model, cnv17c_model, cnv18a_model, ncol = 2, align = "vh")
dev.off()





































##############################################
### junk ###
##############################################

chr = "h1s4"
chr.start = 53800001
chr.end = 54500000



cnv.gt.mod = read.table("h1s4:53800001-54500000_mod_gt.txt", header = T)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1] = NULL

cnv.gt.mod = merge(cnv.gt.mod, combined.pops, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod = merge(cnv.gt.mod, neutral.pcs, by = 0)
rownames(cnv.gt.mod) = cnv.gt.mod[,1]
cnv.gt.mod[,1]=NULL

cnv.gt.mod.group = cnv.gt.mod %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.mod.group$h1count <- (cnv.gt.mod.group$n)*cnv.gt.mod.group$cont
cnv.gt.mod.group$h2count <- (cnv.gt.mod.group$n) - cnv.gt.mod.group$h1count




cnv.gt.hist = read.table("h1s4:53800001-54500000_hist_gt.txt", header = T)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1] = NULL

cnv.gt.hist = merge(cnv.gt.hist, combined.pops, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist = merge(cnv.gt.hist, neutral.pcs, by = 0)
rownames(cnv.gt.hist) = cnv.gt.hist[,1]
cnv.gt.hist[,1]=NULL

cnv.gt.hist.group = cnv.gt.hist %>%
  group_by(lon, lat, range, year) %>%
  summarise(cont = mean(gt), n = n(), pc1 = mean(neut_PC1), pc2 = mean(neut_PC2))

cnv.gt.hist.group$h1count <- (cnv.gt.hist.group$n)*cnv.gt.hist.group$cont
cnv.gt.hist.group$h2count <- (cnv.gt.hist.group$n) - cnv.gt.hist.group$h1count


cnv.gt.hist.group <- as.data.frame(cnv.gt.hist.group)
cnv.gt.mod.group <- as.data.frame(cnv.gt.mod.group)

cnv.gt.group = rbind(cnv.gt.mod.group,cnv.gt.hist.group)

options(contrasts = c("contr.sum", "contr.poly"))

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + year:range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:range + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + lat:range + year:range + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + lat:range + pc1, family=binomial)
an_4b = Anova(glm, type = 3)



means1_an_4b = (emtrends(glm, "range", var = "lat", adjust="fdr"))
summary(means1_an_4b)
means2_an_4b = emtrends(glm, pairwise ~ range, var="lat")






means<-emmeans(glm
               , type = "pairs"
               , specs =  c("range", "lat")
               , by =  c("range", "lat")
               , adjust = 'fdr')
as.data.frame(regrid(means))
pairs(regrid(means), simple = "each", combine = TRUE, adjust="fdr")







means2_an_4b = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=year)

options(contrasts = c("contr.sum", "contr.poly"))
# all interactions
scaf27a.glm<-glm(data = scaf27a_merged, cbind(counth1, counth2) ~ year * range * latitude + V1mean , family=binomial)
Anova(scaf27a.glm, type = 3)
scaf27a.glm<-glm(data = scaf27a_merged, cbind(counth1, counth2) ~ year + range + latitude + year:range + year:latitude + range:latitude + V1mean , family=binomial)
Anova(scaf27a.glm, type = 3)
scaf27a.glm<-glm(data = scaf27a_merged, cbind(counth1, counth2) ~ year + range + latitude + year:range + range:latitude + V1mean , family=binomial)
Anova(scaf27a.glm, type = 3)
scaf27a.glm<-glm(data = scaf27a_merged, cbind(counth1, counth2) ~ year + range + latitude + year:range  + V1mean , family=binomial)
Anova(scaf27a.glm, type = 3)
early_hist<-median(scaf27$year[scaf27$year<1905])
late_hist<-median(scaf27$year[scaf27$year>1904 & scaf27$year < 2000])
med_hist<-median(scaf27$year[scaf27$year<2000])
med_modern<-median(scaf27$year[scaf27$year>1999])
early_hist<-median(scaf27a_merged$year[scaf27a_merged$year<1905])
late_hist<-median(scaf27a_merged$year[scaf27a_merged$year>1904 & scaf27a_merged$year < 2000])
med_hist<-median(scaf27a_merged$year[scaf27a_merged$year<2000])
med_modern<-median(scaf27a_merged$year[scaf27a_merged$year>1999])
early_hist<-median(scaf27a$year[scaf27a$year<1905])
late_hist<-median(scaf27a$year[scaf27a$year>1904 & scaf27a$year < 2000])
med_hist<-median(scaf27a$year[scaf27a$year<2000])
med_modern<-median(scaf27a$year[scaf27a$year>1999])
emtrends(scaf27a.glm, pairwise ~ range, var="year")
# analyse with PC1
scaf27a.me = ggemmeans(scaf27a.glm, terms = c("latitude", "year", "range"))
scaf27a.me = as.data.frame(scaf27a.me)
reduced27a<-as.data.frame(t(Anova(scaf27a.glm, type = 3)))
reduced27a$HB<-"HB27a"
write.table(reduced27a, "~/Dropbox/Documents/haploblock_modelling_kh/reduced_cont27a.txt", quote=FALSE, sep="\t")
summary27a<-as.data.frame(t(summary(scaf27a.glm)$coefficients))
write.table(summary27a, "~/Dropbox/Documents/haploblock_modelling_kh/summary_cont27a.txt", quote=FALSE, sep="\t")
#compare with single sample
scaf27a.glm<-glm(data = scaf27a, cbind(h1count, h2count) ~ year * range * latitude + V1 , family=binomial)
Anova(scaf27a.glm, type = 3)
lat=median(cnv.gt.group$latitude[which(cnv.gt.group$range=="Europe")])
lat
ggemmeans(glm, terms = c("time", "range", "latitude [lat]"), cov.reduce=range)
lat=median(cnv.gt.group$latitude[which(cnv.gt.group$range=="North America")])
lat
ggemmeans(glm, terms = c("time", "range", "latitude [lat]"), cov.reduce=range)
11:14
For instance these are
11:14
> lat=median(cnv.gt.group$latitude[which(cnv.gt.group$range=="Europe")])
> lat
[1] 48.82232
> ggemmeans(glm, terms = c("time", "range", "latitude [lat]"), cov.reduce=range)
# Predicted probabilities of cbind(h1count, h2count)
#    range = North America
time     | Predicted |       95% CI
-----------------------------------
modern   |      0.69 | [0.60, 0.76]
historic |      0.37 | [0.30, 0.45]
11:14
lat=median(cnv.gt.group$latitude[which(cnv.gt.group$range=="North America")])
> lat
[1] 41.17518
> ggemmeans(glm, terms = c("time", "range", "latitude [lat]"), cov.reduce=range)
# Predicted probabilities of cbind(h1count, h2count)
#    range = North America
time     | Predicted |       95% CI
-----------------------------------
modern   |      0.54 | [0.48, 0.60]
historic |      0.33 | [0.27, 0.40]









#emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")

means2_an_4b = emmeans(glm
               , type = "pairs"
               , specs =  c("year", "lat")
               , by =  c("year", "lat")
               , adjust = 'fdr')


cnv3.lm.CONT = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv3.lm.CONT = as.data.frame(cnv3.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv3.predicted.na = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," NA")) +
  labs(x = "",
    y = "") +
  xlim(20, 50) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv3.predicted.eu = ggplot() +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv3.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle(paste0(cnv.name[[as.character(chr.start)]]," EU")) +
  labs(x = "",
    y = "") +
  xlim(40, 70) +
  ylim(0, 1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv4b_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)
