




library(tidyr)
library(data.table)
library(dplyr)
library(stringr)
library(Hmisc)
library(ggplot2)
library(gridExtra)
library(circlize)
library(matrixStats)
library(lme4)
library(cowplot)
library(car)
library(emmeans)
library(ggeffects)
library(ppclust)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(emmeans)
library(ggeffects)
library(pafr)



setwd("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/cnv_results_final")





############################################################################################################################################################
###################################################################### NORMALISED ##########################################################################
############################################################################################################################################################


# RUN ONCE
#
#
file.names = dir(path = "depth_out_mod_Q30", pattern = "_all.depth", full.names = TRUE)

# make df of average depth for each gene and sample
window_depth = read.table(file.names[1], header = T)
window_depth$normdepth = window_depth$depth/median(window_depth$depth)
window_depth.df = cbind(window_depth[1], window_depth[2], window_depth[4])

for (i in 2:length(file.names)) {
	TEMP = read.table(file.names[i], header = T)
	TEMP$normdepth = TEMP$depth/median(TEMP$depth)
	window_depth.df = cbind(window_depth.df, TEMP[,4])
}

window_depth.df = cbind(window_depth[1], window_depth[2], window_depth[4])

for (i in 2:length(file.names)) {
	TEMP = read.table(file.names[i], header = T)
	TEMP$normdepth = TEMP$depth/median(TEMP$depth)
	window_depth.df = cbind(window_depth.df, TEMP[,4])
}



# name columns with samples
column.names = str_replace_all(file.names, "_all.depth", "")
column.names = str_replace_all(column.names, "depth_out_mod_Q30/", "")

colnames(window_depth.df)[3:ncol(window_depth.df)] = column.names

write.table(window_depth.df, file = "10kb_windowed_depth_mod_Q30_nrml.txt", quote = F, row.names = F)



############################################################################################################################################################
###################################################################### UNNORMALISED ########################################################################
############################################################################################################################################################

# RUN ONCE
#
#
file.names = dir(path = "depth_out_mod_Q30", pattern = "_all.depth", full.names = TRUE)

# make df of average depth for each gene and sample
window_depth = read.table(file.names[1], header = T)
window_depth$normdepth = window_depth$depth/10000
window_depth.df = cbind(window_depth[1], window_depth[2], window_depth[4])

for (i in 2:length(file.names)) {
	TEMP = read.table(file.names[i], header = T)
	TEMP$normdepth = TEMP$depth/10000
	window_depth.df = cbind(window_depth.df, TEMP[,4])
}

# name columns with samples
column.names = str_replace_all(file.names, "_all.depth", "")
column.names = str_replace_all(column.names, "depth_out_mod_Q30/", "")

colnames(window_depth.df)[3:ncol(window_depth.df)] = column.names

write.table(window_depth.df, file = "10kb_windowed_depth_mod_Q30_unnrml.txt", quote = F, row.names = F)






############################################################################################################################################################
###################################################################### FILTERING ###########################################################################
############################################################################################################################################################




window_depth.filt = read.table("10kb_windowed_depth_mod_Q30_nrml.txt", header = T)
bla = str_replace_all(colnames(window_depth.filt), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")
colnames(window_depth.filt) = bla
window_depth.filt$window=window_depth.filt$window

rownames(window_depth.filt) = paste0(window_depth.filt$ID, ":", window_depth.filt$window)

window_depth.filt[,1:2]= NULL

window_depth.filt.t = as.data.frame(t(window_depth.filt))

sample.info = read.csv("ragweed_sample_info.csv")
neutral.pcs.env=read.table("neutral-pc-env.txt")

#window_depth.filt.t=window_depth.filt.t[which(rownames(window_depth.filt.t) %in% rownames(neutral.pcs.env)),]

total.wind=ncol(window_depth.filt.t)
ind.numb=nrow(window_depth.filt.t)

#window_depth.filt = as.data.frame(t(window_depth.filt.t))
#
#window_depth.filt$windsd = NA
#window_depth.filt$counts= NA
#
#for (i in 1:total.wind){
#	window_depth.filt$windsd[i] = sd(as.numeric(as.matrix(window_depth.filt[i,1:ind.numb])))
#	#window_depth.filt.t.EU$counts[i] = sum(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])) > mean(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])))+1.5*(as.numeric(as.matrix(window_depth.filt.t.EU$windsd[i]))))
#	#window_depth.filt.t.EU$counts[i] = window_depth.filt.t.EU$counts[i]+sum(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])) < mean(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])))-1.5*(as.numeric(as.matrix(window_depth.filt.t.EU$windsd[i]))))
#}
#
#for (i in 1:total.wind){
#	#window_depth.filt.t.EU$windsd[i] = sd(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])))
#	window_depth.filt$counts[i] = sum(as.numeric(as.matrix(window_depth.filt[i,1:ind.numb])) > mean(as.numeric(as.matrix(window_depth.filt[i,1:ind.numb])))+2*(as.numeric(as.matrix(window_depth.filt$windsd[i]))))
#	#window_depth.filt.t.EU$counts[i] = window_depth.filt.t.EU$counts[i]+sum(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])) < mean(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])))-1.5*(as.numeric(as.matrix(window_depth.filt.t.EU$windsd[i]))))
#}
#
#for (i in 1:total.wind){
#	#window_depth.filt.t.EU$windsd[i] = sd(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])))
#	#window_depth.filt.t.EU$counts[i] = sum(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])) > mean(as.numeric(as.matrix(window_depth.filt.t.EU[i,1:ind.numb])))+1.5*(as.numeric(as.matrix(window_depth.filt.t.EU$windsd[i]))))
#	window_depth.filt$counts[i] = window_depth.filt$counts[i]+sum(as.numeric(as.matrix(window_depth.filt[i,1:ind.numb])) < mean(as.numeric(as.matrix(window_depth.filt[i,1:ind.numb])))-2*(as.numeric(as.matrix(window_depth.filt$windsd[i]))))
#}
#
#window_depth.filt.maf=window_depth.filt[which(window_depth.filt$counts > 16), ]
#
#write.table(window_depth.filt.maf, file = "10kb_windowed_depth_2SD5_modern_Q30.txt", quote = F, row.names = T)
#write.table(rownames(window_depth.filt.maf), file = "maf_2SD5_windows_Q30.list", quote = F, row.names = T)
#
#
#
#window_depth.filt.UP=window_depth.filt
#for (i in 1:total.wind){
#	window_depth.filt.UP$counts[i] = sum(as.numeric(as.matrix(window_depth.filt.UP[i,1:ind.numb])) > mean(as.numeric(as.matrix(window_depth.filt.UP[i,1:ind.numb])))+2*(as.numeric(as.matrix(window_depth.filt.UP$windsd[i]))))
#}
#
#window_depth.filt.UP.maf=window_depth.filt.UP[which(window_depth.filt.UP$counts > 16), ]
#
#
#
#window_depth.filt.DOWN=window_depth.filt
#for (i in 1:total.wind){
#	window_depth.filt.DOWN$counts[i] = sum(as.numeric(as.matrix(window_depth.filt.DOWN[i,1:ind.numb])) < mean(as.numeric(as.matrix(window_depth.filt.DOWN[i,1:ind.numb])))-2*(as.numeric(as.matrix(window_depth.filt.DOWN$windsd[i]))))
#}
#
#window_depth.filt.DOWN.maf=window_depth.filt.DOWN[which(window_depth.filt.DOWN$counts > 16), ]




##### PLOT VARIANCE #####


output=matrix(NA,nrow=0,ncol=2)
for (i in 1:total.wind) {
	# make temp coverage df for each window as well as trait value
	#wd.temp=cbind(window.means.NA.temp[i], window_depth.traits$V2, window_depth.traits[pc1])
	# calculate lm between window and trait
	print(paste0("window",i , "of", 2))
	windvar = var(window_depth.filt.t[,i])
	#among = as.data.frame(VarCorr(res))[1,4]
	#within = as.data.frame(VarCorr(res))[2,4]
	#qst=among/(among+2*(within))
	out.line=c(colnames(window_depth.filt.t[i]), windvar)
	## add lm for trait and window to output, then move onto next window
	output=rbind(output,out.line)
}

write.table(as.data.frame(output), file = "window_variance.txt", quote = F, row.names = F, col.names = F)




pheno=1
trait.associations = read.table("window_variance.txt")
#### manhattan of r^2 ####
# read in chromosome data for manhattan plots
contigs = read.table("~/om62_scratch/ragweed2022/ref/chromo-length.txt", header = F)
contigs$V3 = c(1:nrow(contigs))
colnames(contigs) = c("CONTIG", "LENGTH", "CHR")
contigs$CHR = as.numeric(contigs$CHR)

# get breaks for manhattan plots
manbreaks = c()
s = 0
for (i in 1:nrow(contigs)) {
	nbp = contigs[contigs$CHR == i,]$LENGTH
	manbreaks = c(manbreaks, round(nbp/2) + s)
	s = s + nbp
}
manbreaksnames = 1:nrow(contigs)
colnames(trait.associations)=c("window", "var")
trait.associations = trait.associations %>% separate(window, c("CONTIG", "LOC"), ":")
trait.associations$LOC = as.numeric(trait.associations$LOC)
trait.associations$var = as.numeric(as.character(trait.associations$var))
gwas = merge(contigs, trait.associations, by = "CONTIG")
#gwas.sig=0.05/37202
# add a cumulative BP count
gwas$LOCc = NA
s = 0
for (i in 1:length(unique(gwas$CHR))) {
	nbp = max(gwas[gwas$CHR == i,]$LENGTH)
	gwas[gwas$CHR == i, "LOCc"] = gwas[gwas$CHR == i,"LOC"] + s
	s = s + nbp
}
#top.one=round(0.005*nrow(gwas))
#n=0.05
#gwas.sig = head(gwas[order(abs(gwas$qst),decreasing=T),],.05*nrow(gwas))
#gwas.sig.i = gwas.sig[,which(colnames(gwas.sig) == "rsq" | colnames(gwas.sig) == "LOCc")]
assign(paste0("env.plot_", pheno), ggplot() +
	geom_point(data = gwas, aes(x = LOCc, y = var, color = as.factor(CHR)), size = 0.5, shape = 19) +
	#geom_point(data = gwas.sig, aes(x = LOCc, y = var), color = "#f37ded", size = 0.5, shape = 19) +
	ggtitle("variance of depth windows") +
	scale_color_manual(values = rep(c("#bdbdbd", "#636363"), 500)) +
	scale_x_continuous(name = "chromosome", breaks = manbreaks, labels = manbreaksnames) +
	labs(y = "variance") +
	#geom_hline(yintercept = -log10(gwas.sig), color = "black") +
	ylim(0, 25) +
	theme_classic() +
	theme(legend.position = "none"))



png("window_depth_variance.png", width = 1400, height= 300)
grid.arrange(env.plot_1, ncol = 1)
dev.off()














##### MAF FILTERING #####


window_depth.filt.new = as.data.frame(t(window_depth.filt.t))

window_depth.filt.new$maf=NA

for (i in 1:nrow(window_depth.filt.new)){
wind.mean = mean(as.numeric(window_depth.filt.new[i, ]), na.rm = T)
wind.sd = sd(as.numeric(window_depth.filt.new[i, ]), na.rm = T)
wind.upper = wind.mean + 2 * wind.sd
wind.lower = wind.mean - 2 * wind.sd
wind.upper.maf = length(which(as.numeric(window_depth.filt.new[i, ]) > wind.upper)) / sum(!is.na(as.numeric(window_depth.filt.new[i, ])))
wind.lower.maf = length(which(as.numeric(window_depth.filt.new[i, ]) < wind.lower)) / sum(!is.na(as.numeric(window_depth.filt.new[i, ])))
wind.maf = max(wind.upper.maf, wind.lower.maf)
window_depth.filt.new$maf[i] = wind.maf
}

window_depth.filt.new.maf=window_depth.filt.new[which(window_depth.filt.new$maf >= 0.05),]

write.table(window_depth.filt.new.maf, file = "window_depth_2SD5_modern_Q30.txt", quote = F, row.names = T)
write.table(rownames(window_depth.filt.new.maf), file = "window_depth_2SD5_modern_Q30.list", quote = F, row.names = T)




while read line4 line1 line2 line3; do cat window_depth_Q30_maf.list | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$1 && start >= $2 && start <= $3 || chr==$1 && end >= $2 && end <= $3 || chr==$1 && start >= $2 && end <= $3 || chr==$1 && start <= $2 && end >= $3)print $0}';  done < cnvs.txt
while read line1 line2 line3; do cat window_depth_Q30_maf.list | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$1 && start >= $2 && start <= $3 || chr==$1 && end >= $2 && end <= $3 || chr==$1 && start >= $2 && end <= $3 || chr==$1 && start <= $2 && end >= $3)print $0}';  done < putative_merged_cnvs.list




#############################################################################################################################################################################
#######################################################################  CHR4-INS FIG  ######################################################################################
#############################################################################################################################################################################

cd /home/jwilson/om62_scratch/ragweed2022/lostruct/angsd_analysis

##### find better chr4 coordinates

sbatch -J h1s4:30540001-36000001 -e h1s4:30540001-36000001.e -o h1s4:30540001-36000001.o pcangsd_MOD_array.sh h1s4:30540001-36000001

cd /home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/cnv_results_final

cp /home/jwilson/om62_scratch/ragweed2022/lostruct/angsd_analysis/h1s4:30540001-36000001.cov .

vim h1s4:30540001-36000001.bed
h1s4 30540001 36000001 5460000

mkdir depth_out_h1s4:30540001-36000001

while read sample; do bash chr4_split.sh ${sample}; done < all_modern_samples.list

cd depth_out_h1s4\:30540001-36000001/

while read line; do cat ${line}_h1s.depth | tail -n 1 | awk -v sample=${line} '{print sample, $0}' >> all_chr4del_split_first_histmod.depth; done < ../all_samples.list












setwd("/home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/cnv_results_final")



# MDS plots

mdsdf = read.table("h1s4.mds", header = T)
mdsdf$LOC=mdsdf$WIND*100000
i=1
mds.plot = ggplot() +
	geom_point(data = mdsdf,
		aes(x = LOC, y = mdsdf[, (i + 2)]), colour="darkgrey", size = 1, shape = 19) +
	ggtitle("") +
	labs(x = "", y = paste0("MDS", i)) +
	ylim(-0.5, 1) + 
	theme_classic() +
	theme(legend.position = "none",
	plot.title = element_text(size = 14),
	axis.text = element_text(size = 12),
	axis.title = element_text(size = 14)
	)

#png(file = "MDS_PLOT_TEST.png",
#	width = 1920,
#	height = 3 * 240)
#plot_grid(mds.plot)
#dev.off()


# hap1-hap2 alignment
hap1hap2 = read_paf("~/om62_scratch/ragweed2022/assembly/hap1hap2.chrs.filtered.paf", include_tags = FALSE)

hb4 = subset(hap1hap2, qname == "h1s4" & tname == "h2s6")
#hb4 = subset(hap1hap2, qname == "h1s7" & tname == "h2s7")
#hb8 = subset(hap1hap2, qname == "h1s8" & tname == "h2s8")
#hb11 = subset(hap1hap2, qname == "h1s11" & tname == "h2s11")
#hb16b = subset(hap1hap2, qname == "h1s16" & tname == "h2s15")

hb4.dp = dotplot(hb4,
	label_seqs = FALSE,
	xlab = "haplotype 1 chromosome 4",
	ylab = "haplotype 2 chromosome 4",
	line_size = 5) +
	theme_classic() +
	theme(legend.position = "none",
		plot.title = element_text(size = 14),
		axis.text = element_text(size = 12),
		axis.title = element_text(size = 14)
	)


chr.name = "h1s4"


merged.all.chrs.df = read.table("~/ha22_scratch/lowerstruct-ragweed/depth/window_depth_2SD5_modern_Q30_merged_0.6_1Mb.bed")
colnames(merged.all.chrs.df) = c("chrom", "start", "end")
merged.all.chrs.df$size = merged.all.chrs.df$end - merged.all.chrs.df$start


merging = ggplot() +
geom_rect(data = subset(merged.all.chrs.df, chrom == chr.name), aes(xmin = start, xmax = end, , ymin = 0, ymax = size), color = "black", fill = NA) +
labs(x = "", y = "") +
xlim(0, max(subset(merged.all.chrs.df, chrom == chr.name)$end) + 5000) +
ylim(0, max(merged.all.chrs.df$size)) +
theme(legend.position = "none",
plot.title = element_text(size = 36, face="bold"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
legend.text = element_text(color = "black", size = 14)) +
ggtitle("")

col1 = plot_grid(mds.plot, merging, rel_heights = c(10, 10), rel_widths = c(10, 10), align = "v", ncol=1)

png("chr4_fig.png", width=800, height=400)
plot_grid(col1, hb4.dp, rel_heights = c(10, 20), rel_widths = c(20, 10), align = "vh", ncol=2)
dev.off()



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


setwd("/home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30/figs/genotyping/modelling")

neutral <- as.matrix(read.table("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30_histmod/neutral_reg.cov"))
samples <- read.table("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30_histmod/all_samples.list")
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



chr = "h1s4"
chr.start = 30540001
chr.end = 42390000

window_depth = read.table("../../putative_cnv_normalised_depth.txt", header = T)
samples=read.table("../../modern_samples.list")

bla = str_replace_all(colnames(window_depth), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")
colnames(window_depth) = bla

rownames(window_depth) = window_depth[,1]

window_depth[,1] = NULL
window_depth = window_depth[,which(colnames(window_depth) %in% samples$V1)]

current_depth = window_depth[which(rownames(window_depth) == paste0(chr, ":", chr.start, "-", chr.end)),]
current_depth.t = as.data.frame(t(current_depth))
C <- as.matrix(read.table(paste0("../../", chr, ":", chr.start, "-", chr.end, "_mod.cov")))
e <- eigen(C)

pc <- e$vectors
pcs <- as.data.frame(pc[,1:2])
names(pcs)[names(pcs) == "V1"] <- "PC1"
names(pcs)[names(pcs) == "V2"] <- "PC2"

pcs.samples=as.data.frame(cbind(samples$V1,pcs))
rownames(pcs.samples)=pcs.samples[,1]
pcs.samples[,1]=NULL


pcs.samples=merge(pcs.samples, current_depth.t, by=0)
rownames(pcs.samples) = pcs.samples[,1]
pcs.samples[,1] = NULL
colnames(pcs.samples) = c("PC1", "PC2", "depth")


# Optionally set colours using RColorBrewer
cols = brewer.pal(9, "RdPu")
# Define colour pallete
#pal = colorRampPalette(c("blue", "red"))
# Use the following line with RColorBrewer
pal = colorRampPalette(cols)
# Rank variable for colour assignment
pcs.samples$V2=as.numeric(as.character(pcs.samples$depth))
pcs.samples$order = findInterval(pcs.samples$depth, sort(pcs.samples$depth))

pcs.samples$PC1=as.numeric(as.character(pcs.samples$PC1))

#PC1.depth=ggplot() +
#	geom_point(data = pcs.samples, aes(x = depth, y = PC1), colour="orange", size = 1, shape = 19) +
#	#geom_point(data = trait.associations, aes(x = LOC1, y = -log10(V2)), colour="darkgrey", size = 1, shape = 19) +
#	stat_smooth(method = "lm") +
#	ggtitle(paste0("../", chr, ":", chr.start, "-", chr.end)) +
#	labs(x = "Read Depth", y = "PC1") +
#	theme_classic() +
#	theme(legend.position = "none",
#	plot.title = element_text(size = 14),
#	axis.text = element_text(size = 12),
#	axis.title = element_text(size = 14)
#	)	
	



set.seed(42)

h1s4_kcm = cbind(pcs.samples$PC1, pcs.samples$depth)

res.fcm <- fcm(h1s4_kcm[,1:2], centers=2, nstart = 20)
res.fcm2 <- ppclust2(res.fcm, "kmeans")

samples.gt = as.data.frame(cbind(rownames(pcs.samples), res.fcm2$cluster))
tmax = table(samples.gt$V2)
max.clus = names(tmax[tmax==max(tmax)])

tmin = table(samples.gt$V2)
min.clus = names(tmin[tmin==min(tmin)])

samples.gt$V1 = as.character(samples.gt$V1)
samples.gt$V2 = as.numeric(as.character(samples.gt$V2))
samples.gt$gt = NA
samples.gt[which(samples.gt$V2 == max.clus),]$gt = 0
samples.gt[which(samples.gt$V2 == min.clus),]$gt = 1
samples.gt.final = as.data.frame(cbind(samples.gt$V1, samples.gt$gt))
colnames(samples.gt.final) = c("sample", "gt")



mod.pc.df = as.data.frame(cbind(h1s4_kcm, res.fcm2$cluster))





mod.plot = ggplot() +
  geom_point(data = mod.pc.df, aes(x = V3, y = V4, color = as.factor(samples.gt$gt)), size = 5, shape = 19) +
  labs(x = "", y = "") +
theme(legend.position = "none",
plot.title = element_text(size = 36, face="bold"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
panel.background = element_blank(),
legend.title = element_blank(), 
legend.text = element_text(color = "black", size = 14)) +
ggtitle("")


png("/home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/mod_cnv_gt_final_review.png", width=400, height=400, bg = "transparent")
plot_grid(mod.plot)
dev.off()



















combined.pops=read.csv("../combined_pops.csv")
rownames(combined.pops) = combined.pops[,1]
combined.pops[,1]=NULL
combined.pops.mod = combined.pops[which(combined.pops$type == "modern"),]
combined.pops.hist = combined.pops[which(combined.pops$type == "historic"),]

combined.pops[which(combined.pops$range == "eu" | combined.pops$range == "Europe"),][2] = "eu"
combined.pops[which(combined.pops$range == "East" | combined.pops$range == "middle" | combined.pops$range == "na" | combined.pops$range == "South" | combined.pops$range == "South-west" | combined.pops$range == "west" | combined.pops$range == "West"),][2] = "na"

cnv4split = read.table("all_chr4del_split_first_histmod.depth", header = F)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL

sample_medians=read.table("sample_medians_histmod.list")
rownames(sample_medians) = sample_medians[,1]
sample_medians[,1] = NULL

cnv4split = merge(cnv4split, sample_medians, by = 0)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL

cnv4split$depthnrml = cnv4split$V5/cnv4split$V2.y
cnv4split=merge(cnv4split, combined.pops, by = 0)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL
cnv4aa=cnv4split
cnv4aa[which(cnv4aa$type == "historic"),]$depthnrml = cnv4aa[which(cnv4aa$type == "historic"),]$depthnrml+0.052

stl.h = read.table("stlouis-h.list")
stl.m = read.table("stlouis-m.list")


cnv4aa.h = cnv4aa[which(rownames(cnv4aa) %in% stl.h$V1),]
cnv4aa.m = cnv4aa[which(rownames(cnv4aa) %in% stl.m$V1),]

cnv4aa.sl = rbind(cnv4aa.m, cnv4aa.h)

chr4ins.stl=ggplot() +
	geom_point(data = cnv4aa.sl, aes(y = depthnrml, x = year), colour="orange", size = 1, shape = 19) +
	#geom_point(data = trait.associations, aes(x = LOC1, y = -log10(V2)), colour="darkgrey", size = 1, shape = 19) +
	stat_smooth(method = "lm") +
	ggtitle("") +
	labs(y = "Read Depth", x = "Year") +
	theme_classic() +
	theme(legend.position = "none",
	plot.title = element_text(size = 14),
	axis.text = element_text(size = 12),
	axis.title = element_text(size = 14)
	)	
	
	
col1 = plot_grid(chr4_his, PC1.depth, chr4ins.stl, ncol = 1, align = "v", axis = "l", rel_heights = c(4, 7, 7), labels = c("", ""), label_size = 36, vjust = 1)

col2 = plot_grid(df.gg, hb4.dp, mds.plot, ncol = 1, align = "v", axis = "l", rel_heights = c(4, 7, 7), labels = c("", ""), label_size = 36, vjust = 1)



png(file = "MDS_PLOT_TEST2.png",
	width = 1000,
	height = 1000)
plot_grid(col2, col1, align = "h", axis = "b")
dev.off()


png(file = "MDS_PLOT_TEST.png",
	width = 1920/2,
	height = 3 * 240)
plot_grid(col1, col2, align = "h", axis = "b")
dev.off()








### finding weird region on chr4
chr4weird = mdsdf[which(mdsdf$LOC > 37000001 & mdsdf$LOC<40000001),]
h1s4:37500000-39900000


vim h1s4:37500000-39900000.bed
h1s4 37500000 39900000 2400000

mkdir depth_out_h1s4:37500000-39900000

while read sample; do bash chr4_split_weird.sh ${sample}; done < all_modern_samples.list
while read sample; do bash chr4_split_weird_historic.sh ${sample}; done < historic_samples.list

cd depth_out_h1s4:37500000-39900000/

while read line; do cat ${line}_h1s.depth | tail -n 1 | awk -v sample=${line} '{print sample, $0}' >> all_chr4del_split_weird_histmod.depth; done < ../all_modern_samples.list





#window_depth.filt = read.table("../10kb_windowed_depth_normalised_filtered_100pc_ALL.txt", header = T)
#window_depth.filt[,1:2]=NULL
#window_depth.filt.t=as.data.frame(t(window_depth.filt))
#
#bla = str_replace_all(rownames(window_depth.filt.t), "\\.", "-")
#bla = str_replace_all(bla, "X28", "28")
#rownames(window_depth.filt.t) = bla
#
#which(colnames(window_depth.filt.t)=="h1s4:30000001")
#which(colnames(window_depth.filt.t)=="h1s4:40000001")
#
#chr4.del.wind=window_depth.filt.t[,22811:23811]
#
#chr4.del.mean=as.data.frame(cbind(rownames(chr4.del.wind), rowMeans(chr4.del.wind)))














#############################################################################################################################################################################
########################################################################   QST ANALYSIS   ###################################################################################
#############################################################################################################################################################################




windows.depth=read.table("10kb_windowed_depth_mod_Q30_unnrml.txt", header = T)

bla = str_replace_all(colnames(windows.depth), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")
colnames(windows.depth) = bla

rownames(windows.depth) = paste0(windows.depth$ID, ":", windows.depth$window)
windows.depth[,1:2] = NULL

#total.wind=nrow(windows.depth)
#ind.numb=ncol(windows.depth)

windows.depth.t = as.data.frame(t(windows.depth))

EU.env=read.table("EU-baypass-samples-pops-env.txt")
rownames(EU.env) <- EU.env[,1]
EU.env[,1] <- NULL

pop.info= as.data.frame(cbind(rownames(EU.env), as.character(EU.env$V2)))
rownames(pop.info)=pop.info[,1]
pop.info[,1] <- NULL

windows.EU = merge(windows.depth.t, pop.info, by=0)
rownames(windows.EU)=windows.EU[,1]
windows.EU[,1] <- NULL

sample_medians = read.table("sample_medians.list")
rownames(sample_medians)=sample_medians[,1]
sample_medians[,1] <- NULL


pops.temp = windows.EU %>% group_by(V2) %>% filter(n() > 2)
pops = pops.temp$V2
pops.u = unique(pops)
windows.EU.pops = windows.EU[which(windows.EU$V2 %in% pops.u),]
window.number=ncol(windows.EU.pops)-1

windows.EU.pops = merge(windows.EU.pops, sample_medians, by = 0)
rownames(windows.EU.pops)=windows.EU.pops[,1]
windows.EU.pops[,1] <- NULL



populations=cbind(rownames(windows.EU.pops), as.character(windows.EU.pops$V2))
write.table(populations, file="EU-qst-pops.list", quote=F, row.names=F, col.names=F)
#V2=windows.EU.pops$V2

filtered.windows=read.table("window_depth_2SD5_modern_Q30.list")

windows.depth.maf=windows.EU.pops[,which(colnames(windows.EU.pops) %in% filtered.windows$x)]

eek=cbind(windows.depth.maf, windows.EU.pops$V2.x, windows.EU.pops$V2.y)
names(eek)[ncol(eek)] <- "sample.cvg"
names(eek)[ncol(eek)-1] <- "V2"
windows.EU.pops.maf=eek
window.number=ncol(windows.EU.pops.maf)-2

output=matrix(NA,nrow=0,ncol=2)
for (i in 1:window.number) {
	# make temp coverage df for each window as well as trait value
	#wd.temp=cbind(window.means.NA.temp[i], window_depth.traits$V2, window_depth.traits[pc1])
	# calculate lm between window and trait
	print(paste0("window",i , "of", window.number))
	res = lmer(windows.EU.pops.maf[,i] ~ 1 + (1|V2) + sample.cvg, data = windows.EU.pops.maf)
	among = as.data.frame(VarCorr(res))[1,4]
	within = as.data.frame(VarCorr(res))[2,4]
	qst=among/(among+2*(within))
	out.line=c(colnames(windows.EU.pops.maf[i]), qst)
	# add lm for trait and window to output, then move onto next window
	output=rbind(output,out.line)
}

write.table(as.data.frame(output), file = "qst-EU.txt", quote = F, row.names = F, col.names = F)


qst.vals.EU = read.table("qst-EU.txt")
rownames(qst.vals.EU) = qst.vals.EU[,1]
qst.vals.EU[,1] = NULL

qst.vals.EU.t = as.data.frame(t(qst.vals.EU))


neut.fst.EU=read.table("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/neut-fst/EU_all.weir.fst", header=T)

eee=as.data.frame(t(qst.vals.EU.t))

neut.fst.EU.pc5 = head(neut.fst.EU[order(abs(neut.fst.EU$WEIR_AND_COCKERHAM_FST),decreasing=T),],.05*nrow(neut.fst.EU))
neut.fst.EU.pc1 = head(neut.fst.EU[order(abs(neut.fst.EU$WEIR_AND_COCKERHAM_FST),decreasing=T),],.01*nrow(neut.fst.EU))

pc5.fst = min(neut.fst.EU.pc5$WEIR_AND_COCKERHAM_FST)
pc1.fst = min(neut.fst.EU.pc1$WEIR_AND_COCKERHAM_FST)



pheno=1
trait.associations = as.data.frame(cbind(rownames(qst.vals.EU), qst.vals.EU[,pheno]))
#### manhattan of r^2 ####
# read in chromosome data for manhattan plots
contigs = read.table("~/om62_scratch/ragweed2022/ref/chromo-length.txt", header = F)
contigs$V3 = c(1:nrow(contigs))
colnames(contigs) = c("CONTIG", "LENGTH", "CHR")
contigs$CHR = as.numeric(contigs$CHR)

# get breaks for manhattan plots
manbreaks = c()
s = 0
for (i in 1:nrow(contigs)) {
	nbp = contigs[contigs$CHR == i,]$LENGTH
	manbreaks = c(manbreaks, round(nbp/2) + s)
	s = s + nbp
}
manbreaksnames = 1:nrow(contigs)
colnames(trait.associations)=c("window", "qst")
trait.associations = trait.associations %>% separate(window, c("CONTIG", "LOC"), ":")
trait.associations$LOC = as.numeric(trait.associations$LOC)
trait.associations$qst = as.numeric(as.character(trait.associations$qst))
gwas = merge(contigs, trait.associations, by = "CONTIG")
#gwas.sig=0.05/37202
# add a cumulative BP count
gwas$LOCc = NA
s = 0
for (i in 1:length(unique(gwas$CHR))) {
	nbp = max(gwas[gwas$CHR == i,]$LENGTH)
	gwas[gwas$CHR == i, "LOCc"] = gwas[gwas$CHR == i,"LOC"] + s
	s = s + nbp
}
#top.one=round(0.005*nrow(gwas))
#n=0.05
gwas.sig.q = head(gwas[order(abs(gwas$qst),decreasing=T),],.05*nrow(gwas))
gwas.sig.1pc = head(gwas[order(abs(gwas$qst),decreasing=T),],.01*nrow(gwas))
pc5.qst = min(gwas.sig.q$qst)
pc1.qst = min(gwas.sig.1pc$qst)

gwas.sig = gwas[which(abs(gwas$qst) >= pc1.fst),]

#gwas.sig.i = gwas.sig[,which(colnames(gwas.sig) == "rsq" | colnames(gwas.sig) == "LOCc")]
assign(paste0("env.plot_", pheno), ggplot() +
	geom_point(data = gwas, aes(x = LOCc, y = abs(qst), color = as.factor(CHR)), size = 0.5, shape = 19) +
	geom_point(data = gwas.sig, aes(x = LOCc, y = abs(qst)), color = "#ffacfb", size = 0.5, shape = 19) +
	ggtitle("EU") +
	scale_color_manual(values = rep(c("#bdbdbd", "#636363"), 500)) +
	scale_x_continuous(name = "chromosome", breaks = manbreaks, labels = manbreaksnames) +
	labs(y = "QST") +
	#geom_hline(yintercept = -log10(gwas.sig), color = "black") +
	ylim(0, 1) +
	theme_classic() +
	theme(legend.position = "none"))
write.table(as.data.frame(cbind(as.character(gwas.sig[,1]), gwas.sig[,4])), file = "qst_EU_1fst.txt", quote = F, row.names = F, col.names = F)



qst.EU.hist = ggplot() +
		geom_histogram(data = eee, aes(x = V2), color = "black", fill = "white", bins = 50) +
		#geom_vline(data = eee, aes(xintercept = pc5.fst), color="#35cff9") +
		geom_vline(data = eee, aes(xintercept = pc1.fst), color="#ffacfb") +
		ggtitle("QST EU") +
		xlim(0,1) +
		labs(xlim = "Qst") +
		theme_classic()

fst.EU.hist = ggplot() +
		geom_histogram(data = neut.fst.EU, aes(x = WEIR_AND_COCKERHAM_FST), color = "black", fill = "white", bins = 50) +
		ggtitle("FST EU") +
		xlim(0,1) +
		labs(x = "Fst") +
		#geom_vline(data = neut.fst.EU, aes(xintercept = pc5.fst), color="#35cff9") +
		geom_vline(data = neut.fst.EU, aes(xintercept = pc1.fst), color="#ffacfb") +
		theme_classic()

qst.fst = plot_grid(qst.EU.hist, fst.EU.hist, nrow = 2, align = "v", labels = c("", ""), label_size = 36)


png("qst_EU_final.png", width = 1000, height= 400)
plot_grid(env.plot_1, qst.fst, ncol = 2, align = "h", rel_widths = c(3, 1), labels = c("", ""), label_size = 36)
dev.off()










windows.depth=read.table("10kb_windowed_depth_mod_Q30_unnrml.txt", header = T)

bla = str_replace_all(colnames(windows.depth), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")
colnames(windows.depth) = bla

rownames(windows.depth) = paste0(windows.depth$ID, ":", windows.depth$window)
windows.depth[,1:2] = NULL

#total.wind=nrow(windows.depth)
#ind.numb=ncol(windows.depth)

windows.depth.t = as.data.frame(t(windows.depth))

NA.env=read.table("NA-baypass-samples-pops-env.txt")
rownames(NA.env) <- NA.env[,1]
NA.env[,1] <- NULL

pop.info= as.data.frame(cbind(rownames(NA.env), as.character(NA.env$V2)))
rownames(pop.info)=pop.info[,1]
pop.info[,1] <- NULL

windows.NA = merge(windows.depth.t, pop.info, by=0)
rownames(windows.NA)=windows.NA[,1]
windows.NA[,1] <- NULL

sample_medians = read.table("sample_medians.list")
rownames(sample_medians)=sample_medians[,1]
sample_medians[,1] <- NULL


pops.temp = windows.NA %>% group_by(V2) %>% filter(n() > 2)
pops = pops.temp$V2
pops.u = unique(pops)
windows.NA.pops = windows.NA[which(windows.NA$V2 %in% pops.u),]
window.number=ncol(windows.NA.pops)-1

windows.NA.pops = merge(windows.NA.pops, sample_medians, by = 0)
rownames(windows.NA.pops)=windows.NA.pops[,1]
windows.NA.pops[,1] <- NULL



populations=cbind(rownames(windows.NA.pops), as.character(windows.NA.pops$V2))
write.table(populations, file="NA-qst-pops.list", quote=F, row.names=F, col.names=F)
#V2=windows.NA.pops$V2

filtered.windows=read.table("window_depth_2SD5_modern_Q30.list")

windows.depth.maf=windows.NA.pops[,which(colnames(windows.NA.pops) %in% filtered.windows$x)]

eek=cbind(windows.depth.maf, windows.NA.pops$V2.x, windows.NA.pops$V2.y)
names(eek)[ncol(eek)] <- "sample.cvg"
names(eek)[ncol(eek)-1] <- "V2"
windows.NA.pops.maf=eek
window.number=ncol(windows.NA.pops.maf)-2

output=matrix(NA,nrow=0,ncol=2)
for (i in 1:window.number) {
	# make temp coverage df for each window as well as trait value
	#wd.temp=cbind(window.means.NA.temp[i], window_depth.traits$V2, window_depth.traits[pc1])
	# calculate lm between window and trait
	print(paste0("window",i , "of", window.number))
	res = lmer(windows.NA.pops.maf[,i] ~ 1 + (1|V2) + sample.cvg, data = windows.NA.pops.maf)
	among = as.data.frame(VarCorr(res))[1,4]
	within = as.data.frame(VarCorr(res))[2,4]
	qst=among/(among+2*(within))
	out.line=c(colnames(windows.NA.pops.maf[i]), qst)
	# add lm for trait and window to output, then move onto next window
	output=rbind(output,out.line)
}

write.table(as.data.frame(output), file = "qst-NA.txt", quote = F, row.names = F, col.names = F)


qst.vals.NA = read.table("qst-NA.txt")
rownames(qst.vals.NA) = qst.vals.NA[,1]
qst.vals.NA[,1] = NULL

qst.vals.NA.t = as.data.frame(t(qst.vals.NA))


neut.fst.NA=read.table("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/neut-fst/NA_all.weir.fst", header=T)

eee=as.data.frame(t(qst.vals.NA.t))

neut.fst.NA.pc5 = head(neut.fst.NA[order(abs(neut.fst.NA$WEIR_AND_COCKERHAM_FST),decreasing=T),],.05*nrow(neut.fst.NA))
neut.fst.NA.pc1 = head(neut.fst.NA[order(abs(neut.fst.NA$WEIR_AND_COCKERHAM_FST),decreasing=T),],.01*nrow(neut.fst.NA))

pc5.fst = min(neut.fst.NA.pc5$WEIR_AND_COCKERHAM_FST)
pc1.fst = min(neut.fst.NA.pc1$WEIR_AND_COCKERHAM_FST)



pheno=1
trait.associations = as.data.frame(cbind(rownames(qst.vals.NA), qst.vals.NA[,pheno]))
#### manhattan of r^2 ####
# read in chromosome data for manhattan plots
contigs = read.table("~/om62_scratch/ragweed2022/ref/chromo-length.txt", header = F)
contigs$V3 = c(1:nrow(contigs))
colnames(contigs) = c("CONTIG", "LENGTH", "CHR")
contigs$CHR = as.numeric(contigs$CHR)

# get breaks for manhattan plots
manbreaks = c()
s = 0
for (i in 1:nrow(contigs)) {
	nbp = contigs[contigs$CHR == i,]$LENGTH
	manbreaks = c(manbreaks, round(nbp/2) + s)
	s = s + nbp
}
manbreaksnames = 1:nrow(contigs)
colnames(trait.associations)=c("window", "qst")
trait.associations = trait.associations %>% separate(window, c("CONTIG", "LOC"), ":")
trait.associations$LOC = as.numeric(trait.associations$LOC)
trait.associations$qst = as.numeric(as.character(trait.associations$qst))
gwas = merge(contigs, trait.associations, by = "CONTIG")
#gwas.sig=0.05/37202
# add a cumulative BP count
gwas$LOCc = NA
s = 0
for (i in 1:length(unique(gwas$CHR))) {
	nbp = max(gwas[gwas$CHR == i,]$LENGTH)
	gwas[gwas$CHR == i, "LOCc"] = gwas[gwas$CHR == i,"LOC"] + s
	s = s + nbp
}
#top.one=round(0.005*nrow(gwas))
#n=0.05
gwas.sig.q = head(gwas[order(abs(gwas$qst),decreasing=T),],.05*nrow(gwas))
gwas.sig.1pc = head(gwas[order(abs(gwas$qst),decreasing=T),],.01*nrow(gwas))
pc5.qst = min(gwas.sig.q$qst)
pc1.qst = min(gwas.sig.1pc$qst)

gwas.sig = gwas[which(abs(gwas$qst) >= pc1.fst),]

#gwas.sig.i = gwas.sig[,which(colnames(gwas.sig) == "rsq" | colnames(gwas.sig) == "LOCc")]
assign(paste0("env.plot_", pheno), ggplot() +
	geom_point(data = gwas, aes(x = LOCc, y = abs(qst), color = as.factor(CHR)), size = 0.5, shape = 19) +
	geom_point(data = gwas.sig, aes(x = LOCc, y = abs(qst)), color = "#ffacfb", size = 0.5, shape = 19) +
	ggtitle("NA") +
	scale_color_manual(values = rep(c("#bdbdbd", "#636363"), 500)) +
	scale_x_continuous(name = "chromosome", breaks = manbreaks, labels = manbreaksnames) +
	labs(y = "QST") +
	#geom_hline(yintercept = -log10(gwas.sig), color = "black") +
	ylim(0, 1) +
	theme_classic() +
	theme(legend.position = "none"))
write.table(as.data.frame(cbind(as.character(gwas.sig[,1]), gwas.sig[,4])), file = "qst_NA_1fst.txt", quote = F, row.names = F, col.names = F)



qst.NA.hist = ggplot() +
		geom_histogram(data = eee, aes(x = V2), color = "black", fill = "white", bins = 50) +
		#geom_vline(data = eee, aes(xintercept = pc5.fst), color="#35cff9") +
		geom_vline(data = eee, aes(xintercept = pc1.fst), color="#ffacfb") +
		ggtitle("QST NA") +
		xlim(0,1) +
		labs(xlim = "Qst") +
		theme_classic()

fst.NA.hist = ggplot() +
		geom_histogram(data = neut.fst.NA, aes(x = WEIR_AND_COCKERHAM_FST), color = "black", fill = "white", bins = 50) +
		ggtitle("FST NA") +
		xlim(0,1) +
		labs(x = "Fst") +
		#geom_vline(data = neut.fst.NA, aes(xintercept = pc5.fst), color="#35cff9") +
		geom_vline(data = neut.fst.NA, aes(xintercept = pc1.fst), color="#ffacfb") +
		theme_classic()

qst.fst = plot_grid(qst.NA.hist, fst.NA.hist, nrow = 2, align = "v", labels = c("", ""), label_size = 36)


png("qst_NA_final.png", width = 1000, height= 400)
plot_grid(env.plot_1, qst.fst, ncol = 2, align = "h", rel_widths = c(3, 1), labels = c("", ""), label_size = 36)
dev.off()









#############################################################################################################################################################################
########################################################################    TRAITS    #######################################################################################
#############################################################################################################################################################################




windows.depth=read.table("window_depth_2SD5_modern_Q30.txt", header = T)

bla = str_replace_all(colnames(windows.depth), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")

colnames(windows.depth) = bla
windows.depth$maf=NULL


#rownames(windows.depth) = paste0(windows.depth$ID, ":", windows.depth$window)
#windows.depth[,1:2] = NULL
#filtered.windows=read.table("modern_Q30.txt_2SD5_windows.list")
#
#windows.depth.maf=windows.depth[which(rownames(windows.depth) %in% filtered.windows$x),]

windows.depth.maf.t = as.data.frame(t(windows.depth))

windows.list = colnames(windows.depth.maf.t)
traits.df = read.table("../traits.list", header = F)
traits = as.character(traits.df$V1)

pheno=traits[1]
assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
trait.info = get(pheno)
rownames(trait.info) = trait.info$V1
colnames(trait.info) = c("sample", "sample2", paste0(pheno))

window_depth.traits = merge(windows.depth.maf.t, trait.info[3], by = 0)
rownames(window_depth.traits) <- window_depth.traits[,1]
window_depth.traits[,1] <- NULL

for (pheno in traits[2:length(traits)]) {
		assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
		trait.info = get(pheno)
		rownames(trait.info) = trait.info$V1
		colnames(trait.info) = c("sample", "sample2", paste0(pheno))
		
		window_depth.traits = merge(window_depth.traits, trait.info[3], by = 0)
		rownames(window_depth.traits) <- window_depth.traits[,1]
		window_depth.traits[,1] <- NULL
}

neutral_cov=read.table("/home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30/figs/genotyping/modelling/neutral_gwas.cov", header=F)
neutral_cov_samples=read.table("/home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30/figs/genotyping/modelling/gwas_samples.list", header=F)
e = eigen(neutral_cov)
neutral_pcs <- e$vectors
neutral_pcs = cbind(neutral_cov_samples, neutral_pcs[,1:2])
rownames(neutral_pcs) <- neutral_pcs[,1]
neutral_pcs[,1] <- NULL







window_depth.traits = merge(window_depth.traits, neutral_pcs[1], by=0)
rownames(window_depth.traits) = window_depth.traits[,1]
window_depth.traits[,1] <- NULL

pc1=ncol(window_depth.traits)
tr.len=length(traits)+1
end.ind=ncol(window_depth.traits)-tr.len
windows.list=colnames(window_depth.traits[1:end.ind])

alloutput=matrix(NA,nrow=0,ncol=4)


# iterate through traits list
for (j in 1:length(traits)){
	print(paste0("beginning phenotype ", traits[j]))
	output=matrix(NA,nrow=0,ncol=3)
	outputsig=matrix(NA,nrow=0,ncol=4)
	for (i in 1:(length(windows.list)-1)) {
		# make temp coverage df for each window as well as trait value
		wd.temp=cbind(window_depth.traits[i], window_depth.traits[(length(windows.list)+j)], window_depth.traits[pc1])
		# calculate lm between window and trait
		res = lm(as.numeric(wd.temp[,2]) ~ as.numeric(wd.temp[,1]) + as.numeric(wd.temp[,3]))
		pval = summary(res)$coefficients[2,4]
		rsq = summary(res)$adj.r.squared
		out.line=c(colnames(wd.temp[1]),pval,rsq)
		# add lm for trait and window to output, then move onto next window
		output=rbind(output,out.line)
	}
	write.table(as.data.frame(output), file = paste0(traits[j], "_window_depth_pheno_nrml_Q30.txt"), quote = F, row.names = F, col.names = F)
	write.table(as.data.frame(outputsig), file = paste0(traits[j], "_windows_traits_sig.txt"), quote = F, row.names = F, col.names = F)
}

write.table(as.data.frame(alloutput), file = "windows_all_traits_sig.txt", quote = F, row.names = F, col.names = F)





### plotting 


for (pheno in traits) {
	assign(pheno, read.table(paste0(pheno, "_window_depth_pheno_nrml_Q30.txt")))
	trait.associations = get(pheno)
	#### manhattan of r^2 ####
	# read in chromosome data for manhattan plots
	contigs = read.table("~/om62_scratch/ragweed2022/ref/chromo-length.txt", header = F)
	contigs$V3 = c(1:nrow(contigs))
	colnames(contigs) = c("CONTIG", "LENGTH", "CHR")
	contigs$CHR = as.numeric(contigs$CHR)
	
	# get breaks for manhattan plots
	manbreaks = c()
	s = 0
	for (i in 1:nrow(contigs)) {
		nbp = contigs[contigs$CHR == i,]$LENGTH
		manbreaks = c(manbreaks, round(nbp/2) + s)
		s = s + nbp
	}
	manbreaksnames = 1:nrow(contigs)
	colnames(trait.associations)=c("window", "pval", "rsq")
	trait.associations = trait.associations %>% separate(window, c("CONTIG", "LOC"), ":")
	trait.associations$LOC = as.numeric(trait.associations$LOC)
	
	gwas = merge(contigs, trait.associations, by = "CONTIG")
	gwas.sig.line=0.05/end.ind


	# add a cumulative BP count
	gwas$LOCc = NA
	s = 0
	for (i in 1:length(unique(gwas$CHR))) {
		nbp = max(gwas[gwas$CHR == i,]$LENGTH)
		gwas[gwas$CHR == i, "LOCc"] = gwas[gwas$CHR == i,"LOC"] + s
		s = s + nbp
	}
	gwas.sig = as.data.frame(gwas[which(gwas$pval <= gwas.sig.line),])
	assign(paste0("trait.plot_", pheno), ggplot() +
		geom_point(data = gwas, aes(x = LOCc, y = -log10(pval), color = as.factor(CHR)), size = 0.5, shape = 19) +
		geom_point(data = gwas.sig, aes(x = LOCc, y = -log10(pval)), color = "#f37ded", size = 1, shape = 19) +
		geom_hline(yintercept=-log10(gwas.sig.line)) +
		ggtitle(paste0(pheno)) +
		scale_color_manual(values = rep(c("#bdbdbd", "#636363"), 500)) +
		scale_x_continuous(name = "chromosome", breaks = manbreaks, labels = manbreaksnames) +
		labs(y = expression(-log[10]*italic(p))) +
		ylim(0, 10) +
		theme_classic() +
		theme(legend.position = "none"))
	write.table(as.data.frame(cbind(as.character(gwas.sig[,1]), gwas.sig[,4], gwas.sig[,5], gwas.sig[,6])), file = paste0(pheno, "_outliers_2pt5SD_Q30.txt"), quote = F, row.names = F, col.names = F)
}


png("window_depth_lm_associations_update_nrml_Q30.png", width = 2000, height= 4000)
grid.arrange(trait.plot_T_die_day, trait.plot_T_fem_day, trait.plot_T_fit_all, trait.plot_T_fl_end, trait.plot_T_fl_start, trait.plot_T_InflxPt_Day, trait.plot_T_LinearGrowth_Days, trait.plot_T_Longest_leaf, trait.plot_T_male_day, trait.plot_T_male_weight, trait.plot_T_maxHeight_cm, trait.plot_T_MaxSlope_cm_Day, trait.plot_T_momseed_mg, trait.plot_T_pol_day, trait.plot_T_rac_die_day, trait.plot_T_Raceme_length_longest, trait.plot_T_Racemes, trait.plot_T_repal, trait.plot_T_rootshoot_rat, trait.plot_T_root_weights, trait.plot_T_seed_20, trait.plot_T_seeds_day, trait.plot_T_seeds_tot, trait.plot_T_seed_weight_tot, trait.plot_T_sexmismatch, trait.plot_T_sexrat_weight, trait.plot_T_shoot_biom, trait.plot_T_Stem_width_mum, trait.plot_T_tot_biom, ncol = 2)
dev.off()

png("window_depth_lm_associations_update_nrml_Q30_sig.png", width = 2000, height= 1000)
grid.arrange(trait.plot_T_die_day, trait.plot_T_fem_day, trait.plot_T_fl_end, trait.plot_T_fl_start, trait.plot_T_InflxPt_Day, trait.plot_T_LinearGrowth_Days, trait.plot_T_male_day, trait.plot_T_maxHeight_cm, trait.plot_T_pol_day, trait.plot_T_rac_die_day, trait.plot_T_repal, trait.plot_T_root_weights, trait.plot_T_seeds_day, trait.plot_T_sexmismatch, trait.plot_T_sexrat_weight, trait.plot_T_shoot_biom, trait.plot_T_Stem_width_mum, trait.plot_T_tot_biom, ncol = 3)
dev.off()








#############################################################################################################################################################################
################################################################    ENVIRONMENTAL ASSOCIATIONS    ###########################################################################
#############################################################################################################################################################################



windows.depth=read.table("window_depth_2SD5_modern_Q30.txt", header = T)

bla = str_replace_all(colnames(windows.depth), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")

colnames(windows.depth) = bla

window_depth.filt.t.ii = as.data.frame(t(windows.depth))



########################## EU ##########################



EU.env=read.table("EU-baypass-samples-pops-env.txt")
rownames(EU.env) <- EU.env[,1]
EU.env[,1] <- NULL

window_depth.filt.t.ii.EU=merge(window_depth.filt.t.ii, EU.env, by=0)
rownames(window_depth.filt.t.ii.EU)=window_depth.filt.t.ii.EU[,1]
window_depth.filt.t.ii.EU[,1] <- NULL
window_depth.filt.t.ii.EU$V2 = as.character(window_depth.filt.t.ii.EU$V2)
maf.wind=ncol(window_depth.filt.t.ii.EU)-ncol(EU.env)
window.names.EU = colnames(window_depth.filt.t.ii.EU[,1:(maf.wind)])

for (i in 1:maf.wind){window_depth.filt.t.ii.EU[,i]=as.numeric(as.character(window_depth.filt.t.ii.EU[,i]))}

window.summary.EU = window_depth.filt.t.ii.EU %>% group_by(V2) %>% filter(n() > 2) %>% summarise_at(window.names.EU, mean, na.rm = TRUE)
window.means.EU=as.data.frame(window.summary.EU)
rownames(window.means.EU)=window.means.EU[,1]
window.means.EU[,1] <- NULL

#window.means.EU.t=as.data.frame(window.summary.EU)

write.table(window.means.EU, file = "window_EU_pops_mean_Q30.txt", quote = F, row.names = T)


env.var.names.pop=c("pop", "lat", "lon", "amt", "mdr", "it", "ts", "mtwm", "mtcm", "tar", "mtweq", "mtdq", "mtwaq", "mtcq", "ap", "pwm", "pdm", "ps", "pweq", "pdq", "pwaq", "pcq")
env.var.names.pc=c("lat", "lon", "amt", "mdr", "it", "ts", "mtwm", "mtcm", "tar", "mtweq", "mtdq", "mtwaq", "mtcq", "ap", "pwm", "pdm", "ps", "pweq", "pdq", "pwaq", "pcq")
env.var.names=c("lat", "lon", "amt", "mdr", "it", "ts", "mtwm", "mtcm", "tar", "mtweq", "mtdq", "mtwaq", "mtcq", "ap", "pwm", "pdm", "ps", "pweq", "pdq", "pwaq", "pcq")

colnames(EU.env)=env.var.names.pop

window.means.EU=read.table("window_EU_pops_mean_Q30.txt", header=T)

env.summary.EU = EU.env %>% group_by(pop) %>% filter(n() > 2) %>% summarise_at(env.var.names.pc, mean, na.rm = TRUE)
env.means.EU=as.data.frame(env.summary.EU)
rownames(env.means.EU)=env.means.EU[,1]
env.means.EU[,1]=NULL

window.env.means.EU=merge(window.means.EU, env.means.EU, by=0)
rownames(window.env.means.EU)=window.env.means.EU[,1]
window.env.means.EU[,1] <- NULL

maf.wind = ncol(window.means.EU)
window.names.EU = colnames(window.env.means.EU[,1:maf.wind])

col.tot=ncol(window.env.means.EU)

# iterate through traits list
for (j in 1:length(env.var.names)){
	print(paste0("beginning env variable... ", env.var.names[j]))
	output=matrix(NA,nrow=0,ncol=3)
	for (i in 1:length(window.names.EU)) {
		# make temp coverage df for each window as well as trait value
		wd.temp=cbind(window.env.means.EU[i], window.env.means.EU[(length(window.names.EU)+j)])
		# calculate lm between window and trait
		#res = lm(as.numeric(wd.temp[,1]) ~ as.numeric(wd.temp[,2] + as.numeric(wd.temp[,3])))
		res = cor.test(as.numeric(as.character(wd.temp[,1])), as.numeric(as.character(wd.temp[,2])), method="spearman", exact=FALSE)
		#pval = summary(res)$coefficients[2,4]
		pval = res$p.value[1]
		#rsq = summary(res)$adj.r.squared
		tau = res$estimate[1]
		out.line=c(colnames(wd.temp[1]),pval,tau)
		# add lm for trait and window to output, then move onto next window
		output=rbind(output,out.line)
	}
	#write.table(as.data.frame(output), file = paste0(env.var.names[j], "-window-depth-EU-env-lm.txt"), quote = F, row.names = F, col.names = F)
	write.table(as.data.frame(output), file = paste0(env.var.names[j], "_EU_env_pops_spearmans_Q30.txt"), quote = F, row.names = F, col.names = F)
}


for (pheno in env.var.names) {
	assign(pheno, read.table(paste0(pheno, "_EU_env_pops_spearmans_Q30.txt")))
	trait.associations = get(pheno)
	#### manhattan of r^2 ####
	# read in chromosome data for manhattan plots
	contigs = read.table("~/om62_scratch/ragweed2022/ref/chromo-length.txt", header = F)
	contigs$V3 = c(1:nrow(contigs))
	colnames(contigs) = c("CONTIG", "LENGTH", "CHR")
	contigs$CHR = as.numeric(contigs$CHR)
	
	# get breaks for manhattan plots
	manbreaks = c()
	s = 0
	for (i in 1:nrow(contigs)) {
		nbp = contigs[contigs$CHR == i,]$LENGTH
		manbreaks = c(manbreaks, round(nbp/2) + s)
		s = s + nbp
	}
	manbreaksnames = 1:nrow(contigs)
	colnames(trait.associations)=c("window", "pval", "rsq")
	trait.associations = trait.associations %>% separate(window, c("CONTIG", "LOC"), "\\.")
	trait.associations$LOC = as.numeric(trait.associations$LOC)
	
	gwas = merge(contigs, trait.associations, by = "CONTIG")
	#gwas.sig=0.05/37202

	# add a cumulative BP count
	gwas$LOCc = NA
	s = 0
	for (i in 1:length(unique(gwas$CHR))) {
		nbp = max(gwas[gwas$CHR == i,]$LENGTH)
		gwas[gwas$CHR == i, "LOCc"] = gwas[gwas$CHR == i,"LOC"] + s
		s = s + nbp
	}

	#top.one=round(0.005*nrow(gwas))
	#n=0.05
	gwas.sig = head(gwas[order(abs(gwas$rsq),decreasing=T),],.01*nrow(gwas))
	#gwas.sig.i = gwas.sig[,which(colnames(gwas.sig) == "rsq" | colnames(gwas.sig) == "LOCc")]
	assign(paste0("env.plot_", pheno), ggplot() +
		geom_point(data = gwas, aes(x = LOCc, y = -log10(pval), color = as.factor(CHR)), size = 0.5, shape = 19) +
		geom_point(data = gwas.sig, aes(x = LOCc, y = -log10(pval)), color = "#f37ded", size = 1, shape = 19) +
		ggtitle(paste0(pheno)) +
		scale_color_manual(values = rep(c("#bdbdbd", "#636363"), 500)) +
		scale_x_continuous(name = "chromosome", breaks = manbreaks, labels = manbreaksnames) +
		labs(y = "-log10 pval") +
		#geom_hline(yintercept = -log10(gwas.sig), color = "black") +
		ylim(0, 10) +
		theme_classic() +
		theme(legend.position = "none"))
	write.table(as.data.frame(cbind(as.character(gwas.sig[,1]), gwas.sig[,4])), file = paste0(pheno, "_EU_5pc_sp_Q30.txt"), quote = F, row.names = F, col.names = F)

}

png("EU_depth_env_Q30.png", width = 1500, height= 3000)
grid.arrange(env.plot_lat, env.plot_lon, env.plot_amt, env.plot_mdr, env.plot_it, env.plot_ts, env.plot_mtwm, env.plot_mtcm, env.plot_tar, env.plot_mtweq, env.plot_mtdq, env.plot_mtwaq, env.plot_mtcq, env.plot_ap, env.plot_pwm, env.plot_pdm, env.plot_ps, env.plot_pweq, env.plot_pdq, env.plot_pwaq, env.plot_pcq, ncol = 1)
dev.off()

png("EU_depth_env_Q30_r7.png", width = 1500, height= 2000)
grid.arrange(env.plot_amt, env.plot_mdr, env.plot_mtweq, env.plot_mtdq,env.plot_ap, env.plot_ps, ncol = 1)
dev.off()



########################## NA ##########################



NA.env=read.table("NA-baypass-samples-pops-env.txt")
rownames(NA.env) <- NA.env[,1]
NA.env[,1] <- NULL

window_depth.filt.t.ii.NA=merge(window_depth.filt.t.ii, NA.env, by=0)
rownames(window_depth.filt.t.ii.NA)=window_depth.filt.t.ii.NA[,1]
window_depth.filt.t.ii.NA[,1] <- NULL
window_depth.filt.t.ii.NA$V2 = as.character(window_depth.filt.t.ii.NA$V2)
maf.wind=ncol(window_depth.filt.t.ii.NA)-ncol(NA.env)
window.names.NA = colnames(window_depth.filt.t.ii.NA[,1:(maf.wind)])

for (i in 1:maf.wind){window_depth.filt.t.ii.NA[,i]=as.numeric(as.character(window_depth.filt.t.ii.NA[,i]))}

window.summary.NA = window_depth.filt.t.ii.NA %>% group_by(V2) %>% filter(n() > 2) %>% summarise_at(window.names.NA, mean, na.rm = TRUE)
window.means.NA=as.data.frame(window.summary.NA)
rownames(window.means.NA)=window.means.NA[,1]
window.means.NA[,1] <- NULL

#window.means.NA.t=as.data.frame(window.summary.NA)

write.table(window.means.NA, file = "window_NA_pops_mean_Q30.txt", quote = F, row.names = T)


env.var.names.pop=c("pop", "lat", "lon", "amt", "mdr", "it", "ts", "mtwm", "mtcm", "tar", "mtweq", "mtdq", "mtwaq", "mtcq", "ap", "pwm", "pdm", "ps", "pweq", "pdq", "pwaq", "pcq")
env.var.names.pc=c("lat", "lon", "amt", "mdr", "it", "ts", "mtwm", "mtcm", "tar", "mtweq", "mtdq", "mtwaq", "mtcq", "ap", "pwm", "pdm", "ps", "pweq", "pdq", "pwaq", "pcq")
env.var.names=c("lat", "lon", "amt", "mdr", "it", "ts", "mtwm", "mtcm", "tar", "mtweq", "mtdq", "mtwaq", "mtcq", "ap", "pwm", "pdm", "ps", "pweq", "pdq", "pwaq", "pcq")

colnames(NA.env)=env.var.names.pop

window.means.NA=read.table("window_NA_pops_mean_Q30.txt", header=T)

env.summary.NA = NA.env %>% group_by(pop) %>% filter(n() > 2) %>% summarise_at(env.var.names.pc, mean, na.rm = TRUE)
env.means.NA=as.data.frame(env.summary.NA)
rownames(env.means.NA)=env.means.NA[,1]
env.means.NA[,1]=NULL

window.env.means.NA=merge(window.means.NA, env.means.NA, by=0)
rownames(window.env.means.NA)=window.env.means.NA[,1]
window.env.means.NA[,1] <- NULL

maf.wind = ncol(window.means.NA)
window.names.NA = colnames(window.env.means.NA[,1:maf.wind])

col.tot=ncol(window.env.means.NA)

# iterate through traits list
for (j in 1:length(env.var.names)){
	print(paste0("beginning env variable... ", env.var.names[j]))
	output=matrix(NA,nrow=0,ncol=3)
	for (i in 1:length(window.names.NA)) {
		# make temp coverage df for each window as well as trait value
		wd.temp=cbind(window.env.means.NA[i], window.env.means.NA[(length(window.names.NA)+j)])
		# calculate lm between window and trait
		#res = lm(as.numeric(wd.temp[,1]) ~ as.numeric(wd.temp[,2] + as.numeric(wd.temp[,3])))
		res = cor.test(as.numeric(as.character(wd.temp[,1])), as.numeric(as.character(wd.temp[,2])), method="spearman", exact=FALSE)
		#pval = summary(res)$coefficients[2,4]
		pval = res$p.value[1]
		#rsq = summary(res)$adj.r.squared
		tau = res$estimate[1]
		out.line=c(colnames(wd.temp[1]),pval,tau)
		# add lm for trait and window to output, then move onto next window
		output=rbind(output,out.line)
	}
	#write.table(as.data.frame(output), file = paste0(env.var.names[j], "-window-depth-NA-env-lm.txt"), quote = F, row.names = F, col.names = F)
	write.table(as.data.frame(output), file = paste0(env.var.names[j], "_NA_env_pops_spearmans_Q30.txt"), quote = F, row.names = F, col.names = F)
}


for (pheno in env.var.names) {
	assign(pheno, read.table(paste0(pheno, "_NA_env_pops_spearmans_Q30.txt")))
	trait.associations = get(pheno)
	#### manhattan of r^2 ####
	# read in chromosome data for manhattan plots
	contigs = read.table("~/om62_scratch/ragweed2022/ref/chromo-length.txt", header = F)
	contigs$V3 = c(1:nrow(contigs))
	colnames(contigs) = c("CONTIG", "LENGTH", "CHR")
	contigs$CHR = as.numeric(contigs$CHR)
	
	# get breaks for manhattan plots
	manbreaks = c()
	s = 0
	for (i in 1:nrow(contigs)) {
		nbp = contigs[contigs$CHR == i,]$LENGTH
		manbreaks = c(manbreaks, round(nbp/2) + s)
		s = s + nbp
	}
	manbreaksnames = 1:nrow(contigs)
	colnames(trait.associations)=c("window", "pval", "rsq")
	trait.associations = trait.associations %>% separate(window, c("CONTIG", "LOC"), "\\.")
	trait.associations$LOC = as.numeric(trait.associations$LOC)
	
	gwas = merge(contigs, trait.associations, by = "CONTIG")
	#gwas.sig=0.05/37202

	# add a cumulative BP count
	gwas$LOCc = NA
	s = 0
	for (i in 1:length(unique(gwas$CHR))) {
		nbp = max(gwas[gwas$CHR == i,]$LENGTH)
		gwas[gwas$CHR == i, "LOCc"] = gwas[gwas$CHR == i,"LOC"] + s
		s = s + nbp
	}

	#top.one=round(0.005*nrow(gwas))
	#n=0.05
	gwas.sig = head(gwas[order(abs(gwas$rsq),decreasing=T),],.01*nrow(gwas))
	#gwas.sig.i = gwas.sig[,which(colnames(gwas.sig) == "rsq" | colnames(gwas.sig) == "LOCc")]
	assign(paste0("env.plot_", pheno), ggplot() +
		geom_point(data = gwas, aes(x = LOCc, y = -log10(pval), color = as.factor(CHR)), size = 0.5, shape = 19) +
		geom_point(data = gwas.sig, aes(x = LOCc, y = -log10(pval)), color = "#f37ded", size = 1, shape = 19) +
		ggtitle(paste0(pheno)) +
		scale_color_manual(values = rep(c("#bdbdbd", "#636363"), 500)) +
		scale_x_continuous(name = "chromosome", breaks = manbreaks, labels = manbreaksnames) +
		labs(y = "-log10 pval") +
		#geom_hline(yintercept = -log10(gwas.sig), color = "black") +
		ylim(0, 10) +
		theme_classic() +
		theme(legend.position = "none"))
	write.table(as.data.frame(cbind(as.character(gwas.sig[,1]), gwas.sig[,4])), file = paste0(pheno, "_NA_5pc_sp_Q30.txt"), quote = F, row.names = F, col.names = F)

}

png("NA_depth_env_Q30.png", width = 1500, height= 3000)
grid.arrange(env.plot_lat, env.plot_lon, env.plot_amt, env.plot_mdr, env.plot_it, env.plot_ts, env.plot_mtwm, env.plot_mtcm, env.plot_tar, env.plot_mtweq, env.plot_mtdq, env.plot_mtwaq, env.plot_mtcq, env.plot_ap, env.plot_pwm, env.plot_pdm, env.plot_ps, env.plot_pweq, env.plot_pdq, env.plot_pwaq, env.plot_pcq, ncol = 1)
dev.off()

png("NA_depth_env_Q30_r7.png", width = 1500, height= 2000)
grid.arrange(env.plot_amt, env.plot_mdr, env.plot_mtweq, env.plot_mtdq,env.plot_ap, env.plot_ps, ncol = 1)
dev.off()



#############################################################################################################################################################################
##################################################################       OVERLAPS IN BASH       #############################################################################
#############################################################################################################################################################################




cd ~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/cnv_results_final/overlaps_final



while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < maxheight_windows.txt


















while read line1 line2; do mv ../${line1}_EU_1pc_sp_Q30.txt ${line2}_EU_1pc_sp_Q30.txt; done < env-info.list
while read line1 line2; do mv ../${line1}_NA_1pc_sp_Q30.txt ${line2}_NA_1pc_sp_Q30.txt; done < env-info.list



for i in amt mdr mtweq mtdq ap ps; do grep -f ${i}_EU_5pc_sp_Q30.txt qst_EU_1fst.txt; done | sort -u > EU_qst-env.txt
for i in amt mdr mtweq mtdq ap ps; do grep -f ${i}_NA_5pc_sp_Q30.txt qst_NA_1fst.txt; done | sort -u > NA_qst-env.txt


windows_all_traits_sig.txt
T_fl_start_outliers_2pt5SD_Q30.txt

while read line; do cat ${line}_outliers_2pt5SD_Q30.txt | awk '{print $1, $2}'; done < traits.list | sort -u > all_traits_sig.list

for i in amt mdr mtweq mtdq ap ps; do grep -f ${i}_NA_5pc_sp_Q30.txt all_traits_sig.list; echo $i; done | sort -u 
for i in amt mdr mtweq mtdq ap ps; do grep -f ${i}_EU_5pc_sp_Q30.txt all_traits_sig.list; echo $i; done | sort -u 

grep -f qst_NA_1fst.txt all_traits_sig.list | sort -u 
grep -f qst_EU_1fst.txt all_traits_sig.list | sort -u 

cat all_traits_sig.list | awk '{print $1, $2-9999, $2}' > all_traits_sig_full.list

while read line1 line2 line3;  do cat /home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30/figs/genotyping/modelling/genes_all_traits_sig.locs| awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$1 && start >= $2 && start <= $3 || chr==$1 && end >= $2 && end <= $3 || chr==$1 && start >= $2 && end <= $3 || chr==$1 && start <= $2 && end >= $3)print $0}';  done < all_traits_sig_full.list








cat NA_qst_traits.list  | awk '{print $1, $2-19999, $2}' > NA_qst_traits_windows.bed
cat EU_qst_traits.list  | awk '{print $1, $2-9999, $2}' > EU_qst_traits_windows.bed

while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < NA_qst_traits_windows.bed > NA_qst_traits_windows_genes.txt
while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < EU_qst_traits_windows.bed > EU_qst_traits_windows_genes.txt





############ QST only


cat qst_EU_1fst.txt | awk '{print $1, $2-9999, $2}' > qst_EU_full.txt
cat qst_NA_1fst.txt | awk '{print $1, $2-9999, $2}' > qst_NA_full.txt

grep -f qst_EU_full.txt qst_NA_full.txt > qst_OL_full.txt
# 70 out of 270


while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < qst_EU_full.txt > qst_EU_genes.txt
# 85
while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < qst_NA_full.txt > qst_NA_genes.txt
# 53
while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < qst_OL_full.txt > qst_OL_genes.txt

cat qst_EU_genes.txt | awk '{print $3, $4, $5}' > qst_EU_genes.bed
cat qst_NA_genes.txt | awk '{print $3, $4, $5}' > qst_NA_genes.bed
cat qst_OL_genes.txt | awk '{print $3, $4, $5}' > qst_OL_genes.bed





# create GO input file
cat qst_EU_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > qst_EU_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < qst_EU_genes_minus1.bed | sort -u > qst_EU.topGO


cat qst_NA_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > qst_NA_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < qst_NA_genes_minus1.bed | sort -u > qst_NA.topGO


cat qst_OL_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > qst_OL_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < qst_OL_genes_minus1.bed | sort -u > qst_OL.topGO


library(topGO)
setwd("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/cnv_results_final/")


### EU

geneID2GO = readMappings(file = "allannot.topGO")
geneUniverse = names(geneID2GO)


genesOfInterest = read.table("qst_EU.topGO")
genesOfInterest = as.character(genesOfInterest$V1)
geneList = factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) = geneUniverse
# ontology = 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component)
myGOdata = new("topGOdata",
               description = "My project",
               ontology = "BP",
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)
resultFisher = runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
allResSNP = GenTable(myGOdata,
                  classicFisher = resultFisher,
                  orderBy = "resultFisher",
                  ranksOf = "classicFisher",
                  topNodes = 100)
allResSNP
write.table(x = allResSNP, file = "qst_EU.GO.txt", sep = "\t", quote = F)



### NA

geneID2GO = readMappings(file = "allannot.topGO")
geneUniverse = names(geneID2GO)


genesOfInterest = read.table("qst_NA.topGO")
genesOfInterest = as.character(genesOfInterest$V1)
geneList = factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) = geneUniverse
# ontology = 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component)
myGOdata = new("topGOdata",
               description = "My project",
               ontology = "BP",
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)
resultFisher = runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
allResSNP = GenTable(myGOdata,
                  classicFisher = resultFisher,
                  orderBy = "resultFisher",
                  ranksOf = "classicFisher",
                  topNodes = 100)
allResSNP
write.table(x = allResSNP, file = "qst_NA.GO.txt", sep = "\t", quote = F)


### OL

geneID2GO = readMappings(file = "allannot.topGO")
geneUniverse = names(geneID2GO)


genesOfInterest = read.table("qst_OL.topGO")
genesOfInterest = as.character(genesOfInterest$V1)
geneList = factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) = geneUniverse
# ontology = 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component)
myGOdata = new("topGOdata",
               description = "My project",
               ontology = "BP",
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)
resultFisher = runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
allResSNP = GenTable(myGOdata,
                  classicFisher = resultFisher,
                  orderBy = "resultFisher",
                  ranksOf = "classicFisher",
                  topNodes = 100)
allResSNP
write.table(x = allResSNP, file = "qst_OL.GO.txt", sep = "\t", quote = F)













while read line1 line2; do grep -f ../${line2}_outliers_2pt5SD_Q30.txt ../qst_EU_1fst.txt > ${line2}_env_qst_overlap_EU.txt; done < env-info.list
while read line1 line2; do cat ${line2}_env_qst_overlap_EU.txt; done < env-info.list | sort -u > all_env_qst_overlap_EU.txt
cat all_env_qst_overlap_EU.txt | wc -l
# 91

while read line1 line2; do grep -f ${line2}_NA_1pc_sp_Q30.txt ../qst_NA_1pc_outliers_Q30.txt > ${line2}_env_qst_overlap_NA.txt; done < env-info.list
while read line1 line2; do cat ${line2}_env_qst_overlap_NA.txt; done < env-info.list | sort -u > all_env_qst_overlap_NA.txt
cat all_env_qst_overlap_NA.txt | wc -l
# 178


cat all_env_qst_overlap_EU.txt| awk '{print $1, $2-9999, $2}' > all_env_qst_overlap_EU_full.txt
cat all_env_qst_overlap_NA.txt| awk '{print $1, $2-9999, $2}' > all_env_qst_overlap_NA_full.txt


while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < all_env_qst_overlap_EU_full.txt > all_env_qst_overlap_EU_genes.txt
# 37
while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < all_env_qst_overlap_NA_full.txt > all_env_qst_overlap_NA_genes.txt
# 25

cat all_env_qst_overlap_EU_genes.txt | awk '{print $3, $4, $5}' > all_env_qst_overlap_EU_genes.bed
cat all_env_qst_overlap_NA_genes.txt | awk '{print $3, $4, $5}' > all_env_qst_overlap_NA_genes.bed


# create GO input file
cat all_env_qst_overlap_EU_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > all_env_qst_overlap_EU_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < all_env_qst_overlap_EU_genes_minus1.bed | sort -u > all_env_qst_overlap_EU.topGO


cat all_env_qst_overlap_NA_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > all_env_qst_overlap_NA_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < all_env_qst_overlap_NA_genes_minus1.bed | sort -u > all_env_qst_overlap_NA.topGO


#grep -f all_env_qst_overlap_EU.txt qst_EU_histmod_outliers_Q30.txt | wc -l
#grep -f all_env_qst_overlap_EU.txt qst_EU_histmod_outliers_Q30.txt > hist_env_EU_Q30.txt
#
#grep -f all_env_qst_overlap_NA.txt qst_NA_histmod_outliers_Q30.txt | wc -l
#grep -f all_env_qst_overlap_NA.txt qst_NA_histmod_outliers_Q30.txt > hist_env_NA_Q30.txt
#
#
#cat hist_env_EU_Q30.txt| awk '{print $0, $2+9999}' > hist_env_EU_Q30_full.txt
#cat hist_env_NA_Q30.txt| awk '{print $0, $2+9999}' > hist_env_NA_Q30_full.txt
#
#while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < hist_env_EU_Q30_full.txt > hist_env_EU_Q30_full_genes.txt
#
#while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < hist_env_NA_Q30_full.txt > hist_env_NA_Q30_full_genes.txt




cp ../overlaps_q30/allannot.topGO .




#############################################################################################################################################################################
##################################################################       QST CANDIDATE GO       #############################################################################
#############################################################################################################################################################################


############ QST only


cat ../qst_EU_1pc_outliers_Q30.txt | awk '{print $1, $2-9999, $2}' > qst_EU_full.txt
cat ../qst_NA_1pc_outliers_Q30.txt | awk '{print $1, $2-9999, $2}' > qst_NA_full.txt

grep -f qst_EU_full.txt qst_NA_full.txt | wc -l
# 70 out of 270


while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < qst_EU_full.txt > qst_EU_genes.txt
# 85
while read line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}';  done < qst_NA_full.txt > qst_NA_genes.txt
# 53

cat qst_EU_genes.txt | awk '{print $3, $4, $5}' > qst_EU_genes.bed
cat qst_NA_genes.txt | awk '{print $3, $4, $5}' > qst_NA_genes.bed





# create GO input file
cat qst_EU_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > qst_EU_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < qst_EU_genes_minus1.bed | sort -u > qst_EU.topGO


cat qst_NA_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > qst_NA_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < qst_NA_genes_minus1.bed | sort -u > qst_NA.topGO




library(topGO)
setwd("/home/jwilson/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/overlaps_Q30")



### EU

geneID2GO = readMappings(file = "allannot.topGO")
geneUniverse = names(geneID2GO)


genesOfInterest = read.table("qst_EU.topGO")
genesOfInterest = as.character(genesOfInterest$V1)
geneList = factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) = geneUniverse
# ontology = 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component)
myGOdata = new("topGOdata",
               description = "My project",
               ontology = "BP",
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)
resultFisher = runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
allResSNP = GenTable(myGOdata,
                  classicFisher = resultFisher,
                  orderBy = "resultFisher",
                  ranksOf = "classicFisher",
                  topNodes = 100)
allResSNP
write.table(x = allResSNP, file = "qst_EU.GO.txt", sep = "\t", quote = F)



### NA

geneID2GO = readMappings(file = "allannot.topGO")
geneUniverse = names(geneID2GO)


genesOfInterest = read.table("qst_NA.topGO")
genesOfInterest = as.character(genesOfInterest$V1)
geneList = factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) = geneUniverse
# ontology = 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component)
myGOdata = new("topGOdata",
               description = "My project",
               ontology = "BP",
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)
resultFisher = runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
allResSNP = GenTable(myGOdata,
                  classicFisher = resultFisher,
                  orderBy = "resultFisher",
                  ranksOf = "classicFisher",
                  topNodes = 100)
allResSNP
write.table(x = allResSNP, file = "qst_NA.GO.txt", sep = "\t", quote = F)









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
cnv.name[["24890001"]] = "cnv-chr17a"
cnv.name[["21650001"]] = "cnv-chr18a"
cnv.name[["31240001"]] = "cnv-chr18b"
cnv.name[["12550001"]] = "cnv-chr5a"
cnv.name[["30730001"]] = "cnv-chr17b"
cnv.name[["10700001"]] = "cnv-chr8a"
cnv.name[["19190001"]] = "cnv-chr10a"
cnv.name[["21790001"]] = "cnv-chr13a"
cnv.name[["20190001"]] = "cnv-chr17c"
cnv.name[["29680001"]] = "cnv-chr6a"
cnv.name[["20640001"]] = "cnv-chr3a"











chr = "h1s2"
chr.start = 45530001
chr.end = 46160000


cnv.gt.mod = read.table("h1s2:45530001-46160000_mod_gt.txt", header = T)
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




cnv.gt.hist = read.table("h1s2:45530001-46160000_hist_gt.txt", header = T)
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

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + pc1, family=binomial)
Anova(glm, type = 3)


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
  labs(x = "lat",
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
  labs(x = "lat",
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



cnv2a_model_nonsig = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)







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


emtrends(glm, pairwise ~ range, var="lat", at=list(range="na", year=c(quantile(., 0.25), quantile(., 0.5), quantile(., 0.75),))
means1_an_4a = emtrends(glm, "year", var = "lat", adjust="fdr") #, at=list(variety="A", temp=c(20,40)))

emtrends(scaf27b.glm, pairwise ~ year, var="latitude", at=list(range="North_America", year=c(early_hist,med_hist,late_hist, med_modern)))















means1_an_4a = ggemmeans(glm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
at=list(variety="A", temp=c(20,40)












regrid(as.data.frame(means1_an_4a))


pairs(regrid(means1_an_4a), simple = "each", combine = TRUE, adjust="fdr")







means1_an_4a = (emtrends(glm, "year [1830:2019 by=40]", var = "lat", adjust="fdr"))

summary(means1_an_4b)
means2_an_4b = emtrends(glm, pairwise ~ range, var="lat")

#emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")

means2_an_4a = emmeans(glm
               , type = "pairs"
               , specs =  c("year", "lat")
               , by =  c("year", "lat")
               , adjust = 'fdr')

#pairs(regrid(means1), simple = "each", combine = TRUE, adjust="fdr")
#pairs(regrid(means), simple = "each", combine = TRUE, adjust="fdr")
#emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")



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
  labs(x = "lat",
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
  labs(x = "lat",
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
lat=median(gt.full.range$latitude[which(gt.full.range$range=="Europe")])
lat
ggemmeans(glm, terms = c("time", "range", "latitude [lat]"), cov.reduce=range)
lat=median(gt.full.range$latitude[which(gt.full.range$range=="North America")])
lat
ggemmeans(glm, terms = c("time", "range", "latitude [lat]"), cov.reduce=range)
11:14
For instance these are
11:14
> lat=median(gt.full.range$latitude[which(gt.full.range$range=="Europe")])
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
lat=median(gt.full.range$latitude[which(gt.full.range$range=="North America")])
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
  labs(x = "lat",
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
  labs(x = "lat",
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

### no interactions
#means1 = (emtrends(glm, "range", var = "lat", adjust="fdr"))
#
#
##emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")
#
#means2 = emmeans(glm
#               , type = "pairs"
#               , specs =  c("year", "lat")
#               , by =  c("year", "lat")
#               , adjust = 'fdr')
#
##pairs(regrid(means1), simple = "each", combine = TRUE, adjust="fdr")
##pairs(regrid(means), simple = "each", combine = TRUE, adjust="fdr")
##emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")
#
#
##emtrends(glm, "range", var = "lat", adjust="fdr")
##means<-emmeans(glm
##               , type = "pairs"
##               , specs =  c("time", "range")
##               , by =  c("time", "range")
##               , adjust = 'fdr')


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
  labs(x = "lat",
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
  labs(x = "lat",
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









chr = "h1s12"
chr.start = 11570001
chr.end = 11910000


cnv.gt.mod = read.table("h1s12:11570001-11910000_mod_gt.txt", header = T)
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




cnv.gt.hist = read.table("h1s12:11570001-11910000_hist_gt.txt", header = T)
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
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + pc1, family=binomial)
Anova(glm, type = 3)
an_12a = Anova(glm, type = 3)


means1 = (emtrends(glm, "range", var = "lat", adjust="fdr"))


#emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")

means2 = emmeans(glm
               , type = "pairs"
               , specs =  c("year", "lat")
               , by =  c("year", "lat")
               , adjust = 'fdr')

#pairs(regrid(means1), simple = "each", combine = TRUE, adjust="fdr")
#pairs(regrid(means), simple = "each", combine = TRUE, adjust="fdr")
#emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")


#emtrends(glm, "range", var = "lat", adjust="fdr")
#means<-emmeans(glm
#               , type = "pairs"
#               , specs =  c("time", "range")
#               , by =  c("time", "range")
#               , adjust = 'fdr')





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
  labs(x = "lat",
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
  labs(x = "lat",
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


cnv12a_model_nonsig = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)










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



means1_an_14a = (emtrends(glm, pairwise ~ lat:range, var = "year", adjust="fdr"))
means2_an_14a = (emtrends(glm, "range", var = "lat", adjust="fdr"))
means3_an_14a = (emtrends(glm, "year", var = "lat", adjust="fdr"))





#emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")

means2 = emmeans(glm
               , type = "pairs"
               , specs =  c("year", "lat")
               , by =  c("year", "lat")
               , adjust = 'fdr')

#pairs(regrid(means1), simple = "each", combine = TRUE, adjust="fdr")
#pairs(regrid(means), simple = "each", combine = TRUE, adjust="fdr")
#emtrends(test, pairwise ~ time:range, var = "latitude", adjust="fdr")


#emtrends(glm, "range", var = "lat", adjust="fdr")
#means<-emmeans(glm
#               , type = "pairs"
#               , specs =  c("time", "range")
#               , by =  c("time", "range")
#               , adjust = 'fdr')




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
  labs(x = "lat",
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
  labs(x = "lat",
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
an_17a = Anova(glm, type = 3)



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
  labs(x = "lat",
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
  labs(x = "lat",
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
  labs(x = "lat",
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
  labs(x = "lat",
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










chr = "h1s18"
chr.start = 31240001
chr.end = 31960000



cnv.gt.mod = read.table("h1s18:31240001-31960000_mod_gt.txt", header = T)
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




cnv.gt.hist = read.table("h1s18:31240001-31960000_hist_gt.txt", header = T)
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

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:lat + range:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:lat + pc1, family=binomial)


glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + pc1, family=binomial)
Anova(glm, type = 3)

an = Anova(glm, type = 3)

kbl(an)


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
  labs(x = "lat",
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
  labs(x = "lat",
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


cnv18b_model_nonsig = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)






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
  labs(x = "lat",
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
  labs(x = "lat",
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

kbl(an)


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
  labs(x = "lat",
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
  labs(x = "lat",
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
  labs(x = "lat",
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
  labs(x = "lat",
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

an_17d = Anova(glm, type = 3)



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
  labs(x = "lat",
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
  labs(x = "lat",
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
cnv17d_model = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)











chr = "h1s6"
chr.start = 29680001
chr.end = 30030000



cnv.gt.mod = read.table("h1s6:29680001-30030000_mod_gt.txt", header = T)
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




cnv.gt.hist = read.table("h1s6:29680001-30030000_hist_gt.txt", header = T)
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

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + range:lat + year:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + year:lat + pc1, family=binomial)
Anova(glm, type = 3)

glm<-glm(data = cnv.gt.group, cbind(h1count, h2count) ~ year + range + lat + pc1, family=binomial)
an = Anova(glm, type = 3)



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
  labs(x = "lat",
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
  labs(x = "lat",
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


cnv6a_model_nonsig = plot_grid(cnv3.predicted.na, cnv3.predicted.eu, ncol = 2, align = "h", axis = "bt", labels = c("", ""), label_size = 36, vjust = 1)





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
  labs(x = "lat",
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
  labs(x = "lat",
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


png("cnv_predicted_sig.png", height=2000, width=1500)
plot_grid(cnv4a_model, cnv4b_model, cnv9a_model, cnv14a_model, cnv17a_model, cnv18a_model, cnv5a_model, cnv17c_model, cnv8a_model, cnv17d_model, cnv10a_model, ncol = 2, align = "vh")
dev.off()
















png("cnv_predicted_sig_FIG2.png", height=1000, width=1500)
plot_grid(cnv14a_model,cnv10a_model, ncol = 1, align = "v")
dev.off()


write.table(file = "an_4a.txt", t(an_4a)[c(1,3),], sep = "\t")
write.table(file = "an_4b.txt", t(an_4b)[c(1,3),], sep = "\t")
write.table(file = "an_9a.txt", t(an_9a)[c(1,3),], sep = "\t")
write.table(file = "an_14a.txt", t(an_14a)[c(1,3),], sep = "\t")
write.table(file = "an_17a.txt", t(an_17a)[c(1,3),], sep = "\t")
write.table(file = "an_18a.txt", t(an_18a)[c(1,3),], sep = "\t")
write.table(file = "an_5a.txt", t(an_5a)[c(1,3),], sep = "\t")
write.table(file = "an_17b.txt", t(an_17c)[c(1,3),], sep = "\t")
write.table(file = "an_8a.txt", t(an_8a)[c(1,3),], sep = "\t")
write.table(file = "an_17c.txt", t(an_17d)[c(1,3),], sep = "\t")
write.table(file = "an_10a.txt", t(an_10a)[c(1,3),], sep = "\t")





../all_modern_samples.list
while read sample; do
i=1
cat ${sample}_h1s.depth | awk -v chr=h1s 'BEGIN{print "ID\twindow\tdepth"} NR>1{print chr,$1,$2}' > ${sample}_all.depth
for i in {2..18}; do cat ${sample}_h1s.depth | awk -v chr=h1s${i} 'NR>1{print chr,$1,$2}' >> ${sample}_all.depth; done
done < ../combined_samples.list




setwd("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30_histmod")



file.names = dir(path = "../depth_out_merged_100kb_all_Q30", pattern = "_all.depth", full.names = TRUE)

# make df of average depth for each gene and sample
window_depth = read.table(file.names[1], header = F)
window_depth$normdepth = window_depth$V4/window_depth$V6
window_depth.df = cbind(window_depth[1], window_depth[2], window_depth[3], window_depth[7])

for (i in 2:length(file.names)) {
	TEMP = read.table(file.names[i], header = F)
	TEMP$normdepth = TEMP$V4/TEMP$V6
	window_depth.df = cbind(window_depth.df, TEMP[,7])
}


# name columns with samples
column.names = str_replace_all(file.names, "_all.depth", "")
column.names = str_replace_all(column.names, "../depth_out_merged_100kb_all_Q30/", "")

colnames(window_depth.df)[4:ncol(window_depth.df)] = column.names

write.table(window_depth.df, file = "10kb_window_depth_merged_100kb_histmod_nrml.txt", quote = F, row.names = F)





file.names = dir(path = "../depth_out_merged_100kb_all_Q30", pattern = "_all.depth", full.names = TRUE)

# make df of average depth for each gene and sample
window_depth = read.table(file.names[1], header = F)
window_depth$normdepth = window_depth$V4/window_depth$V6
window_depth.df = cbind(window_depth[1], window_depth[2], window_depth[3], window_depth[4])

for (i in 2:length(file.names)) {
	TEMP = read.table(file.names[i], header = F)
	TEMP$normdepth = TEMP$V4/TEMP$V6
	window_depth.df = cbind(window_depth.df, TEMP[,4])
}


# name columns with samples
column.names = str_replace_all(file.names, "_all.depth", "")
column.names = str_replace_all(column.names, "../depth_out_merged_100kb_all_Q30/", "")

colnames(window_depth.df)[4:ncol(window_depth.df)] = column.names

write.table(window_depth.df, file = "10kb_window_depth_merged_100kb_histmod_unnrml.txt", quote = F, row.names = F)





window_depth.df = read.table("10kb_window_depth_merged_100kb_histmod_nrml.txt", header = T)
bla = str_replace_all(colnames(window_depth.df), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")
colnames(window_depth.df) = bla
window_depth.df = window_depth.df[,!(colnames(window_depth.df) %in% c("MB-104"))]

rownames(window_depth.df) = paste0(window_depth.df$V1, ":", window_depth.df$V2, "-",  window_depth.df$V3)
window_depth.df[,1:3] = NULL


sample.medians=read.table("../depth_out_merged_100kb_all_Q30/sample_medians_histmod.list")
rownames(sample.medians)=sample.medians[,1]
sample.medians[,1] = NULL

windows.depth = as.data.frame(t(window_depth.df))

windows.depth.med = merge(windows.depth, sample.medians, by=0)
rownames(windows.depth.med) = windows.depth.med[,1]
windows.depth.med[,1]=NULL

combined.pops=read.csv("../combined_pops.csv")
rownames(combined.pops) = combined.pops[,1]
combined.pops[,1]=NULL
combined.pops.mod = combined.pops[which(combined.pops$type == "modern"),]
combined.pops.hist = combined.pops[which(combined.pops$type == "historic"),]

combined.pops[which(combined.pops$range == "eu" | combined.pops$range == "Europe"),][2] = "eu"
combined.pops[which(combined.pops$range == "East" | combined.pops$range == "middle" | combined.pops$range == "na" | combined.pops$range == "South" | combined.pops$range == "South-west" | combined.pops$range == "west" | combined.pops$range == "West"),][2] = "na"

wind.mod.ee = windows.depth.med[which(rownames(windows.depth.med) %in% rownames(combined.pops.mod)),]
wind.hist.ee = windows.depth.med[which(rownames(windows.depth.med) %in% rownames(combined.pops.hist)),]

wind.mod.ee$V2 = NULL
wind.hist.ee$V2 = NULL





cnv4split = read.table("all_chr4del_split_first_histmod.depth", header = F)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL

sample_medians=read.table("sample_medians_histmod.list")
rownames(sample_medians) = sample_medians[,1]
sample_medians[,1] = NULL

cnv4split = merge(cnv4split, sample_medians, by = 0)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL

cnv4split$depthnrml = cnv4split$V5/cnv4split$V2.y
cnv4split=merge(cnv4split, combined.pops, by = 0)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL
cnv4ab=cnv4split
cnv4ab[which(cnv4ab$type == "historic"),]$depthnrml = cnv4ab[which(cnv4ab$type == "historic"),]$depthnrml+0.052





cnv4split = read.table("all_chr4del_split_first_histmod.depth", header = F)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL

sample_medians=read.table("sample_medians_histmod.list")
rownames(sample_medians) = sample_medians[,1]
sample_medians[,1] = NULL

cnv4split = merge(cnv4split, sample_medians, by = 0)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL

cnv4split$depthnrml = cnv4split$V5/cnv4split$V2.y
cnv4split=merge(cnv4split, combined.pops, by = 0)
rownames(cnv4split) = cnv4split[,1]
cnv4split[,1] = NULL
cnv4aa=cnv4split
cnv4aa[which(cnv4aa$type == "historic"),]$depthnrml = cnv4aa[which(cnv4aa$type == "historic"),]$depthnrml+0.052



analysis.regions=c("h1s1:20430001-21570000", "h1s2:57390001-57530000", "h1s3:20640001-23650000", "h1s4:30540001-42390000", "h1s4:53800001-54500000", "h1s13:21790001-22680000", "h1s14:18280001-18750000", "h1s17:32800001-33040000", "V2")





windows.depth.med.int = windows.depth.med[,which(colnames(windows.depth.med) %in% analysis.regions)]

all.pops = c("ames", "apel", "berlin", "bord", "caluire", "florida", "inns", "ns", "prague", "quebec", "stgalmier", "stlouis")

wind.mod = windows.depth.med.int[which(rownames(windows.depth.med.int) %in% rownames(combined.pops.mod)),]
wind.hist = windows.depth.med.int[which(rownames(windows.depth.med.int) %in% rownames(combined.pops.hist)),]


for (i in 1:(ncol(wind.hist)-1)) {
	wind.hist[,i] = wind.hist[,i]+0.0514 
}

windows.depth.med.int = rbind(wind.hist,wind.mod)

	
wind.mod$time="modern"
wind.hist$time="historic"
	
nwind=ncol(wind.mod)-2
	
window_depth.pop=as.data.frame(rbind(wind.mod, wind.hist))

window_depth.pop=merge(window_depth.pop, combined.pops, by=0)
rownames(window_depth.pop)=window_depth.pop[,1]
window_depth.pop[,1]=NULL
tcol = ncol(window_depth.pop)
meds=nwind+1
	
window_depth.pop.all = window_depth.pop

for (i in 1:nwind){
	dfi = as.data.frame(cbind(window_depth.pop.all[,i], window_depth.pop.all[,tcol]))
}

window_depth.pop.hist = window_depth.pop.all[which(window_depth.pop.all$type == "historic"),]
window_depth.pop.mod = window_depth.pop.all[which(window_depth.pop.all$type == "modern"),]




































































neutral <- as.matrix(read.table("neutral_reg.cov"))
samples <- read.table("../combined_samples.list")
neutral.pcs <- eigen(neutral)
neutral.pcs <- as.data.frame(neutral.pcs$vectors[,1:2])
neutral.pcs <- cbind(samples, neutral.pcs)
colnames(neutral.pcs) <- c("V1", "neut_PC1", "neut_PC2")
rownames(neutral.pcs) = neutral.pcs[,1]
neutral.pcs[,1] = NULL

cnv3 = read.csv("cnv_3_genotypes.csv")
rownames(cnv3) = cnv3[,1]
cnv3[,1] = NULL
cnv3 = merge(cnv3, neutral.pcs, by = 0)
rownames(cnv3) = cnv3[,1]
cnv3[,1] = NULL

cnv4a = read.csv("cnv_4a_genotypes.csv")
rownames(cnv4a) = cnv4a[,1]
cnv4a[,1] = NULL
cnv4a = merge(cnv4a, neutral.pcs, by = 0)
rownames(cnv4a) = cnv4a[,1]
cnv4a[,1] = NULL

cnv4b = read.csv("cnv_4b_genotypes.csv")
rownames(cnv4b) = cnv4b[,1]
cnv4b[,1] = NULL
cnv4b = merge(cnv4b, neutral.pcs, by = 0)
rownames(cnv4b) = cnv4b[,1]
cnv4b[,1] = NULL

cnv13 = read.csv("cnv_13_genotypes.csv")
rownames(cnv13) = cnv13[,1]
cnv13[,1] = NULL
cnv13 = merge(cnv13, neutral.pcs, by = 0)
rownames(cnv13) = cnv13[,1]
cnv13[,1] = NULL

cnv14 = read.csv("cnv_14_genotypes.csv")
rownames(cnv14) = cnv14[,1]
cnv14[,1] = NULL
cnv14 = merge(cnv14, neutral.pcs, by = 0)
rownames(cnv14) = cnv14[,1]
cnv14[,1] = NULL

cnv3.lm = lm(depth ~ range*year*lat + neut_PC1, data = cnv3)
cnv3.lm = lm(depth ~ range+year+lat+range:lat+lat:year+year:range+ neut_PC1, data = cnv3)
cnv3.lm = lm(depth ~ range+year+lat+lat:year+year:range+ neut_PC1, data = cnv3)
cnv3.lm = lm(depth ~ range+year+lat+year:range+ neut_PC1, data = cnv3)
cnv3.lm = lm(depth ~ range+year+lat+ neut_PC1, data = cnv3)

summary(cnv3.lm)
coef = as.data.frame(t(summary(cnv3.lm)$coefficients[,3:4]))


Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  2.5347290  0.7354480   3.447 0.000609 ***
rangena     -0.3515760  0.0520738  -6.751 3.58e-11 ***
year        -0.0002603  0.0003483  -0.747 0.455191
lat         -0.0176419  0.0041923  -4.208 2.99e-05 ***
neut_PC1     3.7041965  0.5252775   7.052 5.08e-12 ***



cnv4a.lm = lm(depth ~ range*year*lat+ neut_PC1, data = cnv4a)
cnv4a.lm = lm(depth ~ range+year+lat+range:lat+lat:year+range:year+ neut_PC1, data = cnv4a)
cnv4a.lm = lm(depth ~ range+year+lat+range:lat+lat:year+ neut_PC1, data = cnv4a)
summary(cnv4a.lm)
coef = cbind(coef,(as.data.frame(t(summary(cnv4a.lm)$coefficients[,3:4]))))


Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept) -4.538e+00  2.133e+00  -2.127  0.03384 *
rangena      4.540e-01  1.753e-01   2.590  0.00985 **
year         2.516e-03  1.106e-03   2.276  0.02321 *
lat          1.010e-01  4.691e-02   2.154  0.03169 *
neut_PC1    -1.001e+00  2.193e-01  -4.565 6.12e-06 ***
rangena:lat -9.030e-03  3.802e-03  -2.375  0.01788 *
year:lat    -4.981e-05  2.437e-05  -2.044  0.04141 *




cnv4b.lm = lm(depth ~ range*year*lat+ neut_PC1, data = cnv4b)
cnv4b.lm = lm(depth ~ range+year+lat+range:lat+lat:year+range:year+ neut_PC1, data = cnv4b)
cnv4b.lm = lm(depth ~ range+year+lat+range:lat+lat:year+ neut_PC1, data = cnv4b)
cnv4b.lm = lm(depth ~ range+year+lat+range:lat+ neut_PC1, data = cnv4b)
summary(cnv4b.lm)

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept) -0.6375292  0.6634967  -0.961 0.337025
rangena      1.0000644  0.3261229   3.067 0.002267 **
year         0.0004633  0.0002770   1.672 0.094999 .
lat          0.0188156  0.0053608   3.510 0.000483 ***
neut_PC1    -4.7176338  0.4280000 -11.023  < 2e-16 ***
rangena:lat -0.0150163  0.0070605  -2.127 0.033861 *


cnv13.lm = lm(depth ~ range*year*lat+ neut_PC1, data = cnv13)
cnv13.lm = lm(depth ~ range+year+lat+range:lat+lat:year+range:year+ neut_PC1, data = cnv13)
cnv13.lm = lm(depth ~ range+year+lat+lat:year+range:year+ neut_PC1, data = cnv13)
cnv13.lm = lm(depth ~ range+year+lat+lat:year+ neut_PC1, data = cnv13)

summary(cnv13.lm)
summary(cnv13.lm)
coef = cbind(coef,(as.data.frame(t(summary(cnv13.lm)$coefficients[,3:4]))))


Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept) -1.717e+00  1.380e+00  -1.244  0.21415
rangena     -1.438e-03  1.364e-02  -0.105  0.91604
year         1.025e-03  7.057e-04   1.453  0.14682
lat          5.904e-02  3.005e-02   1.965  0.04990 *
neut_PC1     4.556e-01  1.395e-01   3.266  0.00115 **
year:lat    -3.085e-05  1.539e-05  -2.005  0.04548 *




cnv14.lm = lm(depth ~ range*year*lat+ neut_PC1, data = cnv14)
cnv14.lm = lm(depth ~ range+year+lat+range:lat+lat:year+range:year+ neut_PC1, data = cnv14)
cnv14.lm = lm(depth ~ range+year+lat+lat:year+range:year+ neut_PC1, data = cnv14)
cnv14.lm = lm(depth ~ range+year+lat+range:year+ neut_PC1, data = cnv14)
cnv14.lm = lm(depth ~ range+year+lat+ neut_PC1, data = cnv14)
summary(cnv14.lm)
coef = cbind(coef,(as.data.frame(t(summary(cnv14.lm)$coefficients[,3:4]))))

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept) -0.6595723  0.3505416  -1.882 0.060397 .
rangena      0.0816890  0.0248203   3.291 0.001059 **
year         0.0006206  0.0001660   3.738 0.000204 ***
lat          0.0061507  0.0019982   3.078 0.002182 **
neut_PC1     1.3675293  0.2503666   5.462 7.01e-08 ***













cnv4ab = merge(cnv4ab, neutral.pcs, by = 0)
rownames(cnv4ab) = cnv4ab[,1]
cnv4ab[,1] = NULL


cnv4ab.lm = lm(depthnrml ~ range*year*lat + neut_PC1, data = cnv4ab)
summary(cnv4ab.lm)

cnv4ab.lm = lm(depthnrml ~ range+year+lat+range:lat+lat:year+range:year+ neut_PC1, data = cnv4ab)
summary(cnv4ab.lm)

cnv4ab.lm = lm(depthnrml ~ range+year+lat+range:lat+lat:year + neut_PC1, data = cnv4ab)
summary(cnv4ab.lm)
cnv4ab.lm = lm(depthnrml ~ range+year+lat+range:lat + neut_PC1, data = cnv4ab)
summary(cnv4ab.lm)

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept) -1.7018101  0.2765340  -6.154 1.42e-09 ***
rangena      0.4281469  0.1362570   3.142  0.00176 **
year         0.0010411  0.0001154   9.021  < 2e-16 ***
lat          0.0047139  0.0022394   2.105  0.03573 *
neut_PC1    -1.1243736  0.1788311  -6.287 6.38e-10 ***
rangena:lat -0.0082825  0.0029497  -2.808  0.00516 **












cnv4aa = merge(cnv4aa, neutral.pcs, by = 0)
rownames(cnv4aa) = cnv4aa[,1]
cnv4aa[,1] = NULL


cnv4aa.lm = lm(depthnrml ~ range*year*lat + neut_PC1, data = cnv4aa)
summary(cnv4aa.lm)

cnv4aa.lm = lm(depthnrml ~ range+year+lat+range:lat+lat:year+range:year+ neut_PC1, data = cnv4aa)
summary(cnv4aa.lm)

cnv4aa.lm = lm(depthnrml ~ range+year+lat+lat:year+range:year + neut_PC1, data = cnv4aa)
summary(cnv4aa.lm)

cnv4aa.lm = lm(depthnrml ~ range+year+lat+lat:year + neut_PC1, data = cnv4aa)
summary(cnv4aa.lm)


coef = as.data.frame(t(summary(cnv3.lm)$coefficients[,3:4]))
coef = cbind(coef,(as.data.frame(t(summary(cnv4aa.lm)$coefficients[,3:4]))))
coef = cbind(coef,(as.data.frame(t(summary(cnv4b.lm)$coefficients[,3:4]))))
coef = cbind(coef,(as.data.frame(t(summary(cnv13.lm)$coefficients[,3:4]))))
coef = cbind(coef,(as.data.frame(t(summary(cnv14.lm)$coefficients[,3:4]))))




cnv3.lm.CONT = ggemmeans(cnv3.lm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
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
  ggtitle("CNV3 na") +
  labs(x = "lat",
    y = "") +
  xlim(20, 50) +
  ylim(0.5, 1.7) +
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
  ggtitle("CNV3 eu") +
  labs(x = "lat",
    y = "") +
  xlim(40, 70) +
  ylim(0.5, 1.7) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )




cnv4aa.lm.CONT = ggemmeans(cnv4aa.lm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv4aa.lm.CONT = as.data.frame(cnv4aa.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv4aa.predicted.na = ggplot() +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV4aa na") +
  labs(x = "lat",
    y = "") +
  xlim(20, 50) +
  ylim(0, 0.8) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv4aa.predicted.eu = ggplot() +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv4aa.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV4aa eu") +
  labs(x = "lat",
    y = "") +
  xlim(40, 70) +
  ylim(0, 0.8) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )








cnv4a.lm.CONT = ggemmeans(cnv4ab.lm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv4a.lm.CONT = as.data.frame(cnv4a.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv4a.predicted.na = ggplot() +
  geom_line(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV4a na") +
  labs(x = "lat",
    y = "") +
  xlim(20, 50) +
  ylim(0.25, 0.8) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv4a.predicted.eu = ggplot() +
  geom_line(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv4a.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV4a eu") +
  labs(x = "lat",
    y = "") +
  xlim(40, 70) +
  ylim(0.25, 0.8) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )









cnv4b.lm.CONT = ggemmeans(cnv4b.lm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv4b.lm.CONT = as.data.frame(cnv4b.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv4b.predicted.na = ggplot() +
  geom_line(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV4b na") +
  labs(x = "lat",
    y = "") +
  xlim(20, 50) +
  ylim(0.5, 1.8) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv4b.predicted.eu = ggplot() +
  geom_line(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv4b.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV4b eu") +
  labs(x = "lat",
    y = "") +
  xlim(40, 70) +
  ylim(0.5, 1.8) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )









cnv13.lm.CONT = ggemmeans(cnv13.lm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv13.lm.CONT = as.data.frame(cnv13.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv13.predicted.na = ggplot() +
  geom_line(data = subset(cnv13.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV13 na") +
  labs(x = "lat",
    y = "") +
  xlim(20, 50) +
  ylim(0, 0.5) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv13.predicted.eu = ggplot() +
  geom_line(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv13.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV13 eu") +
  labs(x = "lat",
    y = "") +
  xlim(40, 70) +
  ylim(0, 0.5) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )






cnv14.lm.CONT = ggemmeans(cnv14.lm, terms = c("lat", "range", "year [1830:2019 by=40]"), cov.reduce=range)
cnv14.lm.CONT = as.data.frame(cnv14.lm.CONT)

one.col <- "#762a83"
two.col <- "#af8dc3"
three.col <- "#d2aed3"
four.col <- "#7fbf7b"
five.col <- "#1b7837"

cnv14.predicted.na = ggplot() +
  geom_line(data = subset(cnv14.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "na" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "na" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "na" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "na" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "na" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV14 na") +
  labs(x = "lat",
    y = "") +
  xlim(20, 50) +
  ylim(0.5, 1.1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )


cnv14.predicted.eu = ggplot() +
  geom_line(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, y = predicted), color = one.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1830"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = one.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, y = predicted), color = two.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1870"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = two.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, y = predicted), color = three.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1910"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = three.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, y = predicted), color = four.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1950"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = four.col, alpha = 0.1) +
  geom_line(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, y = predicted), color = five.col) +
  geom_ribbon(data = subset(cnv14.lm.CONT, group == "eu" & facet == "1990"),
    aes(x = x, ymin = conf.low, ymax = conf.high), fill = five.col, alpha = 0.1) +
  ggtitle("CNV14 eu") +
  labs(x = "lat",
    y = "") +
  xlim(40, 70) +
  ylim(0.5, 1.1) +
  #ylim(0, 1) +
  theme_classic() +
  theme(legend.position = "none",
    plot.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 18)
  )



png("cnv_predicted.png", height=1000, width=1000)
grid.arrange(cnv3.predicted.na, cnv3.predicted.eu, cnv4aa.predicted.na, cnv4aa.predicted.eu, cnv4a.predicted.na, cnv4a.predicted.eu, cnv4b.predicted.na, cnv4b.predicted.eu, cnv13.predicted.na, cnv13.predicted.eu, cnv14.predicted.na, cnv14.predicted.eu, ncol=2)
dev.off()



#############################################################################################################################################################################
###################################################################   MERGED CNV GO ANALYSIS   ##############################################################################
#############################################################################################################################################################################




cd ~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/cnv_results_final

while read line1 line2 line3;  do cat qst_EU_1fst.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start <= $4 && end >= $4 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}' >> ${line1}:${line2}-${line3}_eu_qst.txt;  done < cnvs.list


while read line1 line2 line3;  do cat qst_EU_1fst.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$1 && start <= $2+5 && end >= $2-5) print $0}' >> ${line1}:${line2}-${line3}_eu_qst.txt;  done < cnvs.list
while read line1 line2 line3;  do cat qst_NA_1fst.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$1 && start <= $2+5 && end >= $2-5) print $0}' >> ${line1}:${line2}-${line3}_na_qst.txt;  done < cnvs.list



while read line4 line1 line2 line3;  do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | awk -v chr=${line1} -v start=${line2} -v end=${line3} '{if(chr==$3 && start >= $4 && start <= $5 || chr==$3 && end >= $4 && end <= $5 || chr==$3 && start >= $4 && end <= $5 || chr==$3 && start <= $4 && end >= $5)print $0}' >> ${line4}_genes.list;  done < cnvs.list 



cat merged_cnvs_genes.list | awk '{print $3, $4, $5}' > merged_cnvs_genes.bed


# create GO input file
cat merged_cnvs_genes.bed | awk '{print $1, "\t", $2-1, "\t", $3}' > merged_cnvs_genes_minus1.bed

# create GO input file
while read col1 col2 col3; 
do cat ~/om62_scratch/ragweed2022/fun/ragweed-dipasm-hap1-annot-GO-1e-6.txt | \
awk -F "\t" -v var1=$col1 -v var2=$col2 -v var3=$col3 '{if ($3 == var1 && $5 == var3) {print $2}}'; done < merged_cnvs_genes_minus1.bed | sort -u > merged_cnvs_genes.topGO



cp ../overlaps_q30/allannot.topGO .
### R


library(topGO)
setwd("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/cnv_results_final")



### EU

geneID2GO = readMappings(file = "allannot.topGO")
geneUniverse = names(geneID2GO)


genesOfInterest = read.table("merged_cnvs_genes.topGO")
genesOfInterest = as.character(genesOfInterest$V1)
geneList = factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) = geneUniverse
# ontology = 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component)
myGOdata = new("topGOdata",
               description = "My project",
               ontology = "BP",
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)
resultFisher = runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
allResSNP = GenTable(myGOdata,
                  classicFisher = resultFisher,
                  orderBy = "resultFisher",
                  ranksOf = "classicFisher",
                  topNodes = 100)
allResSNP
write.table(x = allResSNP, file = "merged_cnvs_genes.GO.txt", sep = "\t", quote = F)







#############################################################################################################################################################################
###################################################################   MERGED CNV TRAIT ANALYSIS   ##############################################################################
#############################################################################################################################################################################



setwd("~/om62_scratch/ragweed2022/cnv/gatk/samtools_depth/results_Q30_histmod")


window_depth.df = read.table("10kb_window_depth_merged_100kb_histmod_nrml.txt", header = T)
bla = str_replace_all(colnames(window_depth.df), "\\.", "-")
bla = str_replace_all(bla, "X28", "28")
colnames(window_depth.df) = bla
window_depth.df = window_depth.df[,!(colnames(window_depth.df) %in% c("MB-104"))]

rownames(window_depth.df) = paste0(window_depth.df$V1, ":", window_depth.df$V2, "-",  window_depth.df$V3)
window_depth.df[,1:3] = NULL


sample.medians=read.table("../depth_out_merged_100kb_all_Q30/sample_medians_histmod.list")
rownames(sample.medians)=sample.medians[,1]
sample.medians[,1] = NULL

windows.depth = as.data.frame(t(window_depth.df))

windows.depth.med = merge(windows.depth, sample.medians, by=0)
rownames(windows.depth.med) = windows.depth.med[,1]
windows.depth.med[,1]=NULL

combined.pops=read.csv("../combined_pops.csv")
rownames(combined.pops) = combined.pops[,1]
combined.pops[,1]=NULL
combined.pops.mod = combined.pops[which(combined.pops$type == "modern"),]
combined.pops.hist = combined.pops[which(combined.pops$type == "historic"),]

wind.mod.ee = windows.depth.med[which(rownames(windows.depth.med) %in% rownames(combined.pops.mod)),]
wind.hist.ee = windows.depth.med[which(rownames(windows.depth.med) %in% rownames(combined.pops.hist)),]

wind.mod.ee$V2 = NULL
wind.hist.ee$V2 = NULL


analysis.regions=c("h1s1:20430001-21570000", "h1s2:57390001-57530000", "h1s3:20640001-23650000", "h1s4:30540001-42390000", "h1s4:53800001-54500000", "h1s13:21790001-22680000", "h1s14:18280001-18750000", "h1s17:32800001-33040000", "V2")


windows.depth.med.int = windows.depth.med[,which(colnames(windows.depth.med) %in% analysis.regions)]


traits.df = read.table("../traits.list", header = F)
traits = as.character(traits.df$V1)


pheno=traits[1]
assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
trait.info = get(pheno)
rownames(trait.info) = trait.info$V1
colnames(trait.info) = c("sample", "sample2", paste0(pheno))

window_depth.traits = merge(windows.depth.med.int, trait.info[3], by = 0)
rownames(window_depth.traits) <- window_depth.traits[,1]
window_depth.traits[,1] <- NULL

for (pheno in traits[2:length(traits)]) {
		assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
		trait.info = get(pheno)
		rownames(trait.info) = trait.info$V1
		colnames(trait.info) = c("sample", "sample2", paste0(pheno))
		
		window_depth.traits = merge(window_depth.traits, trait.info[3], by = 0)
		rownames(window_depth.traits) <- window_depth.traits[,1]
		window_depth.traits[,1] <- NULL
}


neutral.pcs.gwas=read.table("../neutral-pc-gwas.txt", header=T)
window_depth.traits = merge(window_depth.traits, neutral.pcs.gwas, by=0)
rownames(window_depth.traits) = window_depth.traits[,1]
window_depth.traits[,1] <- NULL

pc1=ncol(window_depth.traits)-1
pc2=ncol(window_depth.traits)
tr.len=length(traits)+3
end.ind=ncol(window_depth.traits)-tr.len
windows.list=colnames(window_depth.traits[1:end.ind])


trait.pvals = data.frame(cnv = colnames(window_depth.traits[1:end.ind]), window = rep(NA, 8))


# iterate through traits list
for (j in 10:40){
	print(paste0("beginning phenotype ", colnames(window_depth.traits[j])))

	output=matrix(NA,nrow=0,ncol=3)
	for (i in 1:8) {
		# make temp coverage df for each window as well as trait value
		wd.temp=cbind(window_depth.traits[i], window_depth.traits[j], window_depth.traits[pc1])
		# calculate lm between window and trait
		res = lm(as.numeric(wd.temp[,2]) ~ as.numeric(wd.temp[,1]) + as.numeric(wd.temp[,3]))
		pval = summary(res)$coefficients[2,4]
		rsq = summary(res)$adj.r.squared
		out.line=c(colnames(wd.temp[1]),pval,rsq)
		# add lm for trait and window to output, then move onto next window
		output=rbind(output,out.line)
	}
	output.df = as.data.frame(output)
	k=j-8
	trait.pvals[k] = output.df$V2

	write.table(as.data.frame(output), file = paste0(traits[j], "merged_pheno_assoc.txt"), quote = F, row.names = F, col.names = F)
}


for (i in 2:ncol(trait.pvals)){trait.pvals[,i] = as.numeric(as.character(trait.pvals[,i])) * 8}

traits.cn = c("window", traits)
colnames(trait.pvals) = traits.cn
rownames(trait.pvals) = trait.pvals[,1]
trait.pvals[,1] = NULL

trait.pvals.sig = trait.pvals[,which(colSums(trait.pvals<=0.05)>=1)]


                          T_fl_end T_rac_die_day T_Raceme_length_longest T_seed_weight_tot T_sexmismatch T_sexrat_weight
h1s1:20430001-21570000  6.16121303    1.55514735              0.64706473        0.18172487    0.03914080    3.425421e+00
h1s2:57390001-57530000  6.48544126    3.93628363              1.38191793        0.89312500    0.55196378    7.263228e+00
h1s3:20640001-23650000  7.11852908    7.74322859              1.67801772        3.51492759    0.36411251    1.357283e+00
h1s4:30540001-36000000  6.221998      7.441752                0.0425661         0.4198381     4.953323      3.529761
h1s4:53800001-54500000  4.43642827    1.82236001              0.02784214        0.02701908    0.51664113    1.828378e+00
h1s13:21790001-22680000 0.04435216    0.02354078              2.80867493        4.25278018    0.79578290    4.797252e+00
h1s14:18280001-18750000 0.34593838    0.71946260              2.44947307        5.30057535    0.02320371    1.010075e-05
h1s17:32800001-33040000 4.23592093    7.85818982              3.68598881        0.60478073    0.61360221    3.400611e+00










window_depth = read.table("all_chr4del_split_first.depth", header = F)
rownames(window_depth) = window_depth[,1]

sample_medians=read.table("sample_medians_histmod.list")
rownames(sample_medians) = sample_medians[,1]

window_depth = merge(window_depth, sample_medians, by = 0)
window_depth$depthnrml = window_depth$V5/window_depth$V2.y

windows.depth.med.int = as.data.frame(cbind(as.character(window_depth$V1.x), window_depth$depthnrml))
rownames(windows.depth.med.int) = windows.depth.med.int[,1]
windows.depth.med.int[,1] = NULL



traits.df = read.table("../traits.list", header = F)
traits = as.character(traits.df$V1)









pheno=traits[1]
assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
trait.info = get(pheno)
rownames(trait.info) = trait.info$V1
colnames(trait.info) = c("sample", "sample2", paste0(pheno))

window_depth.traits = merge(windows.depth.med.int, trait.info[3], by = 0)
rownames(window_depth.traits) <- window_depth.traits[,1]
window_depth.traits[,1] <- NULL






for (pheno in traits[2:length(traits)]) {
		assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
		trait.info = get(pheno)
		rownames(trait.info) = trait.info$V1
		colnames(trait.info) = c("sample", "sample2", paste0(pheno))
		
		window_depth.traits = merge(window_depth.traits, trait.info[3], by = 0)
		rownames(window_depth.traits) <- window_depth.traits[,1]
		window_depth.traits[,1] <- NULL
}



neutral.pcs.gwas=read.table("neutral-pc-gwas.txt", header=T)
window_depth.traits = merge(window_depth.traits, neutral.pcs.gwas, by=0)
rownames(window_depth.traits) = window_depth.traits[,1]
window_depth.traits[,1] <- NULL

pc1=ncol(window_depth.traits)-1
pc2=ncol(window_depth.traits)
tr.len=length(traits)+3
end.ind=ncol(window_depth.traits)-tr.len
windows.list=colnames(window_depth.traits[1:end.ind])


#trait.pvals = data.frame(cnv = colnames(window_depth.traits[1:end.ind]), window = rep(NA, 8))
alloutput=matrix(NA,nrow=0,ncol=3)

# iterate through traits list
for (j in 2:32){
	print(paste0("beginning phenotype ", colnames(window_depth.traits[j])))

	output=matrix(NA,nrow=0,ncol=3)

	# make temp coverage df for each window as well as trait value
	wd.temp=cbind(window_depth.traits[1], window_depth.traits[j], window_depth.traits[pc1])
	# calculate lm between window and trait
	res = lm(as.numeric(wd.temp[,2]) ~ as.numeric(wd.temp[,1]) + as.numeric(wd.temp[,3]))
	pval = summary(res)$coefficients[2,4]
	rsq = summary(res)$adj.r.squared
	out.line=c(colnames(window_depth.traits[j]),pval,rsq)
	# add lm for trait and window to output, then move onto next window
	output=rbind(output,out.line)
	alloutput=rbind(alloutput,output)
	output.df = as.data.frame(output)
	k=j-8
	#trait.pvals[k] = output.df$V2

	write.table(as.data.frame(output), file = paste0(traits[j], "merged_pheno_assoc_chr4_first.txt"), quote = F, row.names = F, col.names = F)
}










window_depth = read.table("all_chr4del_split_weird.depth", header = F)
rownames(window_depth) = window_depth[,1]

sample_medians=read.table("sample_medians_histmod.list")
rownames(sample_medians) = sample_medians[,1]

window_depth = merge(window_depth, sample_medians, by = 0)
window_depth$depthnrml = window_depth$V5/window_depth$V2.y

windows.depth.med.int = as.data.frame(cbind(as.character(window_depth$V1.x), window_depth$depthnrml))
rownames(windows.depth.med.int) = windows.depth.med.int[,1]
windows.depth.med.int[,1] = NULL



traits.df = read.table("../traits.list", header = F)
traits = as.character(traits.df$V1)









pheno=traits[1]
assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
trait.info = get(pheno)
rownames(trait.info) = trait.info$V1
colnames(trait.info) = c("sample", "sample2", paste0(pheno))

window_depth.traits = merge(windows.depth.med.int, trait.info[3], by = 0)
rownames(window_depth.traits) <- window_depth.traits[,1]
window_depth.traits[,1] <- NULL






for (pheno in traits[2:length(traits)]) {
		assign(pheno, read.table(paste0("~/om62_scratch/ragweed2022/gwas/", pheno, ".WGS.NEW.pheno")))
		trait.info = get(pheno)
		rownames(trait.info) = trait.info$V1
		colnames(trait.info) = c("sample", "sample2", paste0(pheno))
		
		window_depth.traits = merge(window_depth.traits, trait.info[3], by = 0)
		rownames(window_depth.traits) <- window_depth.traits[,1]
		window_depth.traits[,1] <- NULL
}



neutral.pcs.gwas=read.table("neutral-pc-gwas.txt", header=T)
window_depth.traits = merge(window_depth.traits, neutral.pcs.gwas, by=0)
rownames(window_depth.traits) = window_depth.traits[,1]
window_depth.traits[,1] <- NULL

pc1=ncol(window_depth.traits)-1
pc2=ncol(window_depth.traits)
tr.len=length(traits)+3
end.ind=ncol(window_depth.traits)-tr.len
windows.list=colnames(window_depth.traits[1:end.ind])


#trait.pvals = data.frame(cnv = colnames(window_depth.traits[1:end.ind]), window = rep(NA, 8))
alloutput=matrix(NA,nrow=0,ncol=3)

# iterate through traits list
for (j in 2:32){
	print(paste0("beginning phenotype ", colnames(window_depth.traits[j])))

	output=matrix(NA,nrow=0,ncol=3)

	# make temp coverage df for each window as well as trait value
	wd.temp=cbind(window_depth.traits[1], window_depth.traits[j], window_depth.traits[pc1])
	# calculate lm between window and trait
	res = lm(as.numeric(wd.temp[,2]) ~ as.numeric(wd.temp[,1]) + as.numeric(wd.temp[,3]))
	pval = summary(res)$coefficients[2,4]
	rsq = summary(res)$adj.r.squared
	out.line=c(colnames(window_depth.traits[j]),pval,rsq)
	# add lm for trait and window to output, then move onto next window
	output=rbind(output,out.line)
	alloutput=rbind(alloutput,output)
	output.df = as.data.frame(output)
	k=j-8
	#trait.pvals[k] = output.df$V2

	write.table(as.data.frame(output), file = paste0(traits[j], "merged_pheno_assoc_chr4_weird.txt"), quote = F, row.names = F, col.names = F)
}












