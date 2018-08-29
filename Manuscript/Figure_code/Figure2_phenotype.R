# Yuka Takemon
# 08/29/18
# Phenotype distribution plots for publication and statistics
library(tidyverse)
library(ggsignif)
library(ggsci)
library(grid)

setwd("/projects/korstanje-lab/ytakemon/Col4a5xDO/")
load("./Col4a5xDO_192data_YT.Rdata")

#clean data
rownames(Pheno) <- make.names(Pheno[,1]) #move sample ID to row names
Pheno <- Pheno[rownames(Genoprobs),] #subset Pheno to match 192 samples
#clean up Pheno and add log of ACR
Pheno[Pheno < 0 ] = NA
Pheno[Pheno ==  -Inf] = NA
Pheno[Pheno ==  -Inf] = NA
options(na.action = 'na.pass') #leave in NAs
#keep only interested columns

# Phenotype: ACR
PhenoACR <- Pheno[,c("MouseID", "Sex", "ACR6WK", "ACR10WK", "ACR15WK")] %>%
	rename(Week6 = ACR6WK,
				 Week10 = ACR10WK,
			 	 Week15 = ACR15WK) %>%
	gather(key = ACRWK, value = Measurement,
		Week6, Week10, Week15) %>%
	mutate(Sex = as.factor(Sex),
				 ACRWK = as.character(ACRWK))

PhenoACR[PhenoACR$ACRWK == "Week6",]$ACRWK <- 6
PhenoACR[PhenoACR$ACRWK == "Week10",]$ACRWK <- 10
PhenoACR[PhenoACR$ACRWK == "Week15",]$ACRWK <- 15

PhenoACR_plot <- PhenoACR %>%
	na.omit() %>%
	mutate(ACRWK = factor(ACRWK, levels = c("6","10","15"))) %>%
	ggplot(., aes(x = ACRWK, y = log(Measurement))) +
		geom_violin(alpha = 0.3, aes(fill =Sex)) +
		geom_boxplot(width = 0.1, aes(fill =Sex)) +
		geom_signif(comparisons = list(c("6","10")),
								y_position = 9.5,
								map_signif_level = TRUE,
								test = "wilcox.test") +
		geom_signif(comparisons = list(c("10","15")),
								y_position = 9.9,
								map_signif_level = TRUE,
								test = "wilcox.test") +
		geom_signif(comparisons = list(c("6","15")),
								y_position = 10.5,
								map_signif_level = TRUE,
								test = "wilcox.test") +
		scale_y_continuous("Log (ACR mg/g)", breaks= seq(2,12,2), limit = c(2,11))+
		labs(x = "Weeks") +
		guides(fill = FALSE) +
		facet_wrap(~ Sex, strip.position = "left", nrow =2) +
		theme_bw()+
		scale_fill_aaas()

# ANOVA
anovF <- PhenoACR %>%
	filter(Sex == "F") %>%
	na.omit() %>%
	aov(Measurement ~ ACRWK, data = .) %>%
	TukeyHSD(.,"ACRWK")

anovM <- PhenoACR %>%
	filter(Sex == "M") %>%
	na.omit() %>%
	aov(Measurement ~ ACRWK, data = .) %>%
	TukeyHSD(.,"ACRWK")

PhenoACR %>%
	group_by(Sex, ACRWK) %>%
	na.omit() %>%
	summarise(mean_ACR = mean(Measurement),
						log_mean_ACR = log(mean(Measurement)))

# Phenotype: GFR
PhenoGFR <- Pheno[,c("MouseID", "Sex", "C2")] %>%
	#mutate(Sex = factor(Sex, levels = c("M","F"))) %>%
	ggplot(., aes(x = Sex, y = log(C2), fill = Sex)) +
		geom_violin(alpha = 0.3) +
		geom_boxplot(width = 0.2) +
		labs(x = "Sex")+
		scale_y_continuous(expression("Log (GFR "*mu*"l/min) at 14 weeks"),
			breaks= seq(3,7.5,0.5), limit = c(3,7.5))+
		facet_wrap(~ Sex, strip.position = "left", nrow =2, scales="free") +
		theme_bw()+
		guides(fill = FALSE)+
		scale_fill_aaas()
# Average GFR
Pheno[,c("MouseID", "Sex", "C2")] %>%
	na.omit() %>%
	group_by(Sex) %>%
	summarise(meanGFR = mean(C2),
						log_meanGFR = log(mean(C2)))
# T-test
Pheno[,c("MouseID", "Sex", "C2")] %>%
	na.omit() %>%
	t.test(C2 ~ Sex, data = .)

# Phenotype correlation
cor <- Pheno %>%
	select(MouseID, Sex, C2, ACR15WK) %>%
	na.omit() %>%
	filter(C2 < 1000) %>%
	group_by(Sex) %>%
	summarise(R2 = signif(cor(ACR15WK, C2, method = "pearson")^2,2),
						pval = signif(cor.test(ACR15WK, C2, method = "pearson")$p.value, 3))

# if pvalue is nominally small
for(i in nrow(cor)){
	if (cor[i,"pval"] < 0.01){
		cor[i,"pval"] <- "< 0.01"
	}
}

# plot correlation
PhenoCor_plot <- Pheno %>%
	select(MouseID, Sex, C2, ACR15WK) %>%
	na.omit() %>%
	filter(C2 < 1000) %>%
	ggplot(., aes(x = log(C2), y = log(ACR15WK), colour = Sex)) +
		geom_point()+
		geom_smooth(method = "lm", se = FALSE, colour = "black")+
		labs(x = expression("log (GFR "*mu*"l/min) at 14 weeks"), y = "log (ACR mg/g) at 15 weeks")+
		facet_wrap(~Sex, strip.position = "left", nrow = 2)+
		theme_bw()+
		annotate(geom = "text",
			x = 5,
			y = 10,
			label =  c(paste0("R^2 == ",cor[1,"R2"]), paste0("R^2 == ",cor[2,"R2"])),
			hjust = 0,
			parse = TRUE) +
		annotate(geom = "text",
			x = 5,
			y = 9.7,
			label =  c(paste0("p-value = ",cor[1,"pval"]), paste0("p-value = ",cor[2,"pval"])),
			hjust = 0,
			parse = FALSE)+
		guides(colour = FALSE) +
		scale_color_aaas()

pdf("~/Dropbox/Col4a5xDO_Manuscript/Manuscript/Figures_pdf/Figure2_Phenotypes.pdf", width = 12, height = 6.5)
pushViewport(viewport(layout = grid.layout(1, 3,
											widths = unit(c(25,10,20), c("lines","lines","lines")),
										  height = unit(32, "lines"))))
print(PhenoACR_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(PhenoGFR, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(PhenoCor_plot, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
dev.off()

#ACR Females
#pheno_F <- pheno_F[,c("MouseID", "ACR6WK", "ACR10WK", "ACR15WK")]
#ggdata<- melt(pheno_F)
#names(ggdata) <- c("MouseID", "ACR_time", "Value")
#ggplot_F <- ggplot(ggdata, aes(Value, ..count.., fill = ACR_time, colour = ACR_time)) +
#	geom_density( alpha = 0.1)+
#	scale_x_continuous("Albumin to creatinine ratio (mg/g)", limit = c(0 , 15000)) +
#	scale_y_continuous(" ") +
#	labs( title = "Females") +
#	theme( plot.title = element_text(hjust = 0.5),
#			panel.background = element_rect(fill = "white", colour = "black"),
#			axis.text.x = element_blank(),
#			axis.ticks.x = element_blank()
#			) +
#	coord_flip()

#pheno_M <- pheno_M[,c("MouseID", "ACR6WK", "ACR10WK", "ACR15WK")]
#ggdata<- melt(pheno_M)
#names(ggdata) <- c("MouseID", "ACR_time", "Value")
#ggplot_M <- ggplot(ggdata, aes(Value, ..count.., fill = ACR_time, colour = ACR_time)) +
#	geom_density( alpha = 0.1)+
#	scale_x_continuous("Albumin to creatinine ratio (mg/g)", limit = c(0 , 15000)) +
#	scale_y_continuous(" ") +
#	labs( title = "Males") +
#	theme( plot.title = element_text(hjust = 0.5),
#			panel.background = element_rect(fill = "white", colour = "black"),
#			axis.text.x = element_blank(),
#			axis.ticks.x = element_blank()
#			) +
#	coord_flip()

#pdf("./GBRS_reconstruction/reconstruct/best.compiled.genoprob/plot/Col4a5_fig2_ACR_GFR_dist_YUKA.pdf", width = 12, height = 6)
#pushViewport(viewport(layout = grid.layout(1, 2)))
#print(ggplot_F, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
#print(ggplot_M, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
#dev.off()
