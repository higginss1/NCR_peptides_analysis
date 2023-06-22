library(tidyverse)
library(ggplot2)
library(vegan)
library(growthrates)
library(ggridges)
library(gridExtra)
library(cowplot)
library(FrF2)

# get the path to the R input data files containing Liberibacter crescens strain BT-1 growth data
growthdata <- Sys.glob(file.path(getwd(), "Lcre_growth_assay_182_NCR_peptides", "NCR_Peptides_*", "*_R_data_input.txt"))

# get the biomatik peptide data in long format
biomatiklong <- read.table(file.path(getwd(), "Lcre_growth_assay_182_NCR_peptides", "Biomatik_182_NCR_peptide_data_long.txt"), sep="\t", header=T, stringsAsFactors = T)
# get the biomatik peptide data in wide format
biomatikwide <- read.table(file.path(getwd(), "Lcre_growth_assay_182_NCR_peptides", "Biomatik_182_NCR_peptide_data_wide.txt"), sep="\t", header=T, stringsAsFactors = T)

# john peptide data
johnpeptide <- read.table(file.path(getwd(), "Lcre_growth_assay_182_NCR_peptides", "all_peptides_with_GRAVY_score_edit.txt"), sep="\t", header=T, stringsAsFactors = T)


# read in growth data from NCR peptide MIC assay 96-well plate data
growthdf <- data.frame()
for(g in 1:length(growthdata)){
  tmpdata <- read.table(growthdata[[g]], sep="\t", header=T, stringsAsFactors = T)
  tmpdata$plate <- rep(g, dim(tmpdata)[1])
  growthdf <- rbind(growthdf, tmpdata)
}

# merge data sets by wide and long formats for different analyses
growthdfwide <- merge(growthdf, biomatikwide, by.x = "Sample_ID", by.y = "Catalog_number", all.x = T)
growthdflong <- merge(growthdf, biomatiklong, by.x = "Sample_ID", by.y = "Catalog_number", all.x = T)

# now use dplyr to create a new data column with relative growth compared to positive growth control
growthdflongrel <- growthdflong %>% group_by(Sample_ID, Concentration, Timepoint, plate) %>% mutate(meanAbs = mean(Absorbance), sdAbs = sd(Absorbance)) %>% as.data.frame()

# I figured out how to create the normalized absorbance I was looking for (relative to the negative control)
growthdflongrel %>% 
  group_by(plate, Timepoint) %>% 
  mutate(tmp = mean(Absorbance[Sample_ID == "L_cre_BM7"])) %>% 
  mutate(normAbs = tmp - meanAbs) %>% 
  select(normAbs, everything())

# remove NCR peptides that were likely insoluble (starting O.D. > 0.3)
growthdflongrel2 <- growthdflongrel[which(!growthdflongrel$Sample_ID %in% unique(growthdflongrel[which(growthdflongrel$Timepoint == 0 & growthdflongrel$Absorbance > 0.3),]$Sample_ID)),]
levels(growthdflongrel2$Type) <- c("Media only", "No inhibitor", "inhibitor", "NCR peptides")

# plot all soluble NCR peptides, inhibitors, and media only controls compared to the 'red' no inhibitor control
ggplot(growthdflongrel2) + 
  geom_point(aes(x = Timepoint, y = meanAbs, color = Type, shape = Type)) + 
  geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = Type)) + 
  theme_light() + 
  theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black')) + 
  facet_wrap(~plate, ncol=4) + guides(color='none') + scale_color_manual(values=c('grey', 'red', 'grey', 'grey')) + 
  ylab(expression("Abs"["600nm"])) + 
  xlab('Time (days)')


#############
# FIGURE 1A #
#############

split.data <- multisplit(growthdf, c("Sample_ID", "Type", "Concentration", "Well_ID", "plate"))
dat <- split.data[[1]]
dat2 <- split.data[[31]]
fit <- fit_spline(dat$Timepoint, dat$Absorbance)
summary(fit)
coef(fit)
fit2 <- fit_spline(dat2$Timepoint, dat2$Absorbance)
summary(fit2)
coef(fit2)

# starting parameters provided to the parametric nonlinear growth model
# these parameters were derived from values observed in the total data set
p     <- c(y0 = 0.2, mumax = 0.2, K = 0.8)
lower <- c(y0 = 0.01, mumax = 0,   K = 0)
upper <- c(y0 = 0.5, mumax = 5,   K = 1.6)

# fit and plot fit_growthmodel fits
logdat <- NULL
for(i in 1:length(split.data)){
  tmpdat <- split.data[[i]]
  tmpfit <- fit_growthmodel(FUN = grow_logistic, p = p, tmpdat$Timepoint, tmpdat$Absorbance,
                            lower = lower, upper = upper)
  # can alternatively fit a heuristic method to generate data - this method also provides
  # a lag parameter which could be useful
  #tmpfit <- fit_easylinear(tmpdat$Timepoint, tmpdat$Absorbance, h = 5, quota = 0.95)
  tmpdat$y0 <- rep(coef(tmpfit)[1], dim(tmpdat)[1])
  tmpdat$y0_lm <- rep(coef(tmpfit)[2], dim(tmpdat)[1])
  tmpdat$mumax <- rep(coef(tmpfit)[3], dim(tmpdat)[1])
  tmpdat$rsquared <- rep(tmpfit@rsquared, dim(tmpdat)[1])
  # don't need lag at the moment
  #tmpdat$lag <- rep(coef(tmpfit)[4], dim(tmpdat)[1])
  tmpdat$Absorbance.pred <- as.data.frame(predict(tmpfit))$y
  logdat <- rbind(logdat, tmpdat)
  
}

# plot the data fits
logdatplots <- list()
plates <- unique(logdat$plate)
c <- 1
for(i in plates){
  tmplogdat <- logdat %>% filter(plate == i)
  samps <- unique(tmplogdat$Sample_ID)
  for(s in samps){
    
    tmplogdat2 <- tmplogdat %>% filter(Sample_ID %in% c(s, "L_cre_BM7"))
    tmpplot <- ggplot(tmplogdat2) + 
      geom_point(aes(x = Timepoint, y = Absorbance, color = Concentration, group = Concentration), alpha = 1) + 
      geom_line(aes(x = Timepoint, y = Absorbance.pred, group = Well_ID, color = Concentration, linetype = Concentration), size = 0.1) + 
      scale_color_brewer(palette = 'Set2', direction = -1) +
      theme_light() + 
      theme(axis.text = element_text(color = 'black'), plot.title = element_text(hjust=0.5), legend.position = c(0.25, 0.75)) + 
      xlab("Time (days)") + 
      ylab(expression('Absorbance'['600nm'])) + 
      ggtitle(s)
    logdatplots[[c]] <- tmpplot
    c <- c + 1
  }
}

mplot <- marrangeGrob(logdatplots, nrow=2, ncol=2)
ggsave(mplot, file = paste(getwd(), "Lcre_growth_assay_182_NCR_peptides", "growthrates_fitgrowthmodel_growlogistic_diagnostic_plots.pdf", sep = "/"), dpi = 300, width = 8, height = 6)

# now for each NCR peptide on each plate, get the average relative inhibition in mumax caused
# by each NCR peptide, if any
#logdat2 <- logdat %>% group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% summarise(mean_Abs = mean(Absorbance), sd_Abs = sd(Absorbance), mean_y0 = mean(y0), sd_y0 = sd(y0), meany0_lm = mean(y0_lm), sd_y0_lm = sd(y0_lm), mean_mumax = mean(mumax), sd_mumax = sd(mumax), mean_lag = mean(lag), sd_lag = sd(lag))
logdatB <- logdat %>% group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% summarise(mean_Abs = mean(Absorbance), sd_Abs = sd(Absorbance), mean_y0 = mean(y0), sd_y0 = sd(y0), meany0_lm = mean(y0_lm), sd_y0_lm = sd(y0_lm), mean_mumax = mean(mumax), sd_mumax = sd(mumax))
logdatC <- logdatB %>% group_by(plate, Type, Concentration, Sample_ID) %>% summarise(mean_y0 = mean(mean_y0), sd_y0 = sd(sd_y0), meany0_lm = mean(meany0_lm), sd_y0_lm = sd(sd_y0_lm), mean_mumax = mean(mean_mumax), sd_mumax = sd(sd_mumax))
#logdatD <- logdatC %>% group_by(plate) %>% mutate(percent_mumax = 100.0-((mean_mumax/mean_mumax[Type == "Negative_Control"])*100))
logdatD <- logdatC %>% group_by(plate) %>% mutate(norm_mumax = (mean_mumax[Type == "Negative_Control"] - mean_mumax)/mean_mumax[Type == "Negative_Control"], percent_mumax = 100.0*norm_mumax)

levels(logdatD$Type) <- c("Media only", "No inhibitor", "inhibitor", "NCR peptides")
logdatD <- logdatD[!rownames(logdatD) %in% which(logdatD$Type %in% c("Media only", "inhibitor") & logdatD$percent_mumax < 50),]

# plot data

onemgml <- logdatD %>% filter(Concentration %in% c('1_mg_ml', 'None', '0.5_mg_ml', '0.005_mg_ml')) %>%
  ggplot(aes(x = percent_mumax, y = Type, fill = Type)) + 
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.9, aes(point_fill = Type, point_color = Type)) + 
  theme_light() + 
  xlab(bquote(Growth~rate~"("*mu*")"~inhibition~"(%)")) +
  ggtitle('1 mg/ml NCR peptide') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.title.y = element_blank(), axis.text.x = element_text(color = 'black'))

ohonemgml <- logdatD %>% filter(Concentration %in% c('0.1_mg_ml', 'None', '0.5_mg_ml', '0.005_mg_ml')) %>% 
  ggplot(aes(x = percent_mumax, y = Type, fill = Type)) + 
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.9, aes(point_fill = Type, point_color = Type)) + 
  theme_light() + 
  xlab(bquote(Growth~rate~"("*mu*")"~inhibition~"(%)")) +
  ggtitle('0.1 mg/ml NCR peptide') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(color = 'black'))

figure1a <- plot_grid(onemgml, ohonemgml, nrow = 1, rel_widths = c(1.2,1))

# save the input data to figure 1
write.table(logdat, file = file.path(getwd(), "Lcre_growth_assay_182_NCR_peptides","Figure1_input.csv"), sep =",", quote = F, row.names = F)

# now plot the peptides with >40% inhibition for 1 mg/ml NCR peptides
ggplot(growthdflongrel2[which(growthdflongrel2$Sample_ID %in% c(as.vector(unique(logdatD[which(logdatD$percent_mumax>40.0),]$Sample_ID)), "L_cre_BM7") & growthdflongrel2$Concentration %in% c("None", "1_mg_ml", "0.5_mg_ml", "0.005_mg_ml")),]) + 
  geom_point(aes(x = Timepoint, y = meanAbs, color = Type, shape = Type)) + 
  geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = Type)) + 
  theme_light() + 
  theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black')) + 
  facet_wrap(~plate, ncol=4) + guides(color='none') + scale_color_brewer(palette = "Set2", direction = -1) + 
  ylab(expression("Abs"["600nm"])) + 
  xlab('Time (days)')

# now plot the peptides with >40% inhibition for 0.1 mg/ml NCR peptides
ggplot(growthdflongrel2[which(growthdflongrel2$Sample_ID %in% c(as.vector(unique(logdatD[which(logdatD$percent_mumax>40.0),]$Sample_ID)), "L_cre_BM7") & growthdflongrel2$Concentration %in% c("None", "0.1_mg_ml", "0.5_mg_ml", "0.005_mg_ml")),]) + 
  geom_point(aes(x = Timepoint, y = meanAbs, color = Type, shape = Type)) + 
  geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = Type)) + 
  theme_light() + 
  theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black')) + 
  facet_wrap(~plate, ncol=4) + guides(color='none') + scale_color_brewer(palette = "Set2", direction = -1) + 
  ylab(expression("Abs"["600nm"])) + 
  xlab('Time (days)')

# fit an plot fit_easylinear fits
logdat2 <- NULL
for(i in 1:length(split.data)){
  tmpdat2 <- split.data[[i]]
  matchkey <- unique(as.vector(tmpdat2)[1])
  if(matchkey == "PMB" | matchkey == "BM7_only" | matchkey == "Y4761"){
    simpreg <- lm(tmpdat2$Absorbance~tmpdat2$Timepoint, data=tmpdat2)
    simppred <- as.data.frame(predict(simpreg))
    colnames(simppred) <- c("Absorbance.pred")
    tmpdat2$y0 <- rep(0.2, dim(tmpdat2)[1])
    tmpdat2$y0_lm <- rep(0.01, dim(tmpdat2)[1])
    tmpdat2$mumax <- rep(0.01, dim(tmpdat2)[1])
    tmpdat2$lag <- rep(2, dim(tmpdat2)[1])
    tmpdat2$rsquared <- rep(tmpfit@rsquared, dim(tmpdat)[1])
    tmpdat2$Absorbance.pred <- simppred$Absorbance.pred
    logdat2 <- rbind(logdat2, tmpdat2)
    #logdat2plots[[i]] <- ggplot(tmpdat2) + geom_point(aes(x = Timepoint, y = Absorbance)) + geom_line(aes(x = Timepoint, y = Absorbance.pred, group = Well_ID), linetype = 2, size = 1) + theme_light() + theme(plot.title = element_text(hjust=0.5), axis.text = element_text(color = 'black')) + scale_y_log10() + xlab("Time (days)") + ylab(expression('log'[10]~'Absorbance'['600nm'])) + ggtitle(unique(as.vector(tmpdat2)[1])[1]$Sample_ID)
    
  }else{
    #tmpfit <- fit_growthmodel(FUN = grow_logistic, p = p, tmpdat$Timepoint, tmpdat$Absorbance,
    #                          lower = lower, upper = upper)
    # can alternatively fit a heuristic method to generate data - this method also provides
    # a lag parameter which could be useful
    tmpfit <- fit_easylinear(tmpdat2$Timepoint, tmpdat2$Absorbance, h = 5, quota = 0.95)
    tmpdat2$y0 <- rep(coef(tmpfit)[1], dim(tmpdat2)[1])
    tmpdat2$y0_lm <- rep(coef(tmpfit)[2], dim(tmpdat2)[1])
    tmpdat2$mumax <- rep(coef(tmpfit)[3], dim(tmpdat2)[1])
    tmpdat2$lag <- rep(coef(tmpfit)[4], dim(tmpdat2)[1])
    tmpdat2$rsquared <- rep(tmpfit@rsquared, dim(tmpdat)[1])
    tmpdat2$Absorbance.pred <- predict(tmpfit)$y
    logdat2 <- rbind(logdat2, tmpdat2)
    #logdat2plots[[i]] <- ggplot(tmpdat2) + geom_point(aes(x = Timepoint, y = Absorbance)) + geom_line(aes(x = Timepoint, y = Absorbance.pred, group = Well_ID), linetype = 2, size = 1) + theme_light() + theme(plot.title = element_text(hjust=0.5), axis.text = element_text(color = 'black')) + scale_y_log10() + xlab("Time (days)") + ylab(expression('log'[10]~'Absorbance'['600nm'])) + ggtitle(unique(as.vector(tmpdat2)[1])[1]$Sample_ID)
    
  }
}

logdat2a <- logdat2 %>% group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% summarise(mean_Abs = mean(Absorbance), sd_Abs = sd(Absorbance), mean_y0 = mean(y0), sd_y0 = sd(y0), meany0_lm = mean(y0_lm), sd_y0_lm = sd(y0_lm), mean_mumax = mean(mumax), sd_mumax = sd(mumax), mean_lag = mean(lag), sd_lag = sd(lag))
logdat3a <- logdat2a %>% group_by(plate, Type, Concentration, Sample_ID) %>% summarise(mean_y0 = mean(mean_y0), sd_y0 = sd(sd_y0), meany0_lm = mean(meany0_lm), sd_y0_lm = sd(sd_y0_lm), mean_mumax = mean(mean_mumax), sd_mumax = sd(sd_mumax), mean_lag = mean(mean_lag), sd_lag = sd(sd_lag))
logdat4a <- logdat3a %>% group_by(plate) %>% mutate(norm_mumax = mean_mumax[Type == "Negative_Control"] - mean_mumax, percent_mumax = (1 - (mean_mumax/mean_mumax[Type == "Negative_Control"]))*100.0)
logdat3c <- logdat2a %>% group_by(plate, Timepoint, Type, Concentration, Sample_ID) %>% summarise(mean_y0 = mean(mean_y0), sd_y0 = sd(sd_y0), meany0_lm = mean(meany0_lm), sd_y0_lm = sd(sd_y0_lm), mean_mumax = mean(mean_mumax), sd_mumax = sd(sd_mumax), mean_lag = mean(mean_lag), sd_lag = sd(sd_lag))
logdat4c <- logdat3c %>% group_by(plate) %>% mutate(norm_mumax = mean_mumax[Type == "Negative_Control"] - mean_mumax, percent_mumax = (1 - (mean_mumax/mean_mumax[Type == "Negative_Control"]))*100.0)

#############
# FIGURE 1B #
#############

highAbs <- growthdflongrel %>% filter(Absorbance > 0.35, Timepoint %in% c(0,1))  %>% select(Sample_ID) %>% filter(Sample_ID != "L_cre_BM7") %>% distinct()

growthdfa <- growthdf %>% group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% summarise(meanAbs = mean(Absorbance), sdAbs = sd(Absorbance)) %>% as.data.frame()
growthdfb <- growthdfa %>% filter(!Sample_ID %in% as.vector(highAbs$Sample_ID))
ggplot(growthdfb) + geom_point(aes(x = Timepoint, y = meanAbs, color = Type, shape = Type)) + geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = Type)) + theme_light() + theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black')) + facet_wrap(~plate, ncol=4) + guides(color=F) + scale_color_manual(values=c('black', 'red', 'black', 'grey')) + ylab(expression("Abs"["600nm"])) + xlab('Time (days)')
# remove peptides that don't have a low final mean Absorbance
endAbs <- growthdfb %>% filter((Timepoint == 7 & between(meanAbs, 0, 0.8)) | (Timepoint == 4 & between(meanAbs, 0.2, 0.8)))  %>% select(Sample_ID) %>% distinct()
growthdfd <- growthdfa %>% filter(Sample_ID %in% c(as.vector(endAbs$Sample_ID), "PMB", "L_cre_BM7", "Blank", "Y4761"))
ggplot(growthdfd[which(growthdfd$Concentration %in% c("0.1_mg_ml", "None", "0.5_mg_ml", "0.005_mg_ml")),]) + geom_point(aes(x = Timepoint, y = meanAbs, color = Type, shape = Type)) + geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = Type)) + theme_light() + theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black')) + facet_wrap(~plate, ncol=4) + guides(color=F) + scale_color_manual(values=c('black', 'red', 'black', 'grey')) + ylab(expression("Abs"["600nm"])) + xlab('Time (days)')
ggplot(growthdfd[which(growthdfd$Concentration %in% c("1_mg_ml", "None", "0.5_mg_ml", "0.005_mg_ml")),]) + geom_point(aes(x = Timepoint, y = meanAbs, color = Type, shape = Type)) + geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = Type)) + theme_light() + theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black')) + facet_wrap(~plate, ncol=4) + guides(color=F) + scale_color_manual(values=c('black', 'red', 'black', 'grey')) + ylab(expression("Abs"["600nm"])) + xlab('Time (days)')
growthdff <- growthdfb %>% group_by(plate, Timepoint) %>% mutate(diffAbs = meanAbs[Type == "Negative_Control"] - meanAbs)
endAbs2 <- growthdff %>% ungroup() %>% filter(diffAbs > 0.35, Timepoint == 6) %>% select(Sample_ID) %>% distinct()
growthdff2 <- growthdff %>% filter(Sample_ID %in% c(as.vector(endAbs2$Sample_ID), "PMB", "L_cre_BM7", "Blank", "Y4761"))
endAbs3 <- growthdff2 %>% ungroup() %>% filter(between(diffAbs,0.3, 2), Timepoint == 7) %>% select(Sample_ID) %>% distinct()
growthdff3 <- growthdff %>% filter(Sample_ID %in% c(as.vector(endAbs3$Sample_ID), "PMB", "L_cre_BM7", "Blank", "Y4761"))
endAbs4 <- growthdff3 %>% ungroup() %>% filter(Timepoint == 7, meanAbs < 0.7) %>% select(Sample_ID) %>% distinct()
growthdff4 <- growthdff %>% filter(Sample_ID %in% c(as.vector(endAbs4$Sample_ID), "PMB", "L_cre_BM7", "Blank", "Y4761"))

inhibNCRs <- growthdff4 %>% ungroup() %>% select(Sample_ID) %>% filter(!Sample_ID %in% c("PMB", "L_cre_BM7", "BM7_only", "Y4761")) %>% distinct()
levels.design <- rep(2, 50)


# 74 NCR peptides (~41%) were insoluble or displayed strange growth dynamics
# leaving a remaining 108 NCR peptides that were not insoluble
# get mean and sd of technical replicates of soluble NCR peptides and create a wide format of the data
# by using the timepoint as columns instead of a variable
growthdf2 <- growthdf %>% group_by(Concentration, Type, Timepoint, plate, Sample_ID) %>% summarize(meanAbs = mean(Absorbance), sdAbs = sd(Absorbance)) %>% filter(!Sample_ID %in% as.vector(highAbs$Sample_ID))
growthdf3 <- growthdf2 %>% select(Sample_ID, Timepoint, meanAbs) %>% spread(Timepoint, meanAbs)
growth.pca.input <- growthdf3[,c("0","1","2","3","4","5","6","7")]
growth.pca.input <- growth.pca.input %>% as.matrix()
rownames(growth.pca.input) <- growthdf3$Sample_ID
growth.pca <- prcomp(growth.pca.input, scale. = T)
growth.pca.scores <- scores(growth.pca)
growthdf3$PC1 <- growth.pca.scores[,1]
growthdf3$PC2 <- growth.pca.scores[,2]
growthdf3$PC3 <- growth.pca.scores[,3]
ggplot(growthdf3) + geom_point(aes(x = PC1, y = PC2, color = Type)) + theme_light()

# perform new pca analysis while adding in peptide chemical information
growthdf3a <- merge(growthdf3, biomatikwide, by.x = c("Sample_ID"), by.y = c("Catalog_number"), all.x = T)
growthdf3b <- merge(growthdf3a, johnpeptide, by.x = c("Top_20mer_AMP"), by.y = c("Top.20mer.AMP"), all.x = T)
growthdf3b <- growthdf3b %>% filter(!Sample_ID %in% c("L_cre_BM7", "Y4761", "PMB", "BM7_only"))
growthdf3c <- growthdf3b %>% select(c("0", "1", "2", "3", "4", "5", "6", "7", "GRAVY.score", "Number.Cysteines", "AMP.charge.at.pH.7")) %>% as.matrix()
colnames(growthdf3c) <- c("0", "1", "2", "3", "4", "5", "6", "7", "GRAVY", "Cysteines", "Charge")

for(i in 1:ncol(growthdf3c)){
  growthdf3c[is.na(growthdf3c[,i]), i] <- mean(growthdf3c[,i], na.rm = TRUE)
}
growthdf4a <- merge(growthdf3b, logdat4a, by = c("plate", "Concentration", "Sample_ID"), all.x = T)
growthdf4c <- growthdf4a %>% select(c("percent_mumax", "GRAVY.score", "Number.Cysteines", "AMP.charge.at.pH.7")) %>% as.matrix()
growthdf4c[which(growthdf4c[,1] < 0),][,1] <- 0
colnames(growthdf4c) <- c("Growth inhibition", "Peptide hydrophobicity", "No. Cysteines", "Peptide charge")


# collate the final set of most inhibitory NCR peptides observed out of those
# tested
# some number of these peptides will be tested for interactions using a fractional factorial design
inhibiNCRdf <- growthdf4a[which(growthdf4a$Sample_ID %in% unique(growthdff4$Sample_ID) & growthdf4a$Concentration == "1_mg_ml"),c("Sample_ID", "percent_mumax", "Top_20mer_AMP")]
inhibiNCRdf <- inhibiNCRdf[order(inhibiNCRdf$percent_mumax, decreasing=T),]
write.table(inhibiNCRdf, file = file.path(getwd(), "Lcre_growth_assay_182_NCR_peptides", "Top_47_NCR_peptide_inhibitors.txt"), sep = ',', quote = F, col.names = T, row.names = F)

growthdff4b <- growthdff4[which(growthdff4$Sample_ID %in% c(as.vector(inhibiNCRdf$Sample_ID[1:15]), "L_cre_BM7", "PMB", "BM7_only", "Y4761")),]
growthdff5b <- growthdff4b %>% group_by(Type, Concentration, Timepoint) %>% summarise(totmeanAbs = mean(meanAbs), totsdAbs = sd(meanAbs)) %>% as.data.frame()

growthdff5 <- growthdff4 %>% group_by(Type, Concentration, Timepoint) %>% summarise(totmeanAbs = mean(meanAbs), totsdAbs = sd(meanAbs)) %>% as.data.frame()

figure1b <- ggplot(growthdff5[which(growthdff5$Concentration %in% c("1_mg_ml", "None", "0.5_mg_ml", "0.005_mg_ml") & growthdff5$Type != "Positive_Control"),], aes(x = Timepoint, y = totmeanAbs, color = Type)) + 
  geom_line(aes(group = interaction(Type, Concentration), color = Type), size=2) + 
  geom_ribbon(aes(group = interaction(Type, Concentration), x = Timepoint, ymin=totmeanAbs - totsdAbs, ymax = totmeanAbs + totsdAbs, fill = Type), alpha = 0.25, color = NA) + 
  theme_light() + 
  theme(legend.position = 'none', axis.text = element_text(color = 'black', size=11), axis.title = element_text(size = 11, color = 'black')) + 
  xlab("Time (days)") + 
#  ylab(expression(italic('Liberibacter')~'growth')) + 
  ylab(expression(italic('L. crescens')~'growth ('~'Abs'['600nm']~')')) + 
  scale_color_brewer(palette = "Dark2", direction = -1) + 
  scale_fill_brewer(palette = "Dark2", direction = -1) + 
  annotate("text", x = 6, y = 1.1, label = sprintf("-NCR peptides"), color = '#D95F02', fontface = 2) + 
  annotate("text", x = 6, y = 0.6, label = "+NCR peptides", color = '#1B9E77',fontface = 2) + 
  annotate("text", x = 6, y = 0.25, label = "Media only", color = '#7570B3', fontface = 2) +
  annotate(geom = "curve", x = 3.25, y = 1, xend = 6, yend = 0.7, curvature = -0.5, arrow = arrow(length = unit(4, "mm")), color = '#1B9E77') + annotate(geom = "text", x = 2.3, y = 1, label = sprintf("%s\n %s\n%s", "Up to 73%", "growth rate", "suppression"), color = '#1B9E77', fontface = 2)

figure1 <- plot_grid(figure1a, figure1b, nrow=2, labels = c('A','B'))

ggsave(figure1, file = file.path(getwd(), "figures", "Figure1.pdf"), dpi = 300, width = 8, height = 8)
ggsave(figure1, file = file.path(getwd(), "figures", "Figure1.tiff"), dpi = 300, width = 8, height = 8)
ggsave(figure1, file = file.path(getwd(), "figures", "Figure1.png"), dpi = 300, width = 8, height = 8)
ggsave(figure1, file = file.path(getwd(), "figures", "Figure1.eps"), dpi = 300, width = 8, height = 8)
ggsave(figure1, file = file.path(getwd(), "figures", "Figure1.svg"), dpi = 300, width = 8, height = 8)

logdat2B <- logdat2 %>% group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% summarise(mean_Abs = mean(Absorbance), sd_Abs = sd(Absorbance), mean_y0 = mean(y0), sd_y0 = sd(y0), meany0_lm = mean(y0_lm), sd_y0_lm = sd(y0_lm), mean_mumax = mean(mumax), sd_mumax = sd(mumax))
logdat3B <- logdat2B %>% group_by(plate, Type, Concentration, Sample_ID) %>% summarise(mean_y0 = mean(mean_y0), sd_y0 = sd(sd_y0), meany0_lm = mean(meany0_lm), sd_y0_lm = sd(sd_y0_lm), mean_mumax = mean(mean_mumax), sd_mumax = sd(sd_mumax))
logdat4B <- logdat3B %>% group_by(plate) %>% mutate(norm_mumax = (mean_mumax[Type == "Negative_Control"] - mean_mumax)/mean_mumax[Type == "Negative_Control"], percent_mumax = norm_mumax*100.0)

ggplot(logdat4B[which(logdat4B$Concentration %in% c("1_mg_ml", "0.1_mg_ml")),], aes(x = percent_mumax)) + geom_histogram() + theme_bw() + xlab("Growth rate inhibition (%)") + facet_wrap(~Concentration)

levels(logdat4B$Type) <- c("Media only", "No inhibitor", "inhibitor", "NCR peptides")
logdat4B <- logdat4B[!rownames(logdat4B) %in% which(logdat4B$Type %in% c("Media only", "inhibitor") & logdat4B$percent_mumax < 50),]

# plot data

onemgmlefit <- logdat4B %>% filter(Concentration %in% c('1_mg_ml', 'None', '0.5_mg_ml', '0.005_mg_ml')) %>%
  ggplot(aes(x = percent_mumax, y = Type, fill = Type)) + 
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.9, aes(point_fill = Type, point_color = Type)) + 
  theme_light() + 
  xlab(bquote(Growth~rate~"("*mu*")"~inhibition~"(%)")) +
  ggtitle('1 mg/ml NCR peptide') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.title.y = element_blank())

ohonemgmlefit <- logdat4B %>% filter(Concentration %in% c('0.1_mg_ml', 'None', '0.5_mg_ml', '0.005_mg_ml')) %>% 
  ggplot(aes(x = percent_mumax, y = Type, fill = Type)) + 
  geom_density_ridges(jittered_points = TRUE, position = "raincloud", alpha = 0.7, scale = 0.9, aes(point_fill = Type, point_color = Type)) + 
  theme_light() + 
  xlab(bquote(Growth~rate~"("*mu*")"~inhibition~"(%)")) +
  ggtitle('0.1 mg/ml NCR peptide') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none', axis.title.y = element_blank(), axis.text.y = element_blank())

growthratedistefit <- plot_grid(onemgmlefit, ohonemgmlefit, nrow = 1)


logdat2plots <- list()
plates <- unique(logdat2$plate)
c <- 1
for(i in plates){
  tmplogdat <- logdat2 %>% filter(plate == i)
  samps <- unique(tmplogdat$Sample_ID)
  for(s in samps){
    
    tmplogdat2 <- tmplogdat %>% filter(Sample_ID %in% c(s, "L_cre_BM7"))
    #tmpplot <- ggplot(tmplogdat2) + geom_point(aes(x = Timepoint, y = Absorbance, color = Concentration, group = Concentration), alpha = 0.6) + geom_line(aes(x = Timepoint, y = Absorbance.pred, group = Well_ID, color = Concentration), linetype = 2, size = 1) + theme_light() + theme(axis.text = element_text(color = 'black'), plot.title = element_text(hjust=0.5), legend.position = c(0.25, 0.75)) + scale_y_log10() + xlab("Time (days)") + ylab(expression('log'[10]~'Absorbance'['600nm'])) + ggtitle(s)
    tmpplot <- ggplot(tmplogdat2) + 
      geom_point(aes(x = Timepoint, y = Absorbance, color = Concentration, group = Concentration), alpha = 1) + 
      geom_line(aes(x = Timepoint, y = Absorbance.pred, group = Well_ID, color = Concentration, linetype = Concentration), size = 0.1) + 
      scale_y_log10() +
      scale_color_brewer(palette = 'Set2', direction = -1) +
      theme_light() + 
      theme(axis.text = element_text(color = 'black'), plot.title = element_text(hjust=0.5), legend.position = c(0.25, 0.75)) + 
      xlab("Time (days)") + 
      ylab(expression('log'[10]~'Absorbance'['600nm'])) +
      ggtitle(s)
    logdat2plots[[c]] <- tmpplot
    c <- c + 1
  }
}

mplot <- marrangeGrob(logdat2plots, nrow=2, ncol=2)
#ggsave(mplot, file = "C:/Users/sah382/Documents/NCR_Peptide_MIC_Assay_test/growthrates_fiteasylinear_test.pdf", dpi = 300, width = 8, height = 6)
ggsave(mplot, file = paste(getwd(), "Lcre_growth_assay_182_NCR_peptides", "growthrates_fiteasylinear_diagnostic_plots.pdf", sep = "/"), dpi = 300, width = 8, height = 6)


# create a new label column in the data frame to plot sample_ID of those lines/points with growth inhibition
# in the graph
growthdflongrel3 <- growthdflongrel %>% group_by(plate, Timepoint) %>% mutate(tmp = mean(Absorbance[Sample_ID == "L_cre_BM7"])) %>% mutate(inhibition = ifelse(meanAbs < tmp, 'Yes', 'No'))

plate_names <- list('1' = 'Plate 1', '2' = 'Plate 2', '3' = "Plate 3", '4' = "Plate 4", '5' = "Plate 5", '6' = 'Plate 6', '7' = "Plate 7", '8' = "Plate 8", '9' = 'Plate 9', '10' = 'Plate 10', '11' = 'Plate 11', '12' = 'Plate 12', '13' = 'Plate 13', '14' = 'Plate 14', '15' = 'Plate 15', '16' = 'Plate 16')
labfunc <- function(variable,value){
  return(plate_names[value])
}
pepplot <- ggplot(growthdflongrel3[which(growthdflongrel3$Type %in% c("Treatment", "Negative_Control")),]) + geom_point(aes(x = Timepoint, y = meanAbs, color = inhibition, shape = Type)) + geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = inhibition)) + theme_light() + theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black'), legend.position = "bottom") + facet_wrap(~plate, ncol=4, labeller = labfunc) + guides(color=F) + scale_color_manual(values=c('grey', 'red')) + geom_text(data=growthdflongrel3[which(growthdflongrel3$Timepoint==7 & growthdflongrel3$inhibition == "Yes" & growthdflongrel3$Type != "Positive_Control" & growthdflongrel3$Type != "Blank"),], aes(x = Timepoint, y = meanAbs, label = Sample_ID), nudge_x = 0.6, hjust=0, size=2, fontface='bold', check_overlap = T) + xlab("Time (days)") + ylab(expression("Abs"["600nm"])) + scale_x_continuous(breaks = scales::pretty_breaks(), expand = expand_scale(mult = c(0, 0.5), add=c(2,0))) + scale_shape_manual(values = c(19,1))
ggsave(pepplot, file = file.path(getwd(), "figures", "NCR_Peptide_plates_1-16_inhibition.pdf"), dpi = 300, width = 8, height = 6)
ggsave(pepplot, file = file.path(getwd(), "figures", "NCR_Peptide_plates_1-16_inhibition.tiff"), dpi = 300, width = 8, height = 6)
ggsave(pepplot, file = file.path(getwd(), "figures", "NCR_Peptide_plates_1-16_inhibition.png"), dpi = 300, width = 8, height = 6)
ggsave(pepplot, file = file.path(getwd(), "figures", "NCR_Peptide_plates_1-16_inhibition.eps"), dpi = 300, width = 8, height = 6)
ggsave(pepplot, file = file.path(getwd(), "figures", "NCR_Peptide_plates_1-16_inhibition.svg"), dpi = 300, width = 8, height = 6)

# now filter the data to grab only day 7 time points and create a new normalized absorbance value
growthdflongrel3day7 <- growthdflongrel3 %>% group_by(plate) %>% filter(Timepoint == 7) %>% group_by(plate) %>% mutate(normAbs = meanAbs[Sample_ID == "L_cre_BM7"] - meanAbs)
ggplot(growthdflongrel3day7) + geom_bar(aes(x = Amino_acid, y = normAbs)) + theme_light() + facet_wrap(~Amino_acid_position, ncol=4)

# calculate the area between the negative control growth curve for a particular plate and NCR peptide
test <- growthdflongrel3 %>% group_by(plate, Sample_ID, Concentration, Timepoint) %>% summarise(Abs = mean(Absorbance))
plates <- unique(test$plate)
conc <- unique(test$Concentration)
samps <- unique(test$Sample_ID)
absdf <- NULL
for(p in plates){
  for(c in conc){
    for(s in samps){
      tmpdf <- data.frame('plate' = p, 'Concentration' = c, 'Sample_ID' = s)
      # use a function from geiger package to do this
      # NOTE: if f1 is larger than f2, the value will be positive
      # but for some reason, a flat line (e.g., PMB treatment) as f2 results in a negative value
      # value of f2 larger than f1 are supposed to generate negative values, not the other way around
      #tmpdf$area <- geiger:::.area.between.curves(test[which(test$Sample_ID=="L_cre_BM7" & test$plate == p),]$Timepoint, test[which(test$Sample_ID == "L_cre_BM7" & test$plate == p),]$Abs, test[which(test$plate == p & test$Sample_ID == s & test$Concentration == c),]$Abs, xrange=c(min(test[which(test$Sample_ID=="L_cre_BM7" & test$plate == p),]$Timepoint),max(test[which(test$Sample_ID=="L_cre_BM7" & test$plate == p),]$Timepoint)))
      tmpdf$area <- geiger:::.area.between.curves(test[which(test$Sample_ID=="L_cre_BM7" & test$plate == p),]$Timepoint, test[which(test$Sample_ID == "L_cre_BM7" & test$plate == p),]$Abs, test[which(test$plate == p & test$Sample_ID == s & test$Concentration == c),]$Abs, xrange=c(0,7))
      absdf <- rbind(absdf, tmpdf)
    }
  }
}
absdf <- na.omit(absdf)
absdflong <- merge(absdf, biomatiklong, by.x = "Sample_ID", by.y = "Catalog_number", all.x = T)
absdfwide <- merge(absdf, biomatikwide, by.x = "Sample_ID", by.y = "Catalog_number", all.x = T)

ggplot(absdflong[which(absdflong$Concentration %in% c("1_mg_ml", "0.1_mg_ml") & absdflong$area > 1.0),], aes(x = Amino_acid_position, y = Amino_acid, size = area, color = area)) + geom_point(shape = 1) + theme_light() + scale_color_viridis_c(direction=-1) + facet_grid(Concentration~plate)
ggplot(absdflong[which(absdflong$Concentration %in% c("1_mg_ml", "0.1_mg_ml") & absdflong$area > 1.0),], aes(x = Amino_acid_position, y = area, color = Concentration)) + geom_point(shape = 19) + theme_light() + scale_color_viridis_d(direction=-1) + facet_wrap(~Amino_acid)
ggplot(absdflong[which(absdflong$Concentration %in% c("1_mg_ml", "0.1_mg_ml") & absdflong$area > 1.0),], aes(y = Amino_acid_position, x = area)) + geom_density_ridges() + theme_light() + scale_color_viridis_d(direction=-1) + facet_wrap(~Amino_acid)


# perform linear regression to assess the interaction of amino acid position with a particular amino acid on area
summary(lm(area~Amino_acid_position + Amino_acid + Amino_acid_position*Amino_acid, data=absdflong))

############
# FIGURE 2 #
############

# fractional factorial design and analysis

# note: ggbiplot requires plyr package, which can affect behavior of tidyverse
#library(ggbiplot)

# ggridges also requires plyr package, which can affect behavior of dplyr if loaded after it
#library(ggridges)

library(tidyverse)
library(ggplot2)
library(vegan)
library(growthrates)
library(gridExtra)
library(cowplot)
library(ggrepel)

# get the path to the R input data files containing Liberibacter crescens strain BT-1 growth data
growthdata <- Sys.glob(file.path(getwd(), "NCR_Frac_Fac_128", "NCR_Peptides_*", "*_R_data_input.txt"))

# get the biomatik peptide data in long format
biomatiklong <- read.table(file.path(getwd(), "NCR_Frac_Fac_128", "Biomatik_182_NCR_peptide_data_long.txt"), sep="\t", header=T, stringsAsFactors = T)
# get the biomatik peptide data in wide format
biomatikwide <- read.table(file.path(getwd(), "NCR_Frac_Fac_128", "Biomatik_182_NCR_peptide_data_wide.txt"), sep="\t", header=T, stringsAsFactors = T)
# john peptide data
johnpeptide <- read.table(file.path(getwd(), "NCR_Frac_Fac_128", "all_peptides_with_GRAVY_score_edit.txt"), sep="\t", header=T, stringsAsFactors = T)


# read in growth data from NCR peptide MIC assay 96-well plate data
growthdf <- data.frame()
for(g in 1:length(growthdata)){
  tmpdata <- read.table(growthdata[[g]], sep="\t", header=T, stringsAsFactors = T)
  tmpdata$plate <- rep(g, dim(tmpdata)[1])
  growthdf <- rbind(growthdf, tmpdata)
}


# now use dplyr to create a new data column with relative growth compared to positive growth control
growthdflongrel <- growthdf %>% 
  group_by(Sample_ID, plate, Concentration, Timepoint) %>% 
  mutate(meanAbs = mean(Absorbance), sdAbs = sd(Absorbance)) %>% 
  as.data.frame()

# change level names of factor "Type"
levels(growthdflongrel$Type) <- c('Media only', 'No inhibitor',  'inhibitor', 'NCR peptides')
growthdflongrel$plate2 <- paste('Plate', growthdflongrel$plate, sep = " ")
growthdflongrel$plate2 <- factor(growthdflongrel$plate2, levels = paste("Plate", seq(1,11), sep = " "))

figure2 <- ggplot(growthdflongrel) + 
  geom_point(aes(x = Timepoint, y = meanAbs, color = Type, shape = Type)) + 
  geom_line(aes(x = Timepoint, y = meanAbs, group = interaction(Sample_ID, Concentration), color = Type)) + 
  theme_light() + theme(axis.text = element_text(color = 'black'), strip.background = element_rect(fill='black')) + 
  facet_wrap(~plate2, ncol=4, scales = "free_x") + guides(color=NULL) + 
  scale_color_manual(values=c('black', 'red', 'black', 'grey')) + 
  ylab(expression("Abs"["600nm"])) + 
  xlab('Time (days)')

ggsave(figure2, file = file.path(getwd(), "figures", "Figure2.pdf"), dpi = 300, width = 8, height = 6)
ggsave(figure2, file = file.path(getwd(), "figures", "Figure2.eps"), dpi = 300, width = 8, height = 6)
ggsave(figure2, file = file.path(getwd(), "figures", "Figure2.png"), dpi = 300, width = 8, height = 6)
ggsave(figure2, file = file.path(getwd(), "figures", "Figure2.svg"), dpi = 300, width = 8, height = 6)
ggsave(figure2, file = file.path(getwd(), "figures", "Figure2.tiff"), dpi = 300, width = 8, height = 6)

############
# Figure 3 #
############

library(tidyverse)
library(ggplot2)
library(readxl)
library(ggridges)
library(stringr)
library(cowplot)
library(extrafont)
library(grid)
library(gridExtra)

loadfonts(device = "win")

#####################################################################################################
# Function to add qPCR Well ID to amplification data output from ABI QuantStudio 6 Flex qPCR system #
#####################################################################################################

# for amplification data, I need to add Well_ID instead of a generic "Well" number column
# add a Well_ID column based on the numbers in the Well column of the indivampdata data.frame
c <- 1
nlist <- list()
# create a list with 384-well plate well IDs in it
for(i in 1:length(LETTERS[1:16])){for(j in 1:24){nlist <- append(nlist, paste(LETTERS[i], as.character(j), sep= "")); print(c(c, paste(LETTERS[i], as.character(j), sep= ""))); c <- c + 1}}
# add the well number to the names of this list to associated well ID with the
# well number from the qPCR instrument
names(nlist) <- 1:384


welladd <- function(nl, ampd){
  # create a new Well_ID column in ampd variable
  ampd$Well_ID <- rep("A1", dim(ampd)[1])
  # add Well_ID information
  for(i in 1:length(names(nl))){
    ampd[which(ampd$Well %in% names(nl)[i]),]$Well_ID <- nl[[i]][1]
  }
  # change columns names for simplicity
  colnames(ampd) <- c("Well", "Cycle", "Target_Name", "Rn", "Delta_Rn", "Well_ID")
  return(ampd)
}

# set ggplot theme

th <- theme(axis.text = element_text(color = 'black', size = 10, family = 'Arial'), axis.title = element_text(size = 16, color = 'black', family = 'Arial'), plot.title = element_text(hjust = 0.5), legend.position = 'none', strip.background = element_rect(fill = 'white'))


# estimated genome/plasmid size for gene copy estimation using standard curves

size <- list("wingless" = 485705082, 'GAPC2R' = 372000000, 'CLas_16S' = 3135, "EF1AFR" = 372000000)

# have the volume of DNA pipetted per assay

vol2 <- list("wingless" = 1, 'CLas_16S' = 1)

# read in qPCR sample order and well information (as well as sample metadata)

samples <- read.table(file = file.path(getwd(), "Detached_leaf_assay_data_inputs", "SH_NCR_Peptide_DL_DLA_total_2021_09_28_qPCR_sample_sheet_final.txt", sep = "/"), sep = "\t", header = T, stringsAsFactors = F)

# read in ABI QuantStudio 6 Flex qPCR system data

DNA_plate1 <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_DNA_plate_1_and_2_2021-10-2_12232.xls"), sheet = "Results")
DNA_plate1_amp <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_DNA_plate_1_and_2_2021-10-2_12232.xls"), sheet = "Amplification Data")
DNA_plate1_amp <- welladd(nlist, DNA_plate1_amp)

DNA_plate2 <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_DNA_plate_2 and_3_and_4_2021-10-22_141301.xls"), sheet = "Results")
DNA_plate2_amp <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_DNA_plate_2 and_3_and_4_2021-10-22_141301.xls"), sheet = "Amplification Data")
DNA_plate2_amp <- welladd(nlist, DNA_plate2_amp)

DNA_plate3 <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_DNA_plate_4_and_5_2021-10-22_160345.xls"), sheet = "Results")
DNA_plate3_amp <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_DNA_plate_4_and_5_2021-10-22_160345.xls"), sheet = "Amplification Data")
DNA_plate3_amp <- welladd(nlist, DNA_plate3_amp)

cDNA_plate1 <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_cDNA_plate1_and_2_CLas_16S_2021-10-21_114549.xls"), sheet = "Results")
cDNA_plate1_amp <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_cDNA_plate1_and_2_CLas_16S_2021-10-21_114549.xls"), sheet = "Amplification Data")
cDNA_plate1_amp <- welladd(nlist, cDNA_plate1_amp)

cDNA_plate2 <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_peptide_DL_DLA_cDNA_plates_2_3_and_4_repeat_2021-10-22_101044.xls"), sheet = "Results")
cDNA_plate2_amp <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_peptide_DL_DLA_cDNA_plates_2_3_and_4_repeat_2021-10-22_101044.xls"), sheet = "Amplification Data")
cDNA_plate2_amp <- welladd(nlist, cDNA_plate2_amp)

cDNA_plate3 <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_cDNA_plate4_and_5_CLas_16S_2021-10-21_162717.xls"), sheet = "Results")
cDNA_plate3_amp <- read_excel(file.path(getwd(), "Detached_leaf_assay_data_inputs","SAH382_NCR_Peptide_DL_DLA_cDNA_plate4_and_5_CLas_16S_2021-10-21_162717.xls"), sheet = "Amplification Data")
cDNA_plate3_amp <- welladd(nlist, cDNA_plate3_amp)

# extract specific plate information from total sample data
DNA_samp1 <- samples[which(samples$qPCR_plate_ID == '1' & samples$Molecule == 'DNA'),]
DNA_samp2 <- samples[which(samples$qPCR_plate_ID == '2' & samples$Molecule == 'DNA'),]
DNA_samp3 <- samples[which(samples$qPCR_plate_ID == '3' & samples$Molecule == 'DNA'),]

cDNA_samp1 <- samples[which(samples$qPCR_plate_ID == '1' & samples$Molecule == 'cDNA'),]
cDNA_samp2 <- samples[which(samples$qPCR_plate_ID == '2' & samples$Molecule == 'cDNA'),]
cDNA_samp3 <- samples[which(samples$qPCR_plate_ID == '3' & samples$Molecule == 'cDNA'),]

# merge the sample and the Cq data

DNA_merge1 <- merge(DNA_samp1, DNA_plate1[, c("Well Position", "CT")], by.x = "qPCR_Well_ID", by.y = "Well Position")
DNA_amp1 <- merge(DNA_samp1, DNA_plate1_amp[, c("Well_ID", "Cycle", "Rn", "Delta_Rn")], by.x = "qPCR_Well_ID", by.y = "Well_ID")

DNA_merge2 <- merge(DNA_samp2, DNA_plate2[, c("Well Position", "CT")], by.x = "qPCR_Well_ID", by.y = "Well Position")
DNA_amp2 <- merge(DNA_samp2, DNA_plate2_amp[, c("Well_ID", "Cycle", "Rn", "Delta_Rn")], by.x = "qPCR_Well_ID", by.y = "Well_ID")

DNA_merge3 <- merge(DNA_samp3, DNA_plate3[, c("Well Position", "CT")], by.x = "qPCR_Well_ID", by.y = "Well Position")
DNA_amp3 <- merge(DNA_samp3, DNA_plate3_amp[, c("Well_ID", "Cycle", "Rn", "Delta_Rn")], by.x = "qPCR_Well_ID", by.y = "Well_ID")

cDNA_merge1 <- merge(cDNA_samp1, cDNA_plate1[, c("Well Position", "CT")], by.x = "qPCR_Well_ID", by.y = "Well Position")
cDNA_amp1 <- merge(cDNA_samp1, cDNA_plate1_amp[, c("Well_ID", "Cycle", "Rn", "Delta_Rn")], by.x = "qPCR_Well_ID", by.y = "Well_ID")

cDNA_merge2 <- merge(cDNA_samp2, cDNA_plate2[, c("Well Position", "CT")], by.x = "qPCR_Well_ID", by.y = "Well Position")
cDNA_amp2 <- merge(cDNA_samp2, cDNA_plate2_amp[, c("Well_ID", "Cycle", "Rn", "Delta_Rn")], by.x = "qPCR_Well_ID", by.y = "Well_ID")

cDNA_merge3 <- merge(cDNA_samp3, cDNA_plate3[, c("Well Position", "CT")], by.x = "qPCR_Well_ID", by.y = "Well Position")
cDNA_amp3 <- merge(cDNA_samp3, cDNA_plate3_amp[, c("Well_ID", "Cycle", "Rn", "Delta_Rn")], by.x = "qPCR_Well_ID", by.y = "Well_ID")


# for some reason, the CT column isn't read in as a number, so convert to numeric in R
DNA_merge1$CT <- as.numeric(DNA_merge1$CT)
DNA_merge2$CT <- as.numeric(DNA_merge2$CT)
DNA_merge3$CT <- as.numeric(DNA_merge3$CT)
cDNA_merge1$CT <- as.numeric(cDNA_merge1$CT)
cDNA_merge2$CT <- as.numeric(cDNA_merge2$CT)
cDNA_merge3$CT <- as.numeric(cDNA_merge3$CT)

# convert NA value to 40
DNA_merge1[is.na(DNA_merge1$CT),]$CT <- 40
DNA_merge2[is.na(DNA_merge2$CT),]$CT <- 40
DNA_merge3[is.na(DNA_merge3$CT),]$CT <- 40
cDNA_merge1[is.na(cDNA_merge1$CT),]$CT <- 40
cDNA_merge2[is.na(cDNA_merge2$CT),]$CT <- 40
cDNA_merge3[is.na(cDNA_merge3$CT),]$CT <- 40

# generate CLas 16S rRNA gene standard curve information
DNA_stand1 <- DNA_merge1[which(DNA_merge1$Experiment_name == "Standard_curve"),]
DNA_stand1$gene_copies <- DNA_stand1$Concentration * (1/1000^3) * (1/660) * 6.023e+23 * (1/size[["CLas_16S"]]) * vol2[["CLas_16S"]]

DNAmod1 <- summary(lm(CT~log10(gene_copies), data=DNA_stand1[c(3,4,5,6,7,8,9,10),]))
slopeDNA1 <- DNAmod1$coefficients[2,1]
interceptDNA1 <- DNAmod1$coefficients[1,1]
efficiencyDNA1 <- (-1 + 10^(-1/slopeDNA1))*100.0

DNA_stand2 <- DNA_merge2[which(DNA_merge2$Experiment_name == "Standard_curve"),]
DNA_stand2$gene_copies <- DNA_stand2$Concentration * (1/1000^3) * (1/660) * 6.023e+23 * (1/size[["CLas_16S"]]) * vol2[["CLas_16S"]]

DNAmod2 <- summary(lm(CT~log10(gene_copies), data=DNA_stand2[c(3,4,5,6,7,8,9,10),]))
slopeDNA2 <- DNAmod2$coefficients[2,1]
interceptDNA2 <- DNAmod2$coefficients[1,1]
efficiencyDNA2 <- (-1 + 10^(-1/slopeDNA2))*100.0


DNA_stand3 <- DNA_merge3[which(DNA_merge3$Experiment_name == "Standard_curve"),]
DNA_stand3$gene_copies <- DNA_stand3$Concentration * (1/1000^3) * (1/660) * 6.023e+23 * (1/size[["CLas_16S"]]) * vol2[["CLas_16S"]]

DNAmod3 <- summary(lm(CT~log10(gene_copies), data=DNA_stand3[c(3,4,5,6,7,8,9,10),]))
slopeDNA3 <- DNAmod3$coefficients[2,1]
interceptDNA3 <- DNAmod3$coefficients[1,1]
efficiencyDNA3 <- (-1 + 10^(-1/slopeDNA3))*100.0


cDNA_stand1 <- cDNA_merge1[which(cDNA_merge1$Experiment_name == "Standard_curve"),]
cDNA_stand1$gene_copies <- cDNA_stand1$Concentration * (1/1000^3) * (1/660) * 6.023e+23 * (1/size[["CLas_16S"]]) * vol2[["CLas_16S"]]

cDNAmod1 <- summary(lm(CT~log10(gene_copies), data=cDNA_stand1[c(3,4,5,6,7,8,9,10),]))
slopecDNA1 <- cDNAmod1$coefficients[2,1]
interceptcDNA1 <- cDNAmod1$coefficients[1,1]
efficiencycDNA1 <- (-1 + 10^(-1/slopecDNA1))*100.0

cDNA_stand2 <- cDNA_merge2[which(cDNA_merge2$Experiment_name == "Standard_curve"),]
cDNA_stand2$gene_copies <- cDNA_stand2$Concentration * (1/1000^3) * (1/660) * 6.023e+23 * (1/size[["CLas_16S"]]) * vol2[["CLas_16S"]]

cDNAmod2 <- summary(lm(CT~log10(gene_copies), data=cDNA_stand2[c(3,4,5,6,7,8,9,10),]))
slopecDNA2 <- cDNAmod2$coefficients[2,1]
interceptcDNA2 <- cDNAmod2$coefficients[1,1]
efficiencycDNA2 <- (-1 + 10^(-1/slopecDNA2))*100.0


cDNA_stand3 <- cDNA_merge3[which(cDNA_merge3$Experiment_name == "Standard_curve"),]
cDNA_stand3$gene_copies <- cDNA_stand3$Concentration * (1/1000^3) * (1/660) * 6.023e+23 * (1/size[["CLas_16S"]]) * vol2[["CLas_16S"]]

cDNAmod3 <- summary(lm(CT~log10(gene_copies), data=cDNA_stand3[c(3,4,5,6,7,8,9,10),]))
slopecDNA3 <- cDNAmod3$coefficients[2,1]
interceptcDNA3 <- cDNAmod3$coefficients[1,1]
efficiencycDNA3 <- (-1 + 10^(-1/slopecDNA3))*100.0

# create data.frame of qPCR assay values and write to file
stdcurve <- data.frame(slope = c(slopecDNA1, slopecDNA2, slopecDNA3, slopeDNA1, slopeDNA2, slopeDNA3),
                       intercept = c(interceptcDNA1, interceptcDNA2, interceptcDNA3, interceptDNA1, interceptDNA2, interceptDNA3),
                       efficiency = c(efficiencycDNA1, efficiencycDNA2, efficiencycDNA3, efficiencyDNA1, efficiencyDNA2, efficiencyDNA3),
                       Rsquare = c(cDNAmod1$r.squared, cDNAmod2$r.squared, cDNAmod3$r.squared, DNAmod1$r.squared, DNAmod2$r.squared, DNAmod3$r.squared),
                       molecule = c('cDNA', 'cDNA', 'cDNA', 'DNA', 'DNA', 'DNA'),
                       assay_number = c('cDNA1', 'cDNA2', 'cDNA3', 'DNA1', 'DNA2', 'DNA3'))

write.table(stdcurve, file = file.path(getwd(), "qPCR_stdcurve.txt"), quote = F, sep='\t', row.names=F, col.names = T)


# now add CLas 16S gene copy calculations to each qPCR data plate using the standard curve information for that
# respective plate
DNA_merge1$gene_copies <- rep(0, nrow(DNA_merge1))
DNA_merge2$gene_copies <- rep(0, nrow(DNA_merge2))
DNA_merge3$gene_copies <- rep(0, nrow(DNA_merge3))
cDNA_merge1$gene_copies <- rep(0, nrow(cDNA_merge1))
cDNA_merge2$gene_copies <- rep(0, nrow(cDNA_merge2))
cDNA_merge3$gene_copies <- rep(0, nrow(cDNA_merge3))

DNA_merge1$gene_copies <- 10^((DNA_merge1$CT - interceptDNA1)/slopeDNA1)
DNA_merge2$gene_copies <- 10^((DNA_merge2$CT - interceptDNA2)/slopeDNA2)
DNA_merge3$gene_copies <- 10^((DNA_merge3$CT - interceptDNA3)/slopeDNA3)

cDNA_merge1$gene_copies <- 10^((cDNA_merge1$CT - interceptcDNA1)/slopecDNA1)
cDNA_merge2$gene_copies <- 10^((cDNA_merge2$CT - interceptcDNA2)/slopecDNA2)
cDNA_merge3$gene_copies <- 10^((cDNA_merge3$CT - interceptcDNA3)/slopecDNA3)

# now combine all the data

totdatDNA <- rbind(DNA_merge1, DNA_merge2)
totdatDNA <- rbind(totdatDNA, DNA_merge3)

totdatcDNA <- rbind(cDNA_merge1, cDNA_merge2)
totdatcDNA <- rbind(totdatcDNA, cDNA_merge3)

totdatDNA$Reaction_volume <- rep(1, nrow(totdatDNA))
totdatcDNA$Reaction_volume <- rep(1, nrow(totdatcDNA))

totdatDNA$gene_copies_per_ng <- totdatDNA$gene_copies/(totdatDNA$Reaction_volume * totdatDNA$Concentration)
totdatcDNA$DNAse_concentration <- (totdatcDNA$Concentration * 17)/20
totdatcDNA$gene_copies_per_ng <- totdatcDNA$gene_copies/((totdatcDNA$DNAse_concentration * 7.5)/10)

totdatDNA$gene_copies_per_ng_cDNA <- totdatcDNA$gene_copies_per_ng


totdatDNA$gene_copies_per_extract <- (totdatDNA$gene_copies * (totdatDNA$Reaction_volume * totdatDNA$Concentration )) * (totdatDNA$Concentration * 70)
totdatcDNA$gene_copies_per_extract <- totdatcDNA$gene_copies * (((totdatcDNA$DNAse_concentration * 7.5)/10) * 10)

totdatDNA$log_gene_copies_per_ng <- log10(totdatDNA$gene_copies_per_ng)
totdatcDNA$log_gene_copies_per_ng <- log10(totdatcDNA$gene_copies_per_ng)

totdatDNA$log_gene_copies_per_extract <- log10(totdatDNA$gene_copies_per_extract)
totdatcDNA$log_gene_copies_per_extract <- log10(totdatcDNA$gene_copies_per_extract)

totdatDNA$cells_per_extract <- totdatDNA$gene_copies_per_extract/3
totdatDNA$log_cells_per_extract <- log10(totdatDNA$cells_per_extract)

totdatcDNA$cells_per_extract <- totdatcDNA$gene_copies_per_extract/3
totdatcDNA$log_cells_per_extract <- log10(totdatcDNA$cells_per_extract)

totdatcDNAB <- totdatcDNA[match(totdatDNA$Sample_ID, totdatcDNA$Sample_ID),]

totdatDNA$CT_cDNA <- totdatcDNAB$CT
totdatDNA$gene_copies_cDNA <- totdatcDNAB$gene_copies
totdatDNA$gene_copies_per_ng_cDNA <- totdatcDNAB$gene_copies_per_ng
totdatDNA$log_gene_copies_per_ng_cDNA <- totdatcDNAB$log_gene_copies_per_ng
totdatDNA$gene_copies_per_extract_cDNA <- totdatcDNAB$gene_copies_per_extract
totdatDNA$log_gene_copies_per_extract_cDNA <- totdatcDNAB$log_gene_copies_per_extract
totdatDNA$cells_per_extract_cDNA <- totdatcDNAB$cells_per_extract
totdatDNA$log_cells_per_extract_cDNA <- totdatcDNAB$log_cells_per_extract


detachleafDNA <- totdatDNA %>% filter(Experiment_name == 'Detached_leaf')
acquisitionDNA <- totdatDNA %>% filter(Experiment_name == 'Detached_leaf_acquisition')

detachleafcDNA <- totdatcDNA %>% filter(Experiment_name == 'Detached_leaf')
acquisitioncDNA <- totdatcDNA %>% filter(Experiment_name == 'Detached_leaf_acquisition')

# detached leaf assay with six NCR peptides

detachleafDNA$Day <- as.factor(detachleafDNA$Day)
detachleafcDNA$Day <- as.factor(detachleafcDNA$Day)

# add a "Trial" variable to enable faceting by detached leaf experiment
detachleafDNA$Trial <- rep('3', nrow(detachleafDNA))

# import other NCR peptide detached leaf assays
NCRpep12 <- readRDS(file.path(getwd(), "Detached_leaf_assay_data_inputs", "detach_leaf_NCR_peptides_1_2.rds"))

NCRpep38 <- readRDS(file.path(getwd(), "Detached_leaf_assay_data_inputs","detach_leaf_NCR_peptides_3-8.rds"))

NCRpep12a <- NCRpep12[,c("Sample_ID", "Concentration", "Ratio_260_280", "Ratio_260_230", "Treatment", "Treatment_name", "CT", "Reaction_volume", "gene_copies", "gene_copies_per_ng", "gene_copies_per_extract", "log_gene_copies_per_ng", "log_gene_copies_per_extract", "cells_per_extract", "log_cells_per_extract", "CT_cDNA", "gene_copies_cDNA", "gene_copies_per_ng_cDNA", "gene_copies_per_extract_cDNA", "log_gene_copies_per_extract_cDNA", "Day")]
colnames(NCRpep12a) <- c("Sample_ID", "Concentration", "Ratio_260_280", "Ratio_260_230", "Treatment_ID", "Treatment_name", "CT", "Reaction_volume", "gene_copies", "gene_copies_per_ng", "gene_copies_per_extract", "log_gene_copies_per_ng", "log_gene_copies_per_extract", "cells_per_extract", "log_cells_per_extract", "CT_cDNA", "gene_copies_cDNA", "gene_copies_per_ng_cDNA", "gene_copies_per_extract_cDNA", "log_gene_copies_per_extract_cDNA", "Day")

NCRpep12a$Treatment_ID <- NCRpep12a$Treatment_name
NCRpep12a[which(NCRpep12a$Treatment_ID == "KPO4 buffer"),]$Treatment_ID <- "KPO4_buffer"
NCRpep12a[which(NCRpep12a$Treatment_ID == "Polymyxin B sulfate"),]$Treatment_ID <- "PMB"


NCRpep38$log_gene_copies_per_ng <- log10(NCRpep38$gene_copies_per_ng)
NCRpep38$log_gene_copies_per_ng_cDNA <- log10(NCRpep38$gene_copies_per_ng_cDNA)
NCRpep38a <- NCRpep38[,c("Sample_ID", "Concentration", "Ratio_260_280", "Ratio_260_230", "Treatment", "Treatment_name", "CT", "Reaction_volume", "gene_copies", "gene_copies_per_ng", "gene_copies_per_extract", "log_gene_copies_per_ng", "log_gene_copies_per_extract", "cells_per_extract", "log_cells_per_extract", "CT_cDNA", "gene_copies_cDNA", "gene_copies_per_ng_cDNA", "gene_copies_per_extract_cDNA", "log_gene_copies_per_extract_cDNA", "Day")]
colnames(NCRpep38a) <- c("Sample_ID", "Concentration", "Ratio_260_280", "Ratio_260_230", "Treatment_ID", "Treatment_name", "CT", "Reaction_volume", "gene_copies", "gene_copies_per_ng", "gene_copies_per_extract", "log_gene_copies_per_ng", "log_gene_copies_per_extract", "cells_per_extract", "log_cells_per_extract", "CT_cDNA", "gene_copies_cDNA", "gene_copies_per_ng_cDNA", "gene_copies_per_extract_cDNA", "log_gene_copies_per_extract_cDNA", "Day")
NCRpep38a$Treatment_ID <- NCRpep38a$Treatment_name
NCRpep38a[which(NCRpep38a$Treatment_ID == "KPO4 buffer"),]$Treatment_ID <- "KPO4_buffer"
NCRpep38a[which(NCRpep38a$Treatment_ID == "Polymyxin B sulfate"),]$Treatment_ID <- "PMB"


NCRpep12a$Trial <- rep("1", nrow(NCRpep12a))
NCRpep38a$Trial <- rep("2", nrow(NCRpep38a))

NCRpep714 <- detachleafDNA[,c("Sample_ID", "Concentration", "Ratio_260_280", "Ratio_260_230", "Treatment_ID", "Treatment_name", "CT", "Reaction_volume", "gene_copies", "gene_copies_per_ng", "gene_copies_per_extract", "log_gene_copies_per_ng", "log_gene_copies_per_extract", "cells_per_extract", "log_cells_per_extract", "CT_cDNA", "gene_copies_cDNA", "gene_copies_per_ng_cDNA", "gene_copies_per_extract_cDNA", "log_gene_copies_per_extract_cDNA", "Day", "Trial")]

totNCRpep <- rbind(NCRpep12a, NCRpep38a)
totNCRpep <- rbind(totNCRpep, NCRpep714)
totNCRpep$Sample <- unlist(lapply(totNCRpep$Sample_ID, FUN=function(x){str_split(x, "_")[[1]][1]}))
totNCRpep[which(totNCRpep$Treatment_ID == "838921"),]$Treatment_ID <- "803573"

totNCRpep[c("277"),]$Sample_ID <- "B10_T7"
totNCRpep[c("278"),]$Sample_ID <- "B10_T7"
totNCRpep[c("277"),]$Sample <- "B10"
totNCRpep[c("278"),]$Sample <- "B10"
totNCRpep[c("301"),]$Sample_ID <- "A8_T7"
totNCRpep[c("302"),]$Sample_ID <- "A8_T7"
totNCRpep[c("301"),]$Sample <- "A8"
totNCRpep[c("302"),]$Sample <- "A8"
totNCRpep[c("301"),]$Day <- "7"
totNCRpep[c("302"),]$Day <- "7"
totNCRpep[c("313"),]$Sample_ID <- "A2_T7"
totNCRpep[c("314"),]$Sample_ID <- "A2_T7"
totNCRpep[c("313"),]$Sample <- "A2"
totNCRpep[c("314"),]$Sample <- "A2"

# perform ttest at 95% confidence interval threshold on rDNA

t1DNA <- totNCRpep[which(totNCRpep$Trial == "1" & totNCRpep$Sample %in% totNCRpep[which(totNCRpep$Day==0 & totNCRpep$CT <= 35),]$Sample),]
t1DNA$Treatment_ID <- as.factor(t1DNA$Treatment_ID)

treats <- unique(t1DNA$Treatment_ID)
welch95t1DNA <- data.frame()
for(t in treats){
  
  samps <- intersect(t1DNA[which(t1DNA$Treatment_ID == t & t1DNA$Day == "0"),]$Sample, t1DNA[which(t1DNA$Treatment_ID == t & t1DNA$Day == "7"),]$Sample)
  print(samps)
  s1 <- t1DNA %>% filter(Treatment_ID == t, Day == '0', Sample %in% samps) %>% arrange(Sample) %>% select(log_gene_copies_per_ng)
  s2 <- t1DNA %>% filter(Treatment_ID == t, Day == '7', Sample %in% samps) %>% arrange(Sample) %>% select(log_gene_copies_per_ng)
  print(dim(s1))
  print(dim(s2))
  ttestout <- t.test(s1[,1], s2[,1], alternative = "two.sided", paired = T, var.equal = F, conf.level = 0.95)
  welch95t1DNA <- rbind(welch95t1DNA, c(t, ttestout$estimate, ttestout$statistic, ttestout$p.value, ttestout$parameter, ttestout$conf.int[1], ttestout$conf.int[2]))
}

colnames(welch95t1DNA) <- c('Treatment', 'mean_difference', 't-value', 'p-value', 'df', 'lwr','upr')


trial1DNA <- ggplot(totNCRpep[which(totNCRpep$Trial == "1" & totNCRpep$Sample %in% totNCRpep[which(totNCRpep$Day==0 & totNCRpep$CT <= 35),]$Sample),], aes(x = Day, y = gene_copies_per_ng, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  ggpubr::stat_compare_means( method = "wilcox.test", 
                              label.x = 1.0, label.y = 4, size=4) +
  ylab("CLas 16S rDNA gene copies per ng nucleic acid") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

t2DNA <- totNCRpep[which(totNCRpep$Trial == "2" & totNCRpep$Sample %in% totNCRpep[which(totNCRpep$Day==0 & totNCRpep$CT <= 35),]$Sample),]
t2DNA$Treatment_ID <- as.factor(t2DNA$Treatment_ID)

treats <- unique(t2DNA$Treatment_ID)
welch95t2DNA <- data.frame()
for(t in treats){
  
  samps <- intersect(t2DNA[which(t2DNA$Treatment_ID == t & t2DNA$Day == "0"),]$Sample, t2DNA[which(t2DNA$Treatment_ID == t & t2DNA$Day == "7"),]$Sample)
  print(samps)
  s1 <- t2DNA %>% filter(Treatment_ID == t, Day == '0', Sample %in% samps) %>% arrange(Sample) %>% select(log_gene_copies_per_ng)
  s2 <- t2DNA %>% filter(Treatment_ID == t, Day == '7', Sample %in% samps) %>% arrange(Sample) %>% select(log_gene_copies_per_ng)
  print(dim(s1))
  print(dim(s2))
  ttestout <- t.test(s1[,1], s2[,1], alternative = "two.sided", paired = T, var.equal = F, conf.level = 0.95)
  welch95t2DNA <- rbind(welch95t2DNA, c(t, ttestout$estimate, ttestout$statistic, ttestout$p.value, ttestout$parameter, ttestout$conf.int[1], ttestout$conf.int[2]))
}

colnames(welch95t2DNA) <- c('Treatment', 'mean_difference', 't-value', 'p-value', 'df', 'lwr','upr')


trial2DNA <- ggplot(totNCRpep[which(totNCRpep$Trial == "2" & totNCRpep$Sample %in% totNCRpep[which(totNCRpep$Day==0 & totNCRpep$CT <= 35),]$Sample),], aes(x = Day, y = gene_copies_per_ng, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  ggpubr::stat_compare_means( method = "wilcox.test", 
                              label.x = 1.0, label.y = 6, size=4) +
  ylab("CLas 16S rDNA gene copies per ng nucleic acid") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

trial3DNA <- ggplot(totNCRpep[which(totNCRpep$Trial == "3" & totNCRpep$Sample %in% totNCRpep[which(totNCRpep$Day==0 & totNCRpep$CT <= 35),]$Sample),], aes(x = Day, y = gene_copies_per_ng, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  ggpubr::stat_compare_means( method = "wilcox.test", 
                              label.x = 1.0, label.y = 4, size=4) +
  ylab("CLas 16S rDNA gene copies per ng nucleic acid") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

# cDNA

trial1cDNA <- ggplot(totNCRpep[which(totNCRpep$Trial == "1" & totNCRpep$Sample %in% totNCRpep[which(totNCRpep$Day==0 & totNCRpep$CT_cDNA <= 35),]$Sample),], aes(x = Day, y = gene_copies_per_ng_cDNA, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  ggpubr::stat_compare_means( method = "wilcox.test", 
                              label.x = 1.0, label.y = 4, size=4) +
  ylab("CLas 16S rRNA gene copies per ng nucleic acid") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

trial2cDNA <- ggplot(totNCRpep[which(totNCRpep$Trial == "2" & totNCRpep$Sample %in% totNCRpep[which(totNCRpep$Day==0 & totNCRpep$CT_cDNA <= 40),]$Sample),], aes(x = Day, y = gene_copies_per_ng_cDNA, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  ggpubr::stat_compare_means( method = "wilcox.test", 
                              label.x = 1.0, label.y = 6, size=4) +
  ylab("CLas 16S rRNA gene copies per ng nucleic acid") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

trial3cDNA <- ggplot(totNCRpep[which(totNCRpep$Trial == "3"),], aes(x = Day, y = gene_copies_per_ng_cDNA, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  ggpubr::stat_compare_means( method = "wilcox.test", 
                              label.x = 1.0, label.y = 4, size=4) +
  ylab("CLas 16S rRNA gene copies per ng nucleic acid") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

# now compare 16S rRNA to 16S rDNA ratio

trial1ratio <- ggplot(totNCRpep[which(totNCRpep$Trial == "1"),], aes(x = Day, y = gene_copies_per_ng_cDNA/gene_copies_per_ng, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  theme(axis.title = element_blank(), strip.text.x = element_text(size = 6)) +
  ggpubr::stat_compare_means( aes(label = ..p.signif..), method = "wilcox.test", 
                              label.x = 1.5, label.y = 4, size=5) +
  #ylab("Ratio of CLas 16S rRNA to rDNA gene copies") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

tmpncrpep <- totNCRpep
remrows <- which(tmpncrpep$Sample %in% c('6B', '6C') & tmpncrpep$Treatment_ID == 'KPO4_buffer' & tmpncrpep$Trial == '2')
totNCRpep <- totNCRpep[-remrows,]

trial2ratio <- ggplot(totNCRpep[which(totNCRpep$Trial == "2"),], aes(x = Day, y = gene_copies_per_ng_cDNA/gene_copies_per_ng, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  theme(axis.title = element_blank(), strip.text.x = element_text(size = 6)) +
  ggpubr::stat_compare_means( aes(label = ..p.signif..), method = "wilcox.test", 
                              label.x = 1.5, label.y = 3, size=5) +
  #ylab("Ratio of CLas 16S rRNA to rDNA gene copies") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

tmpncrpep <- totNCRpep
remrows <- which(tmpncrpep$Sample %in% c('W2', 'W3', 'W5', 'W6') & tmpncrpep$Treatment_ID == 'KPO4_buffer' & tmpncrpep$Trial == '3')
totNCRpep <- totNCRpep[-remrows,]


trial3ratio <- ggplot(totNCRpep[which(totNCRpep$Trial == "3" ),], aes(x = Day, y = gene_copies_per_ng_cDNA/gene_copies_per_ng, color = Day, group = Day)) +     
  geom_jitter(width = 0.2, alpha = 0.5) +
  stat_summary(fun.data = mean_cl_boot, alpha = 1, shape = 3, size=1) +
  facet_grid(Trial~Treatment_ID) + 
  scale_y_log10() +
  theme_bw() + 
  th + 
  theme(axis.title = element_blank(), strip.text.x = element_text(size = 6)) +
  ggpubr::stat_compare_means( aes(label = ..p.signif..), method = "wilcox.test", 
                              label.x = 1.5, label.y = 3, size=5) +
  #ylab("Ratio of CLas 16S rRNA to rDNA gene copies") +
  scale_color_manual(values = c("#004949", "#ff6db6"))

# create plot
windowsFonts("Arial" = windowsFont("Arial"))
y.grob <- textGrob("Ratio of CLas 16S rRNA to rDNA gene copies", 
                   gp=gpar(fontface="bold", col="Black", family = 'Arial', fontsize=16), rot=90)

x.grob <- textGrob("Day", 
                   gp=gpar(fontface="bold", col="Black", family = 'Arial', fontsize=16))

#add to plot

trialcomboplot <- plot_grid(trial1ratio, trial2ratio, trial3ratio, labels = 'AUTO', nrow=3)
figure3 <- arrangeGrob(arrangeGrob(trialcomboplot, left = y.grob, bottom = x.grob))

#############
#  FIGURE 3 #
#############

ggsave(figure3, file = file.path(getwd(), "figures", "Figure3.pdf"), dpi = 300, width = 8, height = 6, device = cairo_pdf)
ggsave(figure3, file = file.path(getwd(), "figures", "Figure3.eps"), dpi = 300, width = 8, height = 6)
ggsave(figure3, file = file.path(getwd(), "figures", "Figure3.png"), dpi = 300, width = 8, height = 6)
ggsave(figure3, file = file.path(getwd(), "figures", "Figure3.svg"), dpi = 300, width = 8, height = 6)
ggsave(figure3, file = file.path(getwd(), "figures", "Figure3.tiff"), dpi = 300, width = 8, height = 6)

#############
# FIGURE S2 #
#############

growth.pca3 <- prcomp(growthdf4c, scale. = T)
growth.pca3.scores <- scores(growth.pca3)
growthdf6 <- growthdf4a
growthdf6$PC1 <- growth.pca3.scores[,1]
growthdf6$PC2 <- growth.pca3.scores[,2]
growthdf6$PC3 <- growth.pca3.scores[,3]

cor.test(growthdf4c[,2], growthdf4c[,1], method = "spearman")
cor.test(growthdf4c[,3], growthdf4c[,1], method = "spearman")
cor.test(growthdf4c[,4], growthdf4c[,1], method = "spearman")

ggbiplot(growth.pca2, scale = F, alpha=0) + geom_point(data = growthdf5, aes(x = PC1, y = PC2, color = percent_mumax), size=2) + theme_light() + theme(strip.background = element_rect(fill = 'black')) + guides(color = guide_legend("Growth rate inhibition (%)"), shape = guide_legend("Setup")) + scale_color_viridis_c(direction=-1) + theme_light() 

pc1 <- ggbiplot(growth.pca3, scale = F, alpha=0) + guides(color = F) + geom_point(data = growthdf6, aes(x = PC1, y = PC2, color = percent_mumax), size=2, alpha=0.3) + theme_light() + theme(strip.background = element_rect(fill = 'black'), legend.position = 'none') + scale_color_viridis_c(direction=-1) + theme_light() + ylim(c(-3,3))
pc1a <- ggbiplot(growth.pca3, scale = F, alpha=1, group = growthdf6$percent_mumax)  + theme_light() + theme(legend.direction = 'horizontal', legend.position = 'top', axis.text = element_text(color = 'black', size = 10, family = 'Arial', face = 'bold'), axis.title = element_text(color = 'black', size = 16, family = 'Arial', face = 'bold')) + scale_color_distiller(palette = "PuBu", direction = 1, limits = c(0, 75)) + ylim(c(-5.5,5.5)) + xlim(c(-7, 7)) + guides(color = guide_colourbar(label.position = "bottom", title = "Growth inhibition (%)", title.position = 'top'))
pc2a <- ggplot(growthdf6, aes(x = PC1, y = AMP.charge.at.pH.7)) + geom_point() + theme_light() + geom_smooth(method = 'lm', formula = 'y~x') + ylab("Predicted NCR peptide charge (pH = 7)")
pc2b <- ggplot(growthdf6, aes(x = Number.Cysteines, y = percent_mumax)) + geom_point() + theme_light() + geom_smooth(method = 'lm', formula = 'y~x') + ylab("Growth rate inhibition (%)") + xlab("No. cysteine residues")
pc3 <- ggplot(growthdf6, aes(x = GRAVY.score, y = percent_mumax)) + geom_point() + theme_light() + geom_smooth(method = 'lm', formula = 'y~x') + ylab("Growth rate inhibition (%)") + xlab("GRAVY score")
pc4 <- ggplot(growthdf6, aes(x = AMP.charge.at.pH.7, y = percent_mumax)) + geom_point() + theme_light() + geom_smooth(method = 'lm', formula = 'y~x') + ylab("Growth rate inhibition (%)")
pcplots <- plot_grid(pc1,pc2b,pc3,pc4)

ggsave(pcplots, file = file.path(getwd(), "figures", "FigureS2.pdf"), dpi = 300, width = 8, height = 6)
ggsave(pcplots, file = file.path(getwd(), "figures", "FigureS2.eps"), dpi = 300, width = 8, height = 6)
ggsave(pcplots, file = file.path(getwd(), "figures", "FigureS2.png"), dpi = 300, width = 8, height = 6)
ggsave(pcplots, file = file.path(getwd(), "figures", "FigureS2.svg"), dpi = 300, width = 8, height = 6)
ggsave(pcplots, file = file.path(getwd(), "figures", "FigureS2.tiff"), dpi = 300, width = 8, height = 6)

#############
# FIGURE S3 #
#############

# now create a resolution V fractional factorial design for the 10 most
# inhibitory NCR peptides

# get the path to the R input data files containing Liberibacter crescens strain BT-1 growth data
growthdata <- Sys.glob(file.path(getwd(), "NCR_Frac_Fac_128", "NCR_Peptides_*", "*_R_data_input.txt"))

# read in growth data from NCR peptide Fractional Factorial assay 96-well plate data
growthdf <- data.frame()
for(g in 1:length(growthdata)){
  tmpdata <- read.table(growthdata[[g]], sep="\t", header=T, stringsAsFactors = T)
  tmpdata$plate <- rep(g, dim(tmpdata)[1])
  growthdf <- rbind(growthdf, tmpdata)
}

# use growthrates package to estimate growth rates

growthdf2 <- growthdf %>% 
             group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% 
             summarise(mean_Abs = mean(Absorbance), sd_Abs = sd(Absorbance))

split.data <- multisplit(growthdf2, c("Sample_ID", "Type", "Concentration", "plate"))

# fit an plot fit_easylinear fits
logdat2 <- NULL
logdat2plots <- list()
for(i in 1:length(split.data)){
  print(i)
  tmpdat2 <- split.data[[i]]
  matchkey <- unique(as.vector(tmpdat2)[2])
  if(matchkey == "PMB" | matchkey == "BM7_only" | matchkey == "Y4761"){
    #simpreg <- lm(tmpdat2$Absorbance~tmpdat2$Timepoint, data=tmpdat2)
    simpreg <- lm(tmpdat2$mean_Abs~tmpdat2$Timepoint, data=tmpdat2)
    simppred <- as.data.frame(predict(simpreg))
    colnames(simppred) <- c("Absorbance.pred")
    tmpdat2$y0 <- rep(0.2, dim(tmpdat2)[1])
    tmpdat2$y0_lm <- rep(0.01, dim(tmpdat2)[1])
    tmpdat2$mumax <- rep(0.01, dim(tmpdat2)[1])
    tmpdat2$lag <- rep(2, dim(tmpdat2)[1])
    tmpdat2$Absorbance.pred <- simppred$Absorbance.pred
    logdat2 <- rbind(logdat2, tmpdat2)
    #logdat2plots[[i]] <- ggplot(tmpdat2) + geom_point(aes(x = Timepoint, y = Absorbance)) + geom_line(aes(x = Timepoint, y = Absorbance.pred, group = Well_ID), linetype = 2, size = 1) + theme_light() + theme(plot.title = element_text(hjust=0.5), axis.text = element_text(color = 'black')) + scale_y_log10() + xlab("Time (days)") + ylab(expression('log'[10]~'Absorbance'['600nm'])) + ggtitle(unique(as.vector(tmpdat2)[1])[1]$Sample_ID)
    
  }else{
    #tmpfit <- fit_growthmodel(FUN = grow_logistic, p = p, tmpdat$Timepoint, tmpdat$Absorbance,
    #                          lower = lower, upper = upper)
    # can alternatively fit a heuristic method to generate data - this method also provides
    # a lag parameter which could be useful
    #tmpfit <- fit_easylinear(tmpdat2$Timepoint, tmpdat2$Absorbance, h = 5, quota = 0.95)
    print(matchkey)
    tmpfit <- fit_easylinear(tmpdat2$Timepoint, tmpdat2$mean_Abs, h = 5, quota = 0.95)
    tmpdat2$y0 <- rep(coef(tmpfit)[1], dim(tmpdat2)[1])
    tmpdat2$y0_lm <- rep(coef(tmpfit)[2], dim(tmpdat2)[1])
    tmpdat2$mumax <- rep(coef(tmpfit)[3], dim(tmpdat2)[1])
    tmpdat2$lag <- rep(coef(tmpfit)[4], dim(tmpdat2)[1])
    tmpdat2$Absorbance.pred <- predict(tmpfit)$y
    logdat2 <- rbind(logdat2, tmpdat2)
    #logdat2plots[[i]] <- ggplot(tmpdat2) + geom_point(aes(x = Timepoint, y = Absorbance)) + geom_line(aes(x = Timepoint, y = Absorbance.pred, group = Well_ID), linetype = 2, size = 1) + theme_light() + theme(plot.title = element_text(hjust=0.5), axis.text = element_text(color = 'black')) + scale_y_log10() + xlab("Time (days)") + ylab(expression('log'[10]~'Absorbance'['600nm'])) + ggtitle(unique(as.vector(tmpdat2)[1])[1]$Sample_ID)
    
  }
}

# for fit_easylinear, check to see which peptides were most inhibitory
#logdat2a <- logdat2 %>% group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% summarise(mean_Abs = mean(Absorbance), sd_Abs = sd(Absorbance), mean_y0 = mean(y0), sd_y0 = sd(y0), meany0_lm = mean(y0_lm), sd_y0_lm = sd(y0_lm), mean_mumax = mean(mumax), sd_mumax = sd(mumax), mean_lag = mean(lag), sd_lag = sd(lag))
logdat2a <- logdat2 %>% group_by(plate, Sample_ID, Type, Concentration, Timepoint) %>% summarise(mean_Abs = mean(mean_Abs), sd_Abs = sd(mean_Abs), mean_y0 = mean(y0), sd_y0 = sd(y0), meany0_lm = mean(y0_lm), sd_y0_lm = sd(y0_lm), mean_mumax = mean(mumax), sd_mumax = sd(mumax), mean_lag = mean(lag), sd_lag = sd(lag))
logdat4a <- logdat2a %>% group_by(plate) %>% mutate(norm_mumax = mean_mumax[Type == "Negative_Control"] - mean_mumax, percent_mumax = (1 - (mean_mumax/mean_mumax[Type == "Negative_Control"]))*100.0)
logdat4a <- logdat4a %>% select(Sample_ID, Concentration, percent_mumax) %>% distinct()

# now create a resolution V fractional factorial design for the 10 most
# inhibitory NCR peptides

# read in the previously identified top inhibitor NCR peptides from file
inhibiNCRdf <- read.table(file = file.path(getwd(), "Lcre_growth_assay_182_NCR_peptides", "Top_47_NCR_peptide_inhibitors.txt"), sep = ",", header=T)

da<- FrF2(nfactors = 10, resolution = 5, default.levels = c('-1','1'), factor.names = as.vector(as.character(inhibiNCRdf$Sample_ID[1:10])), randomize = F, method = "VF2")

# extract 1 mg/ml NCR peptide combinations
onemg <- logdat4a[which(logdat4a$Concentration == "1_mg_ml"),]
# order the data frame by the combination used in the fractional factorial assay
onemgb <- onemg[stringr::str_order(onemg$Sample_ID, numeric = T),]

ponemg <- logdat4a[which(logdat4a$Concentration == "0.1_mg_ml"),]
ponemgb <- ponemg[stringr::str_order(ponemg$Sample_ID, numeric = T),]

planonemg <- add.response(da, onemgb$percent_mumax)

planponemg <- add.response(da, ponemgb$percent_mumax)


summary(lm(planonemg))
summary(aov(planonemg))


interlm <- lm(planonemg)
pr <- residuals(interlm)/(1 - lm.influence(interlm)$hat)
PRESS <- sum(pr^2)
my.anova <- anova(interlm)
tss <- sum(my.anova$"Sum Sq")
# predictive R^2
pred.r.squared <- 1 - PRESS/(tss)
pred.r.squared

summary(lm(planponemg))

write.table(planonemg, file = file.path(getwd(), "NCR_Frac_Fac_128", "NCR_Frac_Fac_Lcrescens_growth_inhib_1mgml.txt"), quote = F, sep='\t', row.names=F, col.names = T)
write.table(planponemg, file = file.path(getwd(), "NCR_Frac_Fac_128", "NCR_Frac_Fac_Lcrescens_growth_inhib_0.1mgml.txt"), quote = F, sep='\t', row.names=F, col.names = T)

tmpmep <- MEPlot(planonemg, abbrev = 5, cex.xax = 1.6, cex.main = 2)
tmpmep <- as.data.frame(t(tmpmep))
tmpmep$peptide <- gsub("X", "", rownames(tmpmep))
colnames(tmpmep) <- c('Absent', 'Present', 'Peptide')
tmpmep2 <- tmpmep %>% gather(Present, Growth_Inhibition, Absent:Present, factor_key=T)
tmpmep2$Peptide <- paste("NCR peptide ", unlist(lapply(strsplit(tmpmep2$Peptide, "\\."), FUN = function(x){return(x[1])})), sep = "")
mep1 <- ggplot(tmpmep2, aes(x = Present, y = Growth_Inhibition)) + 
        geom_point() + 
        geom_line(aes(group = Peptide)) + 
        theme_bw() + theme(axis.title.x = element_blank()) + 
        facet_wrap(~Peptide) + 
        ylab("Growth rate inhibition (%)")

ggsave(mep1, file = file.path(getwd(), "figures", "FigureS3.pdf"), dpi = 300, width = 8, height = 6)
ggsave(mep1, file = file.path(getwd(), "figures", "FigureS3.eps"), dpi = 300, width = 8, height = 6)
ggsave(mep1, file = file.path(getwd(), "figures", "FigureS3.png"), dpi = 300, width = 8, height = 6)
ggsave(mep1, file = file.path(getwd(), "figures", "FigureS3.svg"), dpi = 300, width = 8, height = 6)
ggsave(mep1, file = file.path(getwd(), "figures", "FigureS3.tiff"), dpi = 300, width = 8, height = 6)

#############
# FIGURE S4 #
#############

# now plot the main effects and interactions among peptides as a plot

pdf(file = file.path(getwd(), "figures", "FigureS4.pdf"))
op <- par(mar = c(7,7,4,2) + 0.1, mgp = c(3,0.5,0))
IAPlot(planonemg, abbrev = 5, show.alias = FALSE, lwd = 2, cex = 0.6, cex.xax = 0.6, cex.lab = 0.6, main = NULL)
title(main = NULL, ylab = 'Growth inhibition (%)', xlab = 'NCR peptide presence or absence', line = 6)
dev.off()
