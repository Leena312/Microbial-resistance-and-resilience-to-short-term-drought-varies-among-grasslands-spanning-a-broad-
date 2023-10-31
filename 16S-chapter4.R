##alpha diversity

library(mctoolsr)
library(vegan)
library(car)
library(metacoder)
library(emmeans)
library(multcompView)
library(dplyr)
library(tidyverse)

##All data
input <- load_taxa_table("~/Desktop/R/zotutab_16S.txt", "~/Desktop/R/16SMappingFile.txt")

colSums(input$data_loaded)
inputr = single_rarefy(input, 4500)
##Filter down
Control = filter_data(inputr, 'Treatment', keep_vals = 'Control')
ControlSGS = filter_data(Control, 'Site', keep_vals = 'SGS')
ControlCHY = filter_data(Control, 'Site', keep_vals = 'CHY')
ControlHYS = filter_data(Control, 'Site', keep_vals = 'HYS')
ControlKNZ = filter_data(Control, 'Site', keep_vals = 'KNZ')
Drought = filter_data(inputr, 'Treatment', keep_vals = 'Drought - Resistance')
DroughtSGS = filter_data(Drought, 'Site', keep_vals = 'SGS')
DroughtCHY = filter_data(Drought, 'Site', keep_vals = 'CHY')
DroughtHYS = filter_data(Drought, 'Site', keep_vals = 'HYS')
DroughtKNZ = filter_data(Drought, 'Site', keep_vals = 'KNZ')
##Shannon's Diversity
ControlSGSshan <- diversity(t(ControlSGS$data_loaded), index = "shannon")
ControlCHYshan <- diversity(t(ControlCHY$data_loaded), index = "shannon")
ControlHYSshan <- diversity(t(ControlHYS$data_loaded), index = "shannon")
ControlKNZshan <- diversity(t(ControlKNZ$data_loaded), index = "shannon")

DroughtSGSshan <- diversity(t(DroughtSGS$data_loaded), index = "shannon")
DroughtCHYshan <- diversity(t(DroughtCHY$data_loaded), index = "shannon")
DroughtHYSshan <- diversity(t(DroughtHYS$data_loaded), index = "shannon")
DroughtKNZshan <- diversity(t(DroughtKNZ$data_loaded), index = "shannon")

ControlSGSshan<-as.data.frame(ControlSGSshan)
ControlCHYshan<-as.data.frame(ControlCHYshan)
ControlHYSshan<-as.data.frame(ControlHYSshan)
ControlKNZshan<-as.data.frame(ControlKNZshan)
DroughtSGSshan<-as.data.frame(DroughtSGSshan)
DroughtCHYshan<-as.data.frame(DroughtCHYshan)
DroughtHYSshan<-as.data.frame(DroughtHYSshan)
DroughtKNZshan<-as.data.frame(DroughtKNZshan)
write.csv(ControlSGSshan,"ControlSGSshan.csv")
rename(ControlSGSshan,shan,ControlSGSshan)

b<-rbind(ControlSGSshan,ControlCHYshan,ControlHYSshan,ControlKNZshan,DroughtSGSshan,DroughtCHYshan,DroughtHYSshan,DroughtKNZshan)

colors<-c("blue","orange","blue","orange","blue","orange","blue","orange")

pdf(file="Shannon Diversity chapter 4 all 16S",
    width=10,
    height=8)
a<-boxplot(ControlSGSshan, DroughtSGSshan, ControlCHYshan, DroughtCHYshan, ControlHYSshan, DroughtHYSshan, ControlKNZshan, DroughtKNZshan, col=colors, main = "Treatment", ylim=c(6.2,7.2),ylab = "Shannon's Diversity", xlab = "Treament", 
           names = c("SGS-C", "SGS-D", "CHY-C","CHY-D","HYS-C","HYS-D","KNZ-C","KNZ-D"))
dev.off()

##STATS
controldrought = filter_data(inputr, 'Treatment', filter_vals = 'Drought - Recovery 1')
controldrought2 = filter_data(controldrought, 'Treatment', filter_vals = 'Drought - Recovery 2')
controldrought3 = filter_data(controldrought2, 'Treatment', filter_vals = 'Control - Legacy Effects')
controldrought4 = filter_data(controldrought3, 'Treatment', filter_vals = 'Drought - Resistance to Second Drought')
controldrought5 = filter_data(controldrought4, 'Treatment', filter_vals = 'Drought - Second drought recovery 1')
controldrought6 = filter_data(controldrought5, 'Treatment', filter_vals = 'Drought - Second drought recovery 2')

allshan <- diversity(t(controldrought6$data_loaded), index = "shannon")

allshan <- allshan[match(names(allshan), row.names(controldrought6$map_loaded))]
allshan_df <- cbind(controldrought6$map_loaded, allshan)

controldrought6$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment*Site, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")

emmeans<-emmeans(allshan_glm,pairwise~Treatment|Site)
eff_size(emmeans,sigma(allshan_glm),edf=3)

##Richness
Controrich<-specnumber(t(Control$data_loaded))
Droughtrich<-specnumber(t(Drought$data_loaded))
Controrich<-as.data.frame(Controrich)
Droughtrich<-as.data.frame(Droughtrich)

write.csv(Controrich,"Controrich.csv")
write.csv(Droughtrich,"Droughtrich.csv")

ControlSGSrich <- specnumber(t(ControlSGS$data_loaded))
ControlCHYrich <- specnumber(t(ControlCHY$data_loaded))
ControlHYSrich <- specnumber(t(ControlHYS$data_loaded))
ControlKNZrich <- specnumber(t(ControlKNZ$data_loaded))

DroughtSGSrich <- specnumber(t(DroughtSGS$data_loaded))
DroughtCHYrich <- specnumber(t(DroughtCHY$data_loaded))
DroughtHYSrich <- specnumber(t(DroughtHYS$data_loaded))
DroughtKNZrich <- specnumber(t(DroughtKNZ$data_loaded))

ControlSGSrich<-as.data.frame(ControlSGSrich)
ControlCHYrich<-as.data.frame(ControlCHYrich)
ControlHYSrich<-as.data.frame(ControlHYSrich)
ControlKNZrich<-as.data.frame(ControlKNZrich)
DroughtSGSrich<-as.data.frame(DroughtSGSrich)
DroughtCHYrich<-as.data.frame(DroughtCHYrich)
DroughtHYSrich<-as.data.frame(DroughtHYSrich)
DroughtKNZrich<-as.data.frame(DroughtKNZrich)
ControlSGSrich %>% rename(rich = ControlSGSrich)

b<-rbind(ControlSGSshan,ControlCHYshan,ControlHYSshan,ControlKNZshan,DroughtSGSshan,DroughtCHYshan,DroughtHYSshan,DroughtKNZshan)

write.csv(ControlSGSshan,"ControlSGSshan.csv")

colors<-c("blue","orange","blue","orange","blue","orange","blue","orange")

pdf(file="Richness chapter 4 16S",
    width=10,
    height=8)
boxplot(ControlSGSrich, DroughtSGSrich, ControlCHYrich, DroughtCHYrich, ControlHYSrich, DroughtHYSrich, ControlKNZrich, DroughtKNZrich, col=colors, main = "Treatment", ylab = "Richness", xlab = "Treament", 
        names = c("SGS-C", "SGS-D", "CHY-C","CHY-D","HYS-C","HYS-D","KNZ-C","KNZ-D"))
dev.off()

allshan <- specnumber(t(controldrought6$data_loaded))

allshan <- allshan[match(names(allshan), row.names(controldrought6$map_loaded))]
allshan_df <- cbind(controldrought6$map_loaded, allshan)

controldrought6$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment*Site, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")

emmeans<-emmeans(allshan_glm,pairwise~Treatment|Site)
eff_size(emmeans,sigma(allshan_glm),edf=3)




##CAPS


library(dplyr)
library(BiodiversityR)
library (vegan)
library (MASS)

## All data
input <- load_taxa_table("~/Desktop/R/zotutab_16S.txt", "~/Desktop/R/16SMappingFile.txt")

colSums(input$data_loaded)
inputr = single_rarefy(input, 4500)

controldrought = filter_data(inputr, 'Treatment', filter_vals = 'Drought - Recovery 1')
controldrought2 = filter_data(controldrought, 'Treatment', filter_vals = 'Drought - Recovery 2')
controldrought3 = filter_data(controldrought2, 'Treatment', filter_vals = 'Control - Legacy Effects')
controldrought4 = filter_data(controldrought3, 'Treatment', filter_vals = 'Drought - Resistance to Second Drought')
controldrought5 = filter_data(controldrought4, 'Treatment', filter_vals = 'Drought - Second drought recovery 1')
final = filter_data(controldrought5, 'Treatment', filter_vals = 'Drought - Second drought recovery 2')

SGS = filter_data(final, 'Site', keep_vals = 'SGS')
CHY = filter_data(final, 'Site', keep_vals = 'CHY')
HYS = filter_data(final, 'Site', keep_vals = 'HYS')
KNZ = filter_data(final, 'Site', keep_vals = 'KNZ')

#NMDS
library(ggplot2)
library(vegan)

SGS.mds<-metaMDS(t(SGS$data_loaded))
data.scoresSGS <- as.data.frame(scores(SGS.mds)) 
data.scoresSGS$Treatment <- SGS$map_loaded$Treatment

NMDS1SGS <- glm(NMDS2 ~ Treatment, data = data.scoresSGS)
Anova(NMDS1SGS, test.statistic = "F")

emmeans<-emmeans(NMDS1SGS,pairwise~Treatment)
eff_size(emmeans,sigma(NMDS1SGS),edf=3)

a<-ggplot(data=data.scoresSGS,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + 
  stat_ellipse(data=data.scoresSGS,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the convex hulls
  ggtitle("Stress = 0.068")+
  geom_point(data=data.scoresSGS,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the point markers
  scale_colour_manual(values=c("Control" = "blue", "Drought - Resistance" = "orange")) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=8),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))

CHY.mds<-metaMDS(t(CHY$data_loaded))
data.scoresCHY <- as.data.frame(scores(CHY.mds)) 
data.scoresCHY$Treatment <- CHY$map_loaded$Treatment

NMDS1SGS <- glm(NMDS2 ~ Treatment, data = data.scoresCHY)
Anova(NMDS1SGS, test.statistic = "F")

emmeans<-emmeans(NMDS1SGS,pairwise~Treatment)
eff_size(emmeans,sigma(NMDS1SGS),edf=3)

b<-ggplot(data=data.scoresCHY,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + 
  stat_ellipse(data=data.scoresCHY,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the convex hulls
  ggtitle("Stress = 0.076")+
  geom_point(data=data.scoresCHY,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the point markers
  scale_colour_manual(values=c("Control" = "blue", "Drought - Resistance" = "orange")) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=8),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    
    plot.title=element_text(size=20,face="bold",vjust=2))

HYS.mds<-metaMDS(t(HYS$data_loaded))
data.scoresHYS <- as.data.frame(scores(HYS.mds)) 
data.scoresHYS$Treatment <- HYS$map_loaded$Treatment

NMDS1SGS <- glm(NMDS2 ~ Treatment, data = data.scoresHYS)
Anova(NMDS1SGS, test.statistic = "F")

emmeans<-emmeans(NMDS1SGS,pairwise~Treatment)
eff_size(emmeans,sigma(NMDS1SGS),edf=3)

c<-ggplot(data=data.scoresHYS,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + 
  stat_ellipse(data=data.scoresHYS,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the convex hulls
  ggtitle("Stress = 0.096")+
  geom_point(data=data.scoresHYS,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the point markers
  scale_colour_manual(values=c("Control" = "blue", "Drought - Resistance" = "orange")) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=12),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))

KNZ.mds<-metaMDS(t(KNZ$data_loaded))
data.scoresKNZ <- as.data.frame(scores(KNZ.mds)) 
data.scoresKNZ$Treatment <- KNZ$map_loaded$Treatment

NMDS1SGS <- glm(NMDS2 ~ Treatment, data = data.scoresKNZ)
Anova(NMDS1SGS, test.statistic = "F")

emmeans<-emmeans(NMDS1SGS,pairwise~Treatment)
eff_size(emmeans,sigma(NMDS1SGS),edf=3)

d<-ggplot(data=data.scoresKNZ,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + 
  stat_ellipse(data=data.scoresKNZ,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the convex hulls
  ggtitle("Stress = 0.119")+
  geom_point(data=data.scoresKNZ,aes(x=NMDS1,y=NMDS2,group=factor(Treatment),color=factor(Treatment))) + # add the point markers
  scale_colour_manual(values=c("Control" = "blue", "Drought - Resistance" = "orange")) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=12),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))

a<-rbind(data.scoresSGS,data.scoresCHY,data.scoresHYS,data.scoresKNZ)
write.csv(a,"NMDS16S.csv")

library(patchwork)

final<-a+b+c+d

ggsave(filename = "16SNMDS.pdf",
       plot = final,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##CAPS

otuTAB<-KNZ$data_loaded

otuTAB <- as.data.frame (otuTAB)
otuTABtranspose<-t(otuTAB)
mta<-KNZ$map_loaded
mta$Treatment <- as.factor (mta$Treatment)

o.dist <-vegdist (otuTABtranspose)
Ordination.model1 <- CAPdiscrim (o.dist ~ Treatment, data = mta,
                                 dist="bray", axes=2, m=0, add=FALSE)

pdf(file="GHPBETAKNZ-16S",
    width=12,
    height=12)


plot1 <- ordiplot(Ordination.model1,type="none")
ordisymbol(plot1, mta, "Treatment", legend=TRUE,legend.x="topleft")
ordiellipse(Ordination.model1, groups = mta$Treatment, draw = "polygon",col=c("blue","orange"))
dev.off()

library(vegan)
library(devtools)


model1<-adonis(otuTABtranspose~Treatment,data=mta,permutations=999,method="bray")
model1



##nb.glms

library(mctoolsr)
library(edgeR)
library(vegan)
library(MASS)
library(lme4)
library(colorspace)
library(car)
library(doBy)
library(lmerTest)
library(lattice)
library(ggplot2)
library(ggthemes)
library(scales)
library(viridis)
library(RColorBrewer)
library(colorRamps)
library(reshape)
library(Hmisc)
library(emmeans)


input1<-load_taxa_table("~/Desktop/R/zotutab_16S.txt", "~/Desktop/R/16SMappingFile.txt")
controldrought = filter_data(input1, 'Treatment', filter_vals = 'Drought - Recovery 1')
controldrought2 = filter_data(controldrought, 'Treatment', filter_vals = 'Drought - Recovery 2')
controldrought3 = filter_data(controldrought2, 'Treatment', filter_vals = 'Control - Legacy Effects')
controldrought4 = filter_data(controldrought3, 'Treatment', filter_vals = 'Drought - Resistance to Second Drought')
controldrought5 = filter_data(controldrought4, 'Treatment', filter_vals = 'Drought - Second drought recovery 1')
final = filter_data(controldrought5, 'Treatment', filter_vals = 'Drought - Second drought recovery 2')
zotu_ITS_input = filter_data(final, 'Site', keep_vals = 'SGS')
zotu_ITS_input = filter_data(final, 'Site', keep_vals = 'CHY')
zotu_ITS_input = filter_data(final, 'Site', keep_vals = 'HYS')
zotu_ITS_input = filter_data(final, 'Site', keep_vals = 'KNZ')

##for relative abundance
#zotu_ITS_input$data_loaded = apply(zotu_ITS_input$data_loaded,2,function(x){100*x/sum(x)})




#ITS
y <- zotu_ITS_input$data_loaded
# y <- zotu_ITS_input$data_loaded[rowSums(zotu_ITS_input$data_loaded,na.rm = TRUE)>100,]
z <- zotu_ITS_input$taxonomy_loaded[rownames(zotu_ITS_input$taxonomy_loaded) %in% rownames(y),]
x <- DGEList(y)
x <- calcNormFactors(x)
x <- estimateCommonDisp(x)
zotu_ITS_TMM <- list(x$pseudo.counts,zotu_ITS_input$map_loaded, z)
names(zotu_ITS_TMM) <- names(zotu_ITS_input)


zotu_ITS_TMM$map_loaded$Treatment <- as.factor(zotu_ITS_TMM$map_loaded$Treatment)



##Taxonomy

tax_NBGLMs <- function(data.set) {
  results <- vector()
  for(i in 1:2){           #splits dataset for each taxonomic group
    print(paste(i, "of 2"))
    input_tax_level = summarize_taxonomy(data.set, level = i, report_higher_tax = TRUE, relative = TRUE)
    test.data.raw <-  cbind(data.set$map_loaded,t(input_tax_level))
    start.col <- dim(data.set$map_loaded)[2] + 1
    
    for(j in start.col:ncol(test.data.raw)){         #Does the glm and lsmean for each taxonomic group
      #print(j)
      
      test.glm <- tryCatch({
        suppressWarnings(glm.nb(test.data.raw[,j] ~ Treatment, data = test.data.raw))
      }, error = function(e){
        return(rep(NA, 8))
      })
      if(is(test.glm) != "negbin") {
        results <- rbind(results, test.glm);
        rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
        next
      }
      
      test.aov <- tryCatch({
        Anova(test.glm, test.statistic = "F")
      }, error = function(e){
        return(rep(NA, 8))
      })#, warning = function(w){
      #return(rep(NA, 20))
      #})
      if(is(test.aov)[1] != "anova") {
        results <- rbind(results, test.aov);
        rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
        next
      }
      
      test.lsm <- emmeans(test.glm, specs = "Treatment")
      pwc<-pairs(test.lsm)
      a <- c(test.aov$'Sum Sq', test.aov$Df, test.aov$F, test.aov$`Pr(>F)`,summary(test.lsm)$emmean, summary(test.lsm)$SE,summary(pwc)$p.value)
      results <- rbind(results, a)
      rownames(results)[nrow(results)] <- colnames(test.data.raw)[j]
      
    }
    
  }
  colnames(results) <- c(apply(expand.grid(rownames(test.aov), colnames(test.aov)), 1, paste, collapse="."), 
                         apply(expand.grid(summary(test.lsm)$Treatment, colnames(summary(test.lsm))[2:3]), 1, paste,collapse="."),
                         summary(pwc)$contrast)
  return(results)
}

#ITS
zotu_NBGLM_tax_ITS <- tax_NBGLMs(zotu_ITS_TMM)

###Individual OTUs
otu_NBGLMs <- function(data.set) {
  
  test.data.raw <-  cbind(data.set$map_loaded,t(data.set$data_loaded))
  start.col <- dim(data.set$map_loaded)[2] + 1
  
  results <- vector()
  
  for(j in start.col:ncol(test.data.raw)){         #Does the glm and lsmean for each taxonomic group
    
    if(j%%100 == 0) { print(paste(j, "out of", ncol(test.data.raw))) }
    
    if(sum(test.data.raw[,j], na.rm = TRUE) < 100) { next }
    
    #model <- formula(test.data.raw[,j] ~ RepNumber + Accession)
    
    test.glm <- tryCatch({
      suppressWarnings(glm.nb(test.data.raw[,j] ~ Treatment, data = test.data.raw))
    }, error = function(e){
      return(rep(NA, 8))
    })
    if(is(test.glm) != "negbin") {
      #results <- rbind(results, test.glm);
      #rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
      next
    }
    
    test.aov <- tryCatch({
      Anova(test.glm, test.statistic = "F")
    }, error = function(e){
      return(rep(NA, 8))
    })#, warning = function(w){
    #return(rep(NA, 20))
    #}
    if(is(test.aov)[1] != "anova") {
      #results <- rbind(results, test.aov);
      #rownames(results)[nrow(results)] <- colnames(test.data.raw)[j];
      next
    }
    
    test.lsm <- emmeans(test.glm, specs = "Treatment")
    pwc<-pairs(test.lsm)
    a <- c(test.aov$'Sum Sq', test.aov$Df, test.aov$F, test.aov$`Pr(>F)`,summary(test.lsm)$emmean, summary(test.lsm)$SE,summary(pwc)$p.value)
    results <- rbind(results, a)
    rownames(results)[nrow(results)] <- colnames(test.data.raw)[j]
  }
  
  colnames(results) <- c(apply(expand.grid(rownames(test.aov), colnames(test.aov)), 1, paste, collapse="."), 
                         apply(expand.grid(summary(test.lsm)$Treatment, colnames(summary(test.lsm))[2:3]), 1, paste,collapse="."),
                         summary(pwc)$contrast)
  return(results)
}

zotu_NBGLM_raw_ITS <- otu_NBGLMs(zotu_ITS_TMM)




##CHANGE NAMES
#Transform LSMean estimates back to counts instead of ln(counts)

#ITS
zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))] <- exp(zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))])
zotu_NBGLM_raw_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_raw_ITS))] <- exp(zotu_NBGLM_raw_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_raw_ITS))])

# FDR Adjustments


#ITS
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_tax_ITS))
w <- apply(zotu_NBGLM_tax_ITS[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_tax_ITS <- cbind(zotu_NBGLM_tax_ITS, w)
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_raw_ITS))
w <- apply(zotu_NBGLM_raw_ITS[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_raw_ITS <- cbind(zotu_NBGLM_raw_ITS, w)


rm(p.cols, w)


# Eta-Squared Calcs
####May want to combine this process with the above functions.
etaSq.calc <- function(nbglm.data){
  results <- vector()
  SS.cols <- grep("Sum Sq$", colnames(nbglm.data))
  for(i in 1:nrow(nbglm.data)){
    eta.sq <- nbglm.data[i,SS.cols] / sum(nbglm.data[i,SS.cols], na.rm = T)
    results <- rbind(results, eta.sq)
  }
  colnames(results) <- gsub("Sum Sq", "etaSq", colnames(results))
  rownames(results) <- rownames(nbglm.data)
  return(results)
}

#ITS
zotu_NBGLM_tax_ITS <- cbind(zotu_NBGLM_tax_ITS,  etaSq.calc(zotu_NBGLM_tax_ITS))
zotu_NBGLM_raw_ITS <- cbind(zotu_NBGLM_raw_ITS,  etaSq.calc(zotu_NBGLM_raw_ITS))


#Write Results to CSV
write.csv(zotu_NBGLM_tax_ITS, "zotu_NBGLM_otu_16S_KNZ_initial.csv")
write.csv(zotu_NBGLM_raw_ITS, "zotu_NBGLM_raw_16S_KNZ.csv")
