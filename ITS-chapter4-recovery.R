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
input <- load_taxa_table("~/Desktop/R/zotutab_ITS.txt", "~/Desktop/R/ITSMappingFile.txt")

colSums(input$data_loaded)
inputr = single_rarefy(input, 4000)

##Filter down
Control = filter_data(inputr, 'Treatment', keep_vals = 'Control')
ControlSGS = filter_data(Control, 'Site', keep_vals = 'SGS')
ControlCHY = filter_data(Control, 'Site', keep_vals = 'CHY')
ControlHYS = filter_data(Control, 'Site', keep_vals = 'HYS')
ControlKNZ = filter_data(Control, 'Site', keep_vals = 'KNZ')
Drought = filter_data(inputr, 'Treatment', keep_vals = 'Drought - Recovery 1')
DroughtSGS = filter_data(Drought, 'Site', keep_vals = 'SGS')
DroughtCHY = filter_data(Drought, 'Site', keep_vals = 'CHY')
DroughtHYS = filter_data(Drought, 'Site', keep_vals = 'HYS')
DroughtKNZ = filter_data(Drought, 'Site', keep_vals = 'KNZ')

ControlSGSshan <- diversity(t(ControlSGS$data_loaded), index = "shannon")
ControlCHYshan <- diversity(t(ControlCHY$data_loaded), index = "shannon")
ControlHYSshan <- diversity(t(ControlHYS$data_loaded), index = "shannon")
ControlKNZshan <- diversity(t(ControlKNZ$data_loaded), index = "shannon")

DroughtSGSshan <- diversity(t(DroughtSGS$data_loaded), index = "shannon")
DroughtCHYshan <- diversity(t(DroughtCHY$data_loaded), index = "shannon")
DroughtHYSshan <- diversity(t(DroughtHYS$data_loaded), index = "shannon")
DroughtKNZshan <- diversity(t(DroughtKNZ$data_loaded), index = "shannon")

colors<-c("blue","orange","blue","orange","blue","orange","blue","orange")

pdf(file="Shannon Diversity chapter 4 recovery1",
    width=10,
    height=8)
a<-boxplot(ControlSGSshan, DroughtSGSshan, ControlCHYshan, DroughtCHYshan, ControlHYSshan, DroughtHYSshan, ControlKNZshan, DroughtKNZshan, col=colors, main = "Treatment", ylab = "Shannon's Diversity",ylim=c(2.5,5), xlab = "Treament", 
           names = c("SGS-C", "SGS-R1", "CHY-C","CHY-R1","HYS-C","HYS-R1","KNZ-C","KNZ-R1"))
dev.off()

##STATS
controldrought = filter_data(inputr, 'Treatment', filter_vals = 'Drought - Resistance')
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

ControlSGSrich <- specnumber(t(ControlSGS$data_loaded))
ControlCHYrich <- specnumber(t(ControlCHY$data_loaded))
ControlHYSrich <- specnumber(t(ControlHYS$data_loaded))
ControlKNZrich <- specnumber(t(ControlKNZ$data_loaded))

DroughtSGSrich <- specnumber(t(DroughtSGS$data_loaded))
DroughtCHYrich <- specnumber(t(DroughtCHY$data_loaded))
DroughtHYSrich <- specnumber(t(DroughtHYS$data_loaded))
DroughtKNZrich <- specnumber(t(DroughtKNZ$data_loaded))

colors<-c("blue","orange","blue","orange","blue","orange","blue","orange")

pdf(file="Richness chapter 4 all",
    width=10,
    height=8)
boxplot(ControlSGSrich, DroughtSGSrich, ControlCHYrich, DroughtCHYrich, ControlHYSrich, DroughtHYSrich, ControlKNZrich, DroughtKNZrich, col=colors, main = "Treatment", ylab = "Richness", xlab = "Treament", 
        names = c("SGS-C", "SGS-R1", "CHY-C","CHY-R1","HYS-C","HYS-R1","KNZ-C","KNZ-R1"))
dev.off()

allshan <- specnumber(t(controldrought6$data_loaded))

allshan <- allshan[match(names(allshan), row.names(controldrought6$map_loaded))]
allshan_df <- cbind(controldrought6$map_loaded, allshan)

controldrought6$allshan <- allshan

allshan_glm <- glm(allshan ~ Treatment*Site, data = allshan_df)
Anova(allshan_glm, test.statistic = "F")

emmeans<-emmeans(allshan_glm,pairwise~Treatment|Site)
eff_size(emmeans,sigma(allshan_glm),edf=3)

library(dplyr)
library(BiodiversityR)
library (vegan)
library (MASS)

## All data
input <- load_taxa_table("~/Desktop/R/zotutab_ITS.txt", "~/Desktop/R/ITSMappingFile.txt")

colSums(input$data_loaded)
inputr = single_rarefy(input, 4000)

controldrought = filter_data(inputr, 'Treatment', filter_vals = 'Drought - Resistance')
controldrought2 = filter_data(controldrought, 'Treatment', filter_vals = 'Drought - Recovery 2')
controldrought3 = filter_data(controldrought2, 'Treatment', filter_vals = 'Control - Legacy Effects')
controldrought4 = filter_data(controldrought3, 'Treatment', filter_vals = 'Drought - Resistance to Second Drought')
controldrought5 = filter_data(controldrought4, 'Treatment', filter_vals = 'Drought - Second drought recovery 1')
final = filter_data(controldrought5, 'Treatment', filter_vals = 'Drought - Second drought recovery 2')

SGS = filter_data(final, 'Site', keep_vals = 'SGS')
CHY = filter_data(final, 'Site', keep_vals = 'CHY')
HYS = filter_data(final, 'Site', keep_vals = 'HYS')
KNZ = filter_data(final, 'Site', keep_vals = 'KNZ')

otuTAB<-HYS$data_loaded

otuTAB <- as.data.frame (otuTAB)
otuTABtranspose<-t(otuTAB)
mta<-HYS$map_loaded
mta$Treatment <- as.factor (mta$Treatment)

o.dist <-vegdist (otuTABtranspose)
Ordination.model1 <- CAPdiscrim (o.dist ~ Treatment, data = mta,
                                 dist="bray", axes=2, m=0, add=FALSE)

pdf(file="GHPBETAKNZRecovery",
    width=8,
    height=8)


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


input1<-load_taxa_table("~/Desktop/R/zotutab_ITS.txt", "~/Desktop/R/ITSMappingFile.txt")
controldrought = filter_data(input1, 'Treatment', filter_vals = 'Drought - Resistance')
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
zotu_ITS_input$data_loaded = apply(zotu_ITS_input$data_loaded,2,function(x){100*x/sum(x)})




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
    input_tax_level = summarize_taxonomy(data.set, level = i, report_higher_tax = TRUE, relative = FALSE)
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


##CHANGE NAMES
#Transform LSMean estimates back to counts instead of ln(counts)

#ITS
zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))] <- exp(zotu_NBGLM_tax_ITS[,grep(".emmean|.SE", colnames(zotu_NBGLM_tax_ITS))])

# FDR Adjustments


#ITS
p.cols <- grep(".Pr[:(:][:>:]F[:):]", colnames(zotu_NBGLM_tax_ITS))
w <- apply(zotu_NBGLM_tax_ITS[,p.cols], 2, function(x){ p.adjust(x, method = 'fdr') })
colnames(w) <- paste(colnames(w), "FDR", sep = ".")
zotu_NBGLM_tax_ITS <- cbind(zotu_NBGLM_tax_ITS, w)

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



#Write Results to CSV
write.csv(zotu_NBGLM_tax_ITS, "zotu_NBGLM_phylum_ITS_KNZ_recovery1.csv")




