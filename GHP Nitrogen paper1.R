##Nitrogen Samples

inorganicNGHP<-read.csv(file.choose(),header=TRUE)
treatmentlabels<-read.csv(file.choose(),header=TRUE)
library(car)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(grid)
library(gridExtra)
library(patchwork)

##Start with just resistance

##merge samples by sample to get treatments added (this is faster than doing it by hand)
inorganicNwithlabels<-merge(inorganicNGHP,treatmentlabels,by="Sample")

##subset samples so you only have the treatments wanted: Control and Drought - Resistance 

inorganicNfinalpaper1<-subset(inorganicNwithlabels,New.Treatment!="Control - Legacy Effects"&New.Treatment!="Drought - Resistance"
                              &New.Treatment!="Drought - Resistance to Second Drought"&New.Treatment!="Drought - Recovery 1"&New.Treatment!="Drought - Recovery 2")


##Make a graph with all sites in one?

seFunc<-function(x){
  n<-sum(!is.na(x))
  ci<-sd(x,na.rm=T)/sqrt(n) 
  lims<-c(mean(x)-ci, mean(x)+ci)
  names(lims)<-c("ymin","ymax")
  return(lims)
}

##inorganic Nitrogen
a<-ggplot(inorganicNfinalpaper1,aes(Site,Inorganic,fill=New.Treatment)) +
 
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPinorganic.pdf",
       plot = a,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Ammonium

b<-ggplot(inorganicNfinalpaper1,aes(Site,Ammonium,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +

  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))


##Nitrate

c<-ggplot(inorganicNfinalpaper1,aes(Site,Nitrate,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    axis.title.y= element_blank(),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPinorganicNch.5-recovery.pdf",
       plot = b+c,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Stats 

ammonmodel<-lm(Ammonium~New.Treatment*Site,data=inorganicNfinalpaper1)
nitratemodel<-lm(Nitrate~New.Treatment*Site,data=inorganicNfinalpaper1)

Anova(ammonmodel,type=3)
emm<-emmeans(ammonmodel,pairwise~New.Treatment|Site)
eff_size(emm,sigma(ammonmodel),edf=3)


Anova(nitratemodel,type=3)
emm<-emmeans(nitratemodel,pairwise~New.Treatment|Site)
eff_size(emm,sigma(nitratemodel),edf=3)

##sensitivity analysis: need to take each site individually. And then subtract the means of each

SumStats1<-summarise(group_by(inorganicNfinalpaper1,New.Treatment,Site),
                    n=n(),
                    mean=mean(Ammonium),
                    sd=sd(Ammonium),
                    se=sd/sqrt(n))

SumStats2<-summarise(group_by(inorganicNfinalpaper1,New.Treatment,Site),
                     n=n(),
                     mean=mean(Nitrate),
                     sd=sd(Nitrate),
                     se=sd/sqrt(n))

library(effectsize)

sd_pooled()

##make graph by MAP. So make file with means from above by MAP. Then need to make curves... and do stats



##Bacteria and Fungi

qPCRGHP<-read.csv(file.choose(),header=TRUE)

##merge samples by sample to get treatments added (this is faster than doing it by hand)
qPCRwithlabels<-merge(qPCRGHP,treatmentlabels,by="Sample")

##subset samples so you only have the treatments wanted: Control and Drought - Resistance 

qPCRfinalpaper1<-subset(qPCRwithlabels,New.Treatment!="Control - Legacy Effects"&New.Treatment!="Drought - Second drought recovery 1"
                              &New.Treatment!="Drought - Second drought recovery 2"&New.Treatment!="Drought - Recovery 1"&New.Treatment!="Drought - Recovery 2")


qPCRfinalpaper1$Fungi<-as.numeric(qPCRfinalpaper1$Fungi)
qPCRfinalpaper1$Total<-as.numeric(qPCRfinalpaper1$Total)
qPCRfinalpaper1$BacteriaFungiRatio<-as.numeric(qPCRfinalpaper1$BacteriaFungiRatio)



##Bacteria

d<-ggplot(qPCRfinalpaper1,aes(Site,Bacteria,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Resistance","Drought - Resistance to Second Drought"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab("Copy number per gram of soil ") +
  
  theme(
    axis.text = element_text(color='black',size=10),
    legend.position = "none",
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPbacteria.pdf",
       plot = d,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Fungi

e<-ggplot(qPCRfinalpaper1,aes(Site,Fungi,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Resistance","Drought - Resistance to Second Drought"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.1,show.legend=F,position=position_dodge(0.5)) + 
  ylab("Copy number per gram of soil ") + 
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    axis.title.y = element_blank(),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
   
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPfungi.pdf",
       plot = e,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)
library(patchwork)
qpcrgraph<-d+e
ggsave(filename = "GHPfungiandbacteria.pdf",
       plot = qpcrgraph,
       bg = "transparent",
       width = 10, height = 6, units = "in",
       dpi = 600)


##Stats

bacteriamodel<-lm(Bacteria~New.Treatment*Site,data=qPCRfinalpaper1)

Anova(bacteriamodel,type=3)
emmeans(bacteriamodel,pairwise~New.Treatment|Site)

fungimodel<-lm(Fungi~New.Treatment*Site,data=qPCRfinalpaper1)

Anova(fungimodel,type=3)
emmeans(fungimodel,pairwise~New.Treatment|Site)

##Enzymes

EnzymesGHP<-read.csv(file.choose(),header=TRUE)

##merge samples by sample to get treatments added (this is faster than doing it by hand)
EnzymesGHPlabels<-merge(EnzymesGHP,treatmentlabels,by="Sample")

##subset samples so you only have the treatments wanted: Control and Drought - Resistance 

EnzymesGHPfinal<-subset(EnzymesGHPlabels,New.Treatment!="Control - Legacy Effects"&New.Treatment!="Drought - Resistance"
                        &New.Treatment!="Drought - Resistance to Second Drought"&New.Treatment!="Drought - Recovery 1"&New.Treatment!="Drought - Recovery 2")



carbonenzymes<-h+i+j+m

ggsave(filename = "carbonenyzmespanelghpch5-recovery.pdf",
       plot = carbonenzymes,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

##AG

h<-ggplot(EnzymesGHPfinal,aes(Site,AG,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))



##BG

i<-ggplot(EnzymesGHPfinal,aes(Site,BG,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    
    plot.title=element_text(size=20,face="bold",vjust=2))



##CB

j<-ggplot(EnzymesGHPfinal,aes(Site,CB,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))



##NAG

k<-ggplot(EnzymesGHPfinal,aes(Site,NAG,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))



##PHOS

l<-ggplot(EnzymesGHPfinal,aes(Site,PHOS,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))



##XYL

m<-ggplot(EnzymesGHPfinal,aes(Site,XYL,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))



##LAP

n<-ggplot(EnzymesGHPfinal,aes(Site,LAP,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Second drought recovery 1","Drought - Second drought recovery 2"),values=c("blue","orange","red")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))



nitrogenphos<-(k+n)/l
ggsave(filename = "GHPnitrogencombinedch5-recovery.pdf",
       plot = nitrogenphos,
       bg = "transparent",
       width = 10, height = 6, units = "in",
       dpi = 600)

##Carbon

o<-ggplot(EnzymesGHPfinal,aes(Site,Carbon,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Control - Legacy Effects"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPCarboncombinedrecovery2.pdf",
       plot = o,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Nitrogen

p<-ggplot(EnzymesGHPfinal,aes(Site,Nitrogen,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Control - Legacy Effects"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPNitrogenaddedrecovery2.pdf",
       plot = p,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Stats

AGmodel<-lm(AG~New.Treatment*Site,data=EnzymesGHPfinal)
BGmodel<-lm(BG~New.Treatment*Site,data=EnzymesGHPfinal)
CBmodel<-lm(CB~New.Treatment*Site,data=EnzymesGHPfinal)
NAGmodel<-lm(NAG~New.Treatment*Site,data=EnzymesGHPfinal)
LAPmodel<-lm(LAP~New.Treatment*Site,data=EnzymesGHPfinal)
PHOSmodel<-lm(PHOS~New.Treatment*Site,data=EnzymesGHPfinal)
XYLmodel<-lm(XYL~New.Treatment*Site,data=EnzymesGHPfinal)
Carbonmodel<-lm(Carbon~New.Treatment*Site,data=EnzymesGHPfinal)
Nitrogenmodel<-lm(Nitrogen~New.Treatment*Site,data=EnzymesGHPfinal)

Anova(AGmodel,type=3)
emm<-emmeans(AGmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(AGmodel),edf=3)

Anova(BGmodel,type=3)
emm<-emmeans(BGmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(BGmodel),edf=3)

Anova(CBmodel,type=3)
emm<-emmeans(CBmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(CBmodel),edf=3)

Anova(XYLmodel,type=3)
emm<-emmeans(XYLmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(XYLmodel),edf=3)

Anova(NAGmodel,type=3)
emm<-emmeans(NAGmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(NAGmodel),edf=3)

Anova(LAPmodel,type=3)
emm<-emmeans(LAPmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(LAPmodel),edf=3)

Anova(PHOSmodel,type=3)
emm<-emmeans(PHOSmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(PHOSmodel),edf=3)

Anova(Carbonmodel,type=3)
emm<-emmeans(Carbonmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(Carbonmodel),edf=3)

Anova(Nitrogenmodel,type=3)
emm<-emmeans(Nitrogenmodel,pairwise~New.Treatment|Site)
emm
eff_size(emm,sigma(Nitrogenmodel),edf=3)



##Recovery Week 1

##Nitrogen Samples

inorganicNGHP<-read.csv(file.choose(),header=TRUE)
treatmentlabels<-read.csv(file.choose(),header=TRUE)
library(car)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(grid)
library(gridExtra)
library(patchwork)

##Start with just resistance

##merge samples by sample to get treatments added (this is faster than doing it by hand)
inorganicNwithlabels<-merge(inorganicNGHP,treatmentlabels,by="Sample")

##subset samples so you only have the treatments wanted: Control and Drought - Resistance 

inorganicNfinalpaper1<-subset(inorganicNwithlabels,New.Treatment!="Control - Legacy Effects"&New.Treatment!="Drought - Resistance to Second Drought"&New.Treatment!="Drought - Second drought recovery 1"
                              &New.Treatment!="Drought - Second drought recovery 2"&New.Treatment!="Drought - Recovery 2"&New.Treatment!="Drought - Resistance")

##Make a graph with all sites in one?

seFunc<-function(x){
  n<-sum(!is.na(x))
  ci<-sd(x,na.rm=T)/sqrt(n) 
  lims<-c(mean(x)-ci, mean(x)+ci)
  names(lims)<-c("ymin","ymax")
  return(lims)
}

##inorganic Nitrogen
a<-ggplot(inorganicNfinalpaper1,aes(Site,Inorganic,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  scale_fill_manual("",breaks=c("Control","Drought - Resistance"),values=c("blue","orange")) +
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPinorganic.pdf",
       plot = a,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Ammonium

b<-ggplot(inorganicNfinalpaper1,aes(Site,Ammonium,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 2"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))


##Nitrate

c<-ggplot(inorganicNfinalpaper1,aes(Site,Nitrate,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 2"),values=c("blue","orange")) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    axis.title.y= element_blank(),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPinorganicNrecovery2.pdf",
       plot = b+c,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Stats 

ammonmodel<-lm(Ammonium~New.Treatment*Site,data=inorganicNfinalpaper1)
nitratemodel<-lm(Nitrate~New.Treatment*Site,data=inorganicNfinalpaper1)

Anova(ammonmodel,type=3)
emm1<-emmeans(ammonmodel,pairwise~New.Treatment|Site)
eff_size(emm1,sigma(ammonmodel),edf=3)

Anova(nitratemodel,type=3)
emm2<-emmeans(nitratemodel,pairwise~New.Treatment|Site)
eff_size(emm2,sigma(nitratemodel),edf=3)

##Recovery Week 2

##Nitrogen Samples

inorganicNGHP<-read.csv(file.choose(),header=TRUE)
treatmentlabels<-read.csv(file.choose(),header=TRUE)
library(car)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(emmeans)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(Hmisc)
library(grid)
library(gridExtra)
library(patchwork)

##Start with just resistance

##merge samples by sample to get treatments added (this is faster than doing it by hand)
inorganicNwithlabels<-merge(inorganicNGHP,treatmentlabels,by="Sample")

##subset samples so you only have the treatments wanted: Control and Drought - Resistance 

inorganicNfinalpaper1<-subset(inorganicNwithlabels,New.Treatment!="Control - Legacy Effects"&New.Treatment!="Drought - Resistance to Second Drought"&New.Treatment!="Drought - Second drought recovery 1"
                              &New.Treatment!="Drought - Second drought recovery 2"&New.Treatment!="Drought - Recovery 1"&New.Treatment!="Drought - Resistance")

##Make a graph with all sites in one?

seFunc<-function(x){
  n<-sum(!is.na(x))
  ci<-sd(x,na.rm=T)/sqrt(n) 
  lims<-c(mean(x)-ci, mean(x)+ci)
  names(lims)<-c("ymin","ymax")
  return(lims)
}

##inorganic Nitrogen
a<-ggplot(inorganicNfinalpaper1,aes(Site,Inorganic,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  scale_fill_manual("",breaks=c("Control","Drought - Resistance"),values=c("blue","orange")) +
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPinorganic.pdf",
       plot = a,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Ammonium

b<-ggplot(inorganicNfinalpaper1,aes(Site,Ammonium,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))


##Nitrate

c<-ggplot(inorganicNfinalpaper1,aes(Site,Nitrate,fill=New.Treatment)) +
  
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote(''*mu~ 'g N/ g dry soil'  )) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    axis.title.y= element_blank(),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPinorganicNrecovery1.pdf",
       plot = b+c,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Stats 

ammonmodel<-lm(Ammonium~New.Treatment*Site,data=inorganicNfinalpaper1)
nitratemodel<-lm(Nitrate~New.Treatment*Site,data=inorganicNfinalpaper1)

Anova(ammonmodel,type=3)
emm3<-emmeans(ammonmodel,pairwise~New.Treatment|Site)
eff_size(emm3,sigma(ammonmodel),edf=3)

Anova(nitratemodel,type=3)
emm4<-emmeans(nitratemodel,pairwise~New.Treatment|Site)
eff_size(emm4,sigma(nitratemodel),edf=3)

EnzymesGHP<-read.csv(file.choose(),header=TRUE)

##merge samples by sample to get treatments added (this is faster than doing it by hand)
EnzymesGHPlabels<-merge(EnzymesGHP,treatmentlabels,by="Sample")

##subset samples so you only have the treatments wanted: Control and Drought - Resistance 

EnzymesGHPfinal<-subset(EnzymesGHPlabels,New.Treatment!="Control - Legacy Effects"&New.Treatment!="Drought - Resistance to Second Drought"&New.Treatment!="Drought - Second drought recovery 1"
                        &New.Treatment!="Drought - Second drought recovery 2"&New.Treatment!="Drought - Resistance"&New.Treatment!="Drought - Recovery 2")

carbonenzymes<-h+i+j+m

ggsave(filename = "carbonenyzmespanelghp.pdf",
       plot = carbonenzymes,
       bg = "transparent",
       width = 12, height = 6, units = "in",
       dpi = 600)

##AG

h<-ggplot(EnzymesGHPfinal,aes(Site,AG,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPAG.pdf",
       plot = h,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##BG

i<-ggplot(EnzymesGHPfinal,aes(Site,BG,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPBG.pdf",
       plot = i,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##CB

j<-ggplot(EnzymesGHPfinal,aes(Site,CB,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPCB.pdf",
       plot = j,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##NAG

k<-ggplot(EnzymesGHPfinal,aes(Site,NAG,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPNAG.pdf",
       plot = k,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##PHOS

l<-ggplot(EnzymesGHPfinal,aes(Site,PHOS,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPPHOS.pdf",
       plot = l,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##XYL

m<-ggplot(EnzymesGHPfinal,aes(Site,XYL,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPXYL.pdf",
       plot = m,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##LAP

n<-ggplot(EnzymesGHPfinal,aes(Site,LAP,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=10),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPLAP.pdf",
       plot = n,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

nitrogenphos<-(k+n)/l
ggsave(filename = "GHPnitrogencombined.pdf",
       plot = nitrogenphos,
       bg = "transparent",
       width = 10, height = 6, units = "in",
       dpi = 600)

##Carbon

o<-ggplot(EnzymesGHPfinal,aes(Site,Carbon,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPCarboncombined.pdf",
       plot = o,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Nitrogen

p<-ggplot(EnzymesGHPfinal,aes(Site,Nitrogen,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 1"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab(bquote('Total nmol activity/g dry soil/ h')) +
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPNitrogenadded.pdf",
       plot = p,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Stats

AGmodel<-lm(AG~New.Treatment*Site,data=EnzymesGHPfinal)
BGmodel<-lm(BG~New.Treatment*Site,data=EnzymesGHPfinal)
CBmodel<-lm(CB~New.Treatment*Site,data=EnzymesGHPfinal)
NAGmodel<-lm(NAG~New.Treatment*Site,data=EnzymesGHPfinal)
LAPmodel<-lm(LAP~New.Treatment*Site,data=EnzymesGHPfinal)
PHOSmodel<-lm(PHOS~New.Treatment*Site,data=EnzymesGHPfinal)
XYLmodel<-lm(XYL~New.Treatment*Site,data=EnzymesGHPfinal)
Carbonmodel<-lm(Carbon~New.Treatment*Site,data=EnzymesGHPfinal)
Nitrogenmodel<-lm(Nitrogen~New.Treatment*Site,data=EnzymesGHPfinal)

Anova(AGmodel,type=3)
emm1<-emmeans(AGmodel,pairwise~New.Treatment|Site)
eff_size(emm1,sigma(AGmodel),edf=3)

Anova(BGmodel,type=3)
emmeans(BGmodel,pairwise~New.Treatment|Site)

Anova(CBmodel,type=3)
emmeans(CBmodel,pairwise~New.Treatment|Site)

Anova(XYLmodel,type=3)
emmeans(XYLmodel,pairwise~New.Treatment|Site)

Anova(NAGmodel,type=3)
emmeans(NAGmodel,pairwise~New.Treatment|Site)

Anova(LAPmodel,type=3)
emmeans(LAPmodel,pairwise~New.Treatment|Site)

Anova(PHOSmodel,type=3)
emmeans(PHOSmodel,pairwise~New.Treatment|Site)

Anova(Carbonmodel,type=3)
emmeans(Carbonmodel,pairwise~New.Treatment|Site)

Anova(Nitrogenmodel,type=3)
emmeans(Nitrogenmodel,pairwise~New.Treatment|Site)

qPCRGHP<-read.csv(file.choose(),header=TRUE)

##merge samples by sample to get treatments added (this is faster than doing it by hand)
qPCRwithlabels<-merge(qPCRGHP,treatmentlabels,by="Sample")

##subset samples so you only have the treatments wanted: Control and Drought - Resistance 

qPCRfinalpaper1<-subset(qPCRwithlabels,New.Treatment!="Control - Legacy Effects"&New.Treatment!="Drought - Resistance to Second Drought"&New.Treatment!="Drought - Second drought recovery 1"
                        &New.Treatment!="Drought - Second drought recovery 2"&New.Treatment!="Drought - Resistance"&New.Treatment!="Drought - Recovery 1")


qPCRfinalpaper1$Fungi<-as.numeric(qPCRfinalpaper1$Fungi)
qPCRfinalpaper1$Total<-as.numeric(qPCRfinalpaper1$Total)
qPCRfinalpaper1$BacteriaFungiRatio<-as.numeric(qPCRfinalpaper1$BacteriaFungiRatio)



##Bacteria

d<-ggplot(qPCRfinalpaper1,aes(Site,Bacteria,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 2"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.05,show.legend=F,position=position_dodge(0.5)) + 
  ylab("Copy number per gram of soil ") +
  
  theme(
    axis.text = element_text(color='black',size=10),
    legend.position = "none",
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPbacteria.pdf",
       plot = d,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

##Fungi

e<-ggplot(qPCRfinalpaper1,aes(Site,Fungi,fill=New.Treatment)) +
  scale_fill_manual("",breaks=c("Control","Drought - Recovery 2"),values=c("blue","orange")) +
  scale_x_discrete(limits = c("SGS", "CHY", "HYS","KNZ"))+
  stat_summary(geom="bar",fun.y="mean",width=0.5,color="black",position=position_dodge()) +
  stat_summary(geom="errorbar",fun.data=seFunc,width=0.1,show.legend=F,position=position_dodge(0.5)) + 
  ylab("Copy number per gram of soil ") + scale_y_log10()+
  
  theme(
    axis.text = element_text(color='black',size=10),
    axis.title = element_text(color='black',size=12),
    axis.ticks = element_line(color='black'),
    axis.title.y = element_blank(),
    legend.title = element_text(size=18),
    panel.background = element_rect(fill=NA,color='black'),
    panel.grid = element_blank(),
    
    plot.title=element_text(size=20,face="bold",vjust=2))

ggsave(filename = "GHPfungi.pdf",
       plot = e,
       bg = "transparent",
       width = 8, height = 6, units = "in",
       dpi = 600)

qpcrgraph<-d+e
ggsave(filename = "GHPfungiandbacteriarecovery2.pdf",
       plot = qpcrgraph,
       bg = "transparent",
       width = 10, height = 6, units = "in",
       dpi = 600)


##Stats

bacteriamodel<-lm(Bacteria~New.Treatment*Site,data=qPCRfinalpaper1)

Anova(bacteriamodel,type=3)
emmeans(bacteriamodel,pairwise~New.Treatment|Site)

fungimodel<-lm(Fungi~New.Treatment*Site,data=qPCRfinalpaper1)

Anova(fungimodel,type=3)
emmeans(fungimodel,pairwise~New.Treatment|Site)




