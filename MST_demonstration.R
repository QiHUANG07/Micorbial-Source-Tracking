##Microbial source tracking

#R codes for data analyses and plotting
#NOTE: the R script below has been simplified for illustration

###################
#library preparation#
###################
#install.packages("ggplot2")
#install.packages("reshape")
#install.packages("yardstick")
#install.packages("dplyr")
library(ggplot2)
library(reshape2)
library(yardstick)
library(dplyr)

setwd("D:/Huang Qi/Fight/Microbial Source Tracking/EMP data/paper draft/Rescue/Additional file 3/")

source('sourcetracker/src/SourceTracker.r') #load SourceTracker

################################
#Leave-one-out cross-validation#----
################################
#data preparation
inte_otu_abu = read.table(file="3654_source_dataset.csv", header=T, row.names=1, sep=",",check.names=FALSE, quote="")
inte_otu_abu<-as.matrix(inte_otu_abu) 
class(inte_otu_abu)<-"numeric"
save(inte_otu_abu, file="inte_otu_abu.RData")
load(file="inte_otu_abu.RData")

otu_3654 <- data.frame(t(inte_otu_abu),check.names=FALSE)#save the otu information for later prediction for real dataset
otu_3654$OTU_ID = rownames(otu_3654)
save(otu_3654, file="otu_3654.RData")
load(file="otu_3654.RData")

design = read.table("design.csv",  header=T, sep=",",check.names=FALSE)

#cross-validation
pred_list<-list()
proptab_list<-list()

for (i in 1:nrow(inte_otu_abu))
{
  st<-sourcetracker(inte_otu_abu[-i,], design$empo_3[-i])
  st_predic1<- predict(st,inte_otu_abu[i,], alpha1=0.001, alpha2=0.001)
  st_predic2 <- as.data.frame(st_predic1[[2]])
  st_predic2$sampleID<-rownames(inte_otu_abu)[i]
  st_predic2<-merge(st_predic2,design,by="sampleID")
  pred_list[[i]]<-st_predic1
  proptab_list[[i]]<-st_predic2
}

#two functions for data processing of sourcetracking results
maxname<-function(data_frame){
  order<-c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human Fecal","Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)","Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)","Water (non-saline)","Water (saline)","WWTPs","Unknown")
  max_name<-order[which.max(data_frame[2:17])]
  return(max_name)}

maxratio<-function(data_frame){
  max_ratio<-data_frame[2:17][which.max(data_frame[2:17])]
  return(max_ratio)}

proptab<-do.call("rbind",proptab_list)
proptab$predicted<-apply(proptab,1,maxname)
proptab$matched<-proptab$predicted == proptab$empo_3
proptab$predictedratio<-apply(proptab,1,maxratio)

save(proptab,file="Traindataset_proptab_empo3.RData")
load(file="Traindataset_proptab_empo.RData")

#plot the barchart
proptab_T<-subset(proptab,matched%in%"TRUE")
class(proptab_T$predictedratio)<-"double"
plotdata_stprob <- aggregate(proptab_T$predictedratio,by = list(eco_type=proptab_T$empo_3), FUN = function(x) c(mean = mean(x), sd = sd(x),n = length(x)))
plotdata_stprob <- do.call(data.frame, plotdata_stprob)
colnames(plotdata_stprob) <- c("empo_3", "mean", "sd", "n")

limits <- aes(ymax = mean + sd,ymin = mean - sd)
ggplot(data = plotdata_stprob, aes(x = empo_3, y = mean, fill = empo_3))+
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(limits, width = 0.3,color="gray50") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill='gray96', colour='gray'),
        text=element_text(size=25,color = "black"))+
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_fill_manual(values=c("#696969","#ffff02","#8b4512","#c98c23","#8ac8ed",
                             "#7cfc02","#006400","#ffa07a","#ee0000","#8FBC8F",
                             "#e46bef","#8a1894","#4169e0","#0a0a9f","#5c6b06"), name="Environmental types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","WWTPs"))+
  scale_x_discrete(labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                            "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                            "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                            "Water (non-saline)","Water (saline)","WWTPs"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust = 1))+
  ylab("Probability")+
  xlab("")

#correct ratio calculation
correct_ratio<-matrix(nrow=3,ncol=15)
correct_ratio<-data.frame(correct_ratio)
colnames(correct_ratio)<-(colnames(proptab))[2:16]
rownames(correct_ratio)<-c("samplesize","correct_sample_size","correct_ratio")

for (i in 2:16)
{
  correct_ratio[1,i-1]=nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
  correct_ratio[2,i-1]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$matched=="TRUE"),])
  correct_ratio[3,i-1]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$matched=="TRUE"),])/nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
}

correct_ratio_t=data.frame(t(correct_ratio),check.names = F)
correct_ratio_t$env_type=rownames(correct_ratio_t)

ggplot(correct_ratio_t, aes(x="", y = correct_ratio, fill=env_type)) +
  geom_bar(width = 1, stat = "identity")+
  ylim(0,1)+
  scale_fill_manual(values=c("#696969","#ffff02","#8b4512","#c98c23","#8ac8ed",
                             "#7cfc02","#006400","#ffa07a","#ee0000","#8FBC8F",
                             "#e46bef","#8a1894","#4169e0","#0a0a9f","#5c6b06"), name="Environmental types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","WWTPs"))+
  scale_x_discrete(labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                            "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                            "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                            "Water (non-saline)","Water (saline)","WWTPs"))+
  coord_polar(theta = "y",direction = -1)+
  facet_grid(.~ env_type)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())
graphics.off()

#Sensitivity and specificity
SS<-matrix(nrow=16,ncol=3)
SS<-data.frame(SS)
rownames(SS)<-(colnames(proptab))[2:17]
colnames(SS)<-c("Sensitivity","Specificity","Sample size")
for (i in 2:17)
{
  SS[i-1,3]=nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
  SS[i-1,1]=nrow(proptab[(proptab$predicted==(colnames(proptab))[i])&(proptab$empo_3==(colnames(proptab))[i]),])/nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
  SS[i-1,2]=nrow(proptab[(proptab$predicted!=(colnames(proptab))[i])&(proptab$empo_3!=(colnames(proptab))[i]),])/nrow(proptab[proptab$empo_3!=(colnames(proptab))[i],])
}

APSS<-matrix(nrow=16,ncol=10)
APSS<-data.frame(APSS)
rownames(APSS)<-(colnames(proptab))[2:17]
colnames(APSS)<-c("TP","FP","TN","FN","Precision","Recall","F1","Sensitivity","Specificity","Sample size")
for (i in 2:17)
{
  APSS[i-1,1]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$predicted==(colnames(proptab))[i]),])#TP
  APSS[i-1,2]=nrow(proptab[(proptab$empo_3!=(colnames(proptab))[i])&(proptab$predicted==(colnames(proptab))[i]),])#FP
  APSS[i-1,3]=nrow(proptab[(proptab$empo_3!=(colnames(proptab))[i])&(proptab$predicted!=(colnames(proptab))[i]),])#TN
  APSS[i-1,4]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$predicted!=(colnames(proptab))[i]),])#FN
  APSS[i-1,5]=APSS[i-1,1]/(APSS[i-1,1]+ APSS[i-1,2]) #P = TP/(TP+FP)
  APSS[i-1,6]=APSS[i-1,1]/(APSS[i-1,1]+ APSS[i-1,4]) #R = TP/(TP+FN)
  APSS[i-1,7]=2*APSS[i-1,5]* APSS[i-1,6]/( APSS[i-1,5]+ APSS[i-1,6]) # F1 = 2P*R/(P+R)
  APSS[i-1,8]=APSS[i-1,1]/(APSS[i-1,1]+APSS[i-1,4]) #sensitivity
  APSS[i-1,9]=APSS[i-1,3]/(APSS[i-1,2]+APSS[i-1,3]) #specificity
  APSS[i-1,10]=nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
}

PRF<-matrix()
PRF<-data.frame(PRF)
PRF$macro_precision =mean(APSS$Precision[1:15]) 
PRF$marcro_reall = mean(APSS$Recall[1:15])
PRF$macro_F1 = mean(APSS$F1[1:15])
PRF$micro_precision=sum(APSS$TP[1:15])/sum(sum(APSS$TP[1:15]),sum(APSS$FP[1:15]))
PRF$micro_recall=sum(APSS$TP[1:15])/sum(sum(APSS$TP[1:15]),sum(APSS$FN[1:15]))
PRF$micro_F1=2*micro_precision*micro_recall/(micro_precision+micro_recall)

#ROC curve
library(yardstick)
library(dplyr)
autoplot(roc_curve(proptab, empo_3, 2:16))
roc_auc(proptab, empo_3, 2:16) 
roc_auc(proptab, empo_3, 2:16, estimator = "macro_weighted")

proptab_roc=roc_curve(proptab, empo_3, 2:16)
proptab_roc$.level=factor(proptab_roc$.level,
                          levels = c("Animal corpus","Animal fecal", "Human Fecal","WWTPs",
                                     "Plant corpus","Plant rhizosphere",
                                     "Aerosol (non-saline)","Sediment (non-saline)","Soil (non-saline)","Surface (non-saline)","Water (non-saline)",
                                     "Hypersaline (saline)","Sediment (saline)","Surface (saline)","Water (saline)"))

proptab_roc$.level=factor(proptab_roc$.level,
                          levels = c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human Fecal",
                                     "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                                     "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                                     "Water (non-saline)","Water (saline)","WWTPs"))


ggplot(proptab_roc,aes(x = 1 - specificity, y = sensitivity,color=.level)) +
  geom_path(size=1.2) +
  geom_abline(lty = 3) +
  coord_equal() +
  theme_bw() +
  facet_wrap( ~ .level, ncol=4)+
  scale_color_manual(values=c("#696969","#ffff02","#8b4512","#c98c23","#8ac8ed",
                              "#7cfc02","#006400","#ffa07a","#ee0000","#8FBC8F",
                              "#e46bef","#8a1894","#4169e0","#0a0a9f","#5c6b06"), name="Environmental types",
                     labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                              "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                              "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                              "Water (non-saline)","Water (saline)","WWTPs"))+
  theme(axis.text.x = element_text(size=11, ),
        axis.text.y = element_text(size=11),
        axis.title = element_text(size = 14),
        legend.text=element_text(size=13, family="Arial"),
        legend.title = element_text(size=13),
        strip.background = element_rect(fill = 'white'),
        strip.text=element_text(size = 11),
        panel.grid=element_blank())+
  xlab("1-Specificity")+
  ylab("Sensitivity")


#pack the trained classifier
st_empo3<-sourcetracker(inte_otu_abu, design$empo_3)
save(st_empo3, file="st_empo3.RData")
load(file="st_empo3.RData")

#################################
##Parameter tuning#----
#################################
#Orthogonal experiment of the parameters using a 75-sample subdataset#
#extract 75 datasets
load(file="st_empo3.RData")
load(file="inte_otu_abu.RData")

setwd("./parameter_tuning_orthogonal_design")
subdb <- c()
for (i in 1:15)
{
  empo=c("Animal corpus","Animal fecal", "Human Fecal","WWTPs",
         "Plant corpus","Plant rhizosphere",
         "Aerosol (non-saline)","Sediment (non-saline)","Soil (non-saline)","Surface (non-saline)","Water (non-saline)",
         "Hypersaline (saline)","Sediment (saline)","Surface (saline)","Water (saline)")
  subdb=c(subdb,sample(rownames(subset(design,empo_3==empo[i])),size=5,replace = FALSE))
}


subdb <- data.frame(subdb,stringsAsFactors = F)
save(subdb,file = "subdb.RData")

load("subdb.RData")

pred_list<-list()
proptab_list<-list()

for (j in 1:nrow(subdb))
{
  index=which(design==subdb[j,1],arr.ind=TRUE)
  i=index[1]
  st<-sourcetracker(inte_otu_abu[-i,], design$empo_3[-i])
  st_predic1<- predict(st,inte_otu_abu[i,], alpha1=0.001, alpha2=0.001,beta=10,burnin=50,rarefaction_depth=3000,nrestarts=10) #set different parameters here
  st_predic2 <- as.data.frame(st_predic1[[2]])
  st_predic2$sampleID<-rownames(inte_otu_abu)[i]
  st_predic2<-merge(st_predic2,design,by="sampleID")
  pred_list[[i]]<-st_predic1
  proptab_list[[i]]<-st_predic2
}

proptab_list_1 <- proptab_list
pred_list_1 <- pred_list

save(proptab_list_1, file="proptab_list_1.RData")
save(pred_list_1, file="pred_list_1.RData")

load("proptab_list_1.RData")
#two functions for data processing of sourcetracking results
maxname<-function(data_frame){
  order<-c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human Fecal","Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)","Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)","Water (non-saline)","Water (saline)","WWTPs","Unknown")
  max_name<-order[which.max(data_frame[1:16])]
  return(max_name)}

maxratio<-function(data_frame){
  max_ratio<-data_frame[2:17][which.max(data_frame[2:17])]
  return(max_ratio)}

proptab<-do.call("rbind",proptab_list)
proptab$predicted<-apply(proptab,1,maxname)
proptab<-merge(proptab,design_test,by="sampleID")
proptab$matched<-proptab$predicted == proptab$empo_3
proptab$predictedratio<-apply(proptab,1,maxratio)

#Parameter adjustment result using three simulated samples#
setwd("parameter_tuning_with_three_simulated_samples/dsA")
load("dsA.RData")
#load("dsB.RData")
#load("dsC.RData")

pred_list_5<-list()
proptab_list_5<-list()

for (i in 1:3)
{
  st_predic1<- predict(st_empo3,ds[1,], alpha1=0.00001, alpha2=0.001)
  st_predic2 <- as.data.frame(st_predic1[[2]])
  pred_list_5[[i]]<-st_predic1
  proptab_list_5[[i]]<-st_predic2
}

save(pred_list_5,file="pred_list_5_a00001.RData")
save(proptab_list_5,file="proptab_list_5_a00001.RData")

prop_name=read.table("prop_name.txt",header = 1, sep = "\t", check.names = F)

for (j in 1:25)
{
  load(as.character(prop_name$filename[j]))
  
  proptab <-do.call("rbind",proptab_list_5)
  #proptab["EP",]=c(0.00,0.00,0.00,0.00,0.50,0.00,0.00,0.00,0.15,0.00,0.00,0.15,0.00,0.20,0.00,0)#dsC
  #proptab["EP",]=c(0,0,0.2,0.25,0,0,0,0,0,0,0,0,0,0.55,0,0)#dsB
  proptab["EP",]=c(0,0,0,0,0,0,0,0.3,0,0.4,0,0,0.3,0,0,0)#dsA
  
  for (i in 1:16)
  {
    proptab["rmse",i]=sqrt(((proptab[1,i]-proptab[4,i])^2+(proptab[2,i]-proptab[4,i])^2+(proptab[3,i]-proptab[4,i])^2)/3)
  }
  
  proptab["rmse",17]=(cor(as.numeric(proptab[1,1:15]),as.numeric(proptab[4,1:15]))+
                        cor(as.numeric(proptab[2,1:15]),as.numeric(proptab[4,1:15]))+ cor(as.numeric(proptab[3,1:15]),as.numeric(proptab[4,1:15])))/3
  
  write.table(proptab[5,], file = "proptab_mse.csv", sep = ",",row.names = TRUE,col.names = F, append = T)
}

#################################
#Generalization evaluation#----
#################################
load(file="st_empo3.RData")

setwd("./generalization")
load(file="Generalization_sink_otu_index_nor.RData")

pred_list<-list()
proptab_list<-list()
l=1

for (i in 1:nrow(sink_otu_index_nor))
{
  st_predic1<- predict(st_empo3, sink_otu_index_nor[i,], alpha1=0.001, alpha2=0.001)
  st_predic2 <- as.data.frame(st_predic1[[2]])
  st_predic2$sampleID<-rownames(sink_otu_index_nor)[i]
  pred_list[[l]]<-st_predic1
  proptab_list[[l]]<-st_predic2
  l=l+1
}

proptab=read.table("testdata_plus_as.csv", sep=",", header=T, check.names = F)

proptab_T<-subset(proptab,matched%in%"TRUE")
class(proptab_T$predictedratio)<-"double"
plotdata_stprob <- aggregate(proptab_T$predictedratio,by = list(eco_type=proptab_T$empo_3), FUN = function(x) c(mean = mean(x), sd = sd(x),n = length(x)))
plotdata_stprob <- do.call(data.frame, plotdata_stprob)
colnames(plotdata_stprob) <- c("empo_3", "mean", "sd", "n")
write.table(plotdata_stprob, file="plotdata_stprob.csv", sep=",", col.names=T)

limits <- aes(ymax = mean + sd,ymin = mean - sd)
ggplot(data = plotdata_stprob, aes(x = empo_3, y = mean, fill = empo_3))+
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(limits, width = 0.3,color="gray50") +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill='gray96', colour='gray'),
        text=element_text(size=25,color = "black"))+
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_fill_manual(values=c("#9ddfd3", "#17706e", "#edc988", "#c19065", "#19456b", "#62760c", "#9dad7f", "#d6b0b1", 
                             "#8b5e83","#f08a5d", "#a6b1e1", "#c05555",  "#59A5D8", "#07689f", "#90303d"), name="Environmental types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","WWTPs"))+
  scale_x_discrete(labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                            "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                            "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                            "Water (non-saline)","Water (saline)","WWTPs"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust = 1))+
  ylab("Probability")+
  xlab("")

#---correct ration calculation--
correct_ratio<-matrix(nrow=3,ncol=15)
correct_ratio<-data.frame(correct_ratio)
colnames(correct_ratio)<-(colnames(proptab))[2:16]
rownames(correct_ratio)<-c("samplesize","correct_sample_size","correct_ratio")

for (i in 2:16)
{
  correct_ratio[1,i-1]=nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
  correct_ratio[2,i-1]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$matched=="TRUE"),])
  correct_ratio[3,i-1]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$matched=="TRUE"),])/nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
}

write.table(correct_ratio, file = "correct_ratio.csv", append = FALSE, quote= FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)

nrow(proptab[proptab$matched=="TRUE",])#correct sample size for 474 samples,415

#pie chart for the correct ratio
correct_ratio=read.table("correct_ratio.csv", header=T, row.names=1, sep=",",check.names=FALSE)

correct_ratio_t=data.frame(t(correct_ratio),check.names = F)
correct_ratio_t$env_type=rownames(correct_ratio_t)

png("correct_barchart_2.png", res = 600,width = 12, height =2, units = 'in')
ggplot(correct_ratio_t, aes(x="", y = correct_ratio, fill=env_type)) +
  geom_bar(width = 1, stat = "identity")+
  ylim(0,1)+
  scale_fill_manual(values=c("#9ddfd3", "#17706e", "#edc988", "#c19065", "#19456b", "#62760c", "#9dad7f", "#d6b0b1", 
                             "#8b5e83","#f08a5d", "#a6b1e1", "#c05555",  "#59A5D8", "#07689f", "#90303d"), name="Environmental types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","WWTPs"))+
  scale_x_discrete(labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                            "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                            "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                            "Water (non-saline)","Water (saline)","WWTPs"))+
  coord_polar(theta = "y",direction = -1)+
  facet_grid(.~ env_type)+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())
graphics.off()

#Sensitivity and specificity
SS<-matrix(nrow=16,ncol=3)
SS<-data.frame(SS)
rownames(SS)<-(colnames(proptab))[2:17]
colnames(SS)<-c("Sensitivity","Specificity","Sample size")
for (i in 2:17)
{
  SS[i-1,3]=nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
  SS[i-1,1]=nrow(proptab[(proptab$predicted==(colnames(proptab))[i])&(proptab$empo_3==(colnames(proptab))[i]),])/nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
  SS[i-1,2]=nrow(proptab[(proptab$predicted!=(colnames(proptab))[i])&(proptab$empo_3!=(colnames(proptab))[i]),])/nrow(proptab[proptab$empo_3!=(colnames(proptab))[i],])
}
save(SS,file="SS.RData")
write.table(SS, file = "SS.csv", append = FALSE, quote= FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)


APSS<-matrix(nrow=16,ncol=10)
APSS<-data.frame(APSS)
rownames(APSS)<-(colnames(proptab))[2:17]
colnames(APSS)<-c("TP","FP","TN","FN","Precision","Recall","F1","Sensitivity","Specificity","Sample size")
for (i in 2:17)
{
  APSS[i-1,1]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$predicted==(colnames(proptab))[i]),])#TP
  APSS[i-1,2]=nrow(proptab[(proptab$empo_3!=(colnames(proptab))[i])&(proptab$predicted==(colnames(proptab))[i]),])#FP
  APSS[i-1,3]=nrow(proptab[(proptab$empo_3!=(colnames(proptab))[i])&(proptab$predicted!=(colnames(proptab))[i]),])#TN
  APSS[i-1,4]=nrow(proptab[(proptab$empo_3==(colnames(proptab))[i])&(proptab$predicted!=(colnames(proptab))[i]),])#FN
  APSS[i-1,5]=APSS[i-1,1]/(APSS[i-1,1]+ APSS[i-1,2]) #P = TP/(TP+FP)
  APSS[i-1,6]=APSS[i-1,1]/(APSS[i-1,1]+ APSS[i-1,4]) #R = TP/(TP+FN)
  APSS[i-1,7]=2*APSS[i-1,5]* APSS[i-1,6]/( APSS[i-1,5]+ APSS[i-1,6]) # F1 = 2P*R/(P+R)
  APSS[i-1,8]=APSS[i-1,1]/(APSS[i-1,1]+APSS[i-1,4]) #sensitivity
  APSS[i-1,9]=APSS[i-1,3]/(APSS[i-1,2]+APSS[i-1,3]) #specificity
  APSS[i-1,10]=nrow(proptab[proptab$empo_3==(colnames(proptab))[i],])
}
write.table(APSS, file = "APSS.csv", append = FALSE, quote= FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
save(APSS, file = "APSS.RData")

PRF<-matrix()
PRF<-data.frame(PRF)
PRF$macro_precision =mean(APSS$Precision[1:15]) 
PRF$marcro_reall = mean(APSS$Recall[1:15])
PRF$macro_F1 = mean(APSS$F1[1:15])
PRF$micro_precision=sum(APSS$TP[1:15])/sum(sum(APSS$TP[1:15]),sum(APSS$FP[1:15]))
PRF$micro_recall=sum(APSS$TP[1:15])/sum(sum(APSS$TP[1:15]),sum(APSS$FN[1:15]))
PRF$micro_F1=2*PRF$micro_precision*PRF$micro_recall/(PRF$micro_precision+PRF$micro_recall)
write.table(PRF, file = "PRF.csv", append = FALSE, quote= FALSE, sep = ",",eol = "\n", na = "NA", dec = ".", row.names = TRUE,col.names = TRUE)
save(PRF, file = "PRF.RData")

#################################
#Simulated datasets validation#----
#################################
load(file="st_empo3.RData")

setwd("./simulated_dataset_validation")

load("ds1.RData") #simulated datset 1
pred_list_1<- list()
proptab_list_1 <- list()

st_predic1<- predict(st_empo3,ds[1,], alpha1=0.001, alpha2=0.001)
st_predic2 <- as.data.frame(st_predic1[[2]])

pred_list_1[[1]]<-st_predic1
proptab_list_1[[1]]<-st_predic2
save(pred_list_1, file="pred_list_1.RData")
save(proptab_list_1, file="proptab_list_1.RData")

#after generate all 10 proptab_list and pred_list, merge them together
for (i in 1:10)
{
  pred_name=paste("pred_list_",i,".RData",sep=(""))
  load(file=pred_name)
  prop_name=paste("proptab_list_",i,".RData",sep=(""))
  load(file=prop_name)
}

pred_list<-list()
for (i in 1:10)
{
  x_name=paste("pred_list_",i,sep=(""))
  pred_list<-c(pred_list,get(x_name))
}
save(pred_list,file="pred_list.RData")

proptab_list<-list()
for (i in 1:10)
{
  x_name=paste("proptab_list_",i,sep=(""))
  proptab_list<-c(proptab_list,get(x_name))
}
save(proptab_list,file="proptab_list.RData")

load(file="proptab_list.RData")
load(file = "pred_list.RData")
proptab <-do.call("rbind",proptab_list)
rownames(proptab)<-c("ds1","ds2","ds3","ds4","ds5","ds6","ds7","ds8","ds9","ds10")

exp_ra<-matrix(nrow=10,ncol=15)
exp_ra<-data.frame(exp_ra)
colnames(exp_ra)=colnames(exp_ra)[1:15]
rownames(exp_ra)<-c("ds1","ds2","ds3","ds4","ds5","ds6","ds7","ds8","ds9","ds10")
exp_ra[1,]= c(0,0.1,0.3,0,0,0.3,0.2,0,0,0,0,0,0.1,0,0) #a vector of ratio
exp_ra[2,]= c(0,0,0,0,0,0,0.4,0.15,0.2,0,0,0.1,0,0.15,0)
exp_ra[3,]= c(0.1,0,0,0,0,0,0,0,0,0,0,0.15,0.4,0.35,0)
exp_ra[4,]=c(0,0,0,0,0.7,0,0,0,0,0,0,0.1,0,0.2,0)
exp_ra[5,]=c(0,0,0,0,0,0,0,0.3,0,0.4,0,0,0.3,0,0)
exp_ra[6,]=c(0,0,0.2,0.25,0,0,0,0,0,0,0,0,0,0.55,0)
exp_ra[7,]=c(0.1,0,0,0,0,0.3,0.15,0,0,0.45,0,0,0,0,0)
exp_ra[8,]=c(0,0,0,0.5,0,0,0,0,0,0,0,0,0.5,0,0)
exp_ra[9,]=c(0,0,0,0,0.8,0,0,0,0.2,0,0,0,0,0,0)
exp_ra[10,]=c(0,0,0,0,0.5,0,0,0,0.15,0,0,0.15,0,0.2,0)
save(exp_ra,file="exp_ra.RData")

load(file="exp_ra.RData")
std<-matrix(nrow=10,ncol=15)
std<-data.frame(std)
colnames(std)=colnames(proptab)[1:15]
for (i in 1:10)
{
  std[i,]=apply(pred_list[[i]]$draws,2,sd)[1:15]
}

rstd<-matrix(nrow=10,ncol=15)
rstd<-data.frame(rstd)
colnames(rstd)=colnames(proptab)[1:15]
for (i in 1:10)
{
  rstd[i,]=std[i,]/exp_ra[i,]
}

mse<-matrix(nrow=10,ncol=15)
mse<-data.frame(mse)
colnames(mse)=colnames(proptab)[1:15]
for (i in 1:10)
{
  mse[i,]=(proptab[i,1:15]-exp_ra[i,])^2
}


sim_re=matrix(nrow=50,ncol=15)
sim_re<-data.frame(sim_re)
colnames(sim_re)=colnames(proptab)[1:15]
for (i in 1:10)
{
  sim_re[5*(i-1)+1,]=exp_ra[i,]
  sim_re[5*(i-1)+2,]=proptab[i,1:15]
  sim_re[5*(i-1)+3,]=std[i,]
  sim_re[5*(i-1)+4,]=rstd[i,]
  sim_re[5*(i-1)+5,]=mse[i,]
}
write.table(sim_re, file = "sim_re.csv", sep = ",",row.names = TRUE,col.names = TRUE)

#################################################
#WWTPs influent and effluent validation tracking#----
################################################
load(file="empo3.RData")

setwd("./wwtps")
load(file="influent_sink_otu_index_nor.RData")
#load(file="effluent_sink_otu_index_nor.RData")

pred_list<-list()
proptab_list<-list()
l=1

for (i in 1:nrow(sink_otu_index_nor))
{
  st_predic1<- predict(st_empo3, sink_otu_index_nor[i,], alpha1=0.001, alpha2=0.001)
  st_predic2 <- as.data.frame(st_predic1[[2]])
  st_predic2$sampleID<-rownames(sink_otu_index_nor)[i]
  pred_list[[l]]<-st_predic1
  proptab_list[[l]]<-st_predic2
  l=l+1
}

#plot the results
plot_data_Inf_Eff=read.table("plot_data_Inf_Eff_modified.csv", sep=",", header=T)
plot_data_Inf_Eff$env_type<-factor(plot_data_Inf_Eff$env_type,levels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human Fecal",
                                                                       "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                                                                       "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                                                                       "Water (non-saline)","Water (saline)","WWTPs","Unknown"))
plot_data_Inf_Eff$sample_name<-factor(plot_data_Inf_Eff$sample_name, levels=c("SK-IN-1","SK-EFF-1","SK-IN-2","SK-EFF-2",
                                                                              "ST-IN-1","ST-EFF-1","ST-IN-2","ST-EFF-2",
                                                                              "STL-IN-1","STL-EFF-1",
                                                                              "SWH-IN-1","SWH-EFF-1","SWH-IN-2","SWH-EFF-2","SWH-IN-3","SWH-EFF-3",
                                                                              "TP-IN-1","TP-EFF-1","TP-IN-2","TP-EFF-2",
                                                                              "YL-IN-1","YL-EFF-1","YL-IN-2","YL-EFF-2","YL-IN-3","YL-EFF-3"))


plot_data_Inf=subset(plot_data_Inf_Eff, sewage_type=="Inf")
plot_data_Eff=subset(plot_data_Inf_Eff, sewage_type=="Eff")

order_inf = (subset(plot_data_Inf, env_type=="WWTPs")[,7])[order(subset(plot_data_Inf, env_type=="WWTPs")[,3],decreasing=TRUE)]

plot_data_Inf$sample_name<-factor(plot_data_Inf$sample_name, levels=c("STL-IN-1","YL-IN-2","SK-IN-1","ST-IN-2","ST-IN-1",
                                                                      "YL-IN-3","SWH-IN-3","SK-IN-2","SWH-IN-1","SWH-IN-2",
                                                                      "TP-IN-2","YL-IN-1","TP-IN-1"))

plot_data_Eff$sample_name<-factor(plot_data_Eff$sample_name, levels=c("STL-EFF-1","YL-EFF-2","SK-EFF-1","ST-EFF-2","ST-EFF-1",
                                                                      "YL-EFF-3","SWH-EFF-3","SK-EFF-2","SWH-EFF-1","SWH-EFF-2",
                                                                      "TP-EFF-2","YL-EFF-1","TP-EFF-1"))

png("Inf_barplot.png", res = 600,width = 8, height = 6, units = 'in')
ggplot(data = plot_data_Inf, aes(x = sample_name, y = proportion, fill = env_type))+
  geom_bar(stat = "identity", width = 0.7) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill='gray96', colour='gray'),
        text=element_text(size=25,color = "black"))+
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_fill_manual(values=c("#696969","#ffff02","#8b4512","#c98c23","#8ac8ed",
                             "#7cfc02","#006400","#ffa07a","#ee0000","#8FBC8F",
                             "#e46bef","#8a1894","#4169e0","#0a0a9f","#5c6b06",
                             "#d2b48c"), name="Environmental Types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","WWTPs","Unknown"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust = 1))+
  ylab("Predicted source proportion")+
  xlab("Influent sample")
graphics.off()

png("Eff_barplot.png", res = 600,width = 8, height = 6, units = 'in')
ggplot(data = plot_data_Eff, aes(x = sample_name, y = proportion, fill = env_type))+
  geom_bar(stat = "identity", width = 0.7) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill='gray96', colour='gray'),
        text=element_text(size=25,color = "black"))+
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_fill_manual(values=c("#696969","#ffff02","#8b4512","#c98c23","#8ac8ed",
                             "#7cfc02","#006400","#ffa07a","#ee0000","#8FBC8F",
                             "#e46bef","#8a1894","#4169e0","#0a0a9f","#5c6b06",
                             "#d2b48c"), name="Environmental Types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","WWTPs","Unknown"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust = 1))+
  ylab("Predicted source proportion")+
  xlab("Effluent sample")
graphics.off()

#################################
#Marine sediment source tracking#----
#################################
#import the source dataset
load(file="otu_3654.RData")
load(file="st_empo3.RData")
#import real dataset
setwd("./marine_sediment")
real_dataset <- "sediment-cr-97.txt"

sink_otu <- read.table(real_dataset, row.names = 1,  header=T, sep="\t",check.names=FALSE, fill=TRUE)
sink_otu$OTU_ID <- rownames(sink_otu)

#merge the real source dataset and real dataset
library(dplyr)
index <- left_join(otu_3654, sink_otu, by="OTU_ID")
index[is.na(index)] <- 0 #convert the NA into 0
rownames(index) <- index$OTU_ID

sink_otu_index <- data.frame(t(index[,3656:3690]), check.names = F) #choose the sink sample

#converted from relative abundance
sink_otu_index_nor= sink_otu_index/(rowSums(sink_otu_index))
sink_otu_index_nor <-ceiling(sink_otu_index_nor*10000)
#direct normalization by vegan
#wwtps_otu_index_nor=rrarefy(wwtps_otu_index,10000) 

pred_list<-list()
proptab_list<-list()

for (i in 1:nrow(sink_otu_index_nor))
{
  st_predic1<- predict(st_empo3, sink_otu_index_nor[i,], alpha1=0.001, alpha2=0.001)
  st_predic2 <- as.data.frame(st_predic1[[2]])
  st_predic2$sampleID<-rownames(sink_otu_index_nor)[i]
  pred_list[[i]]<-st_predic1
  proptab_list[[i]]<-st_predic2
}

load("sediment_proptab_empo3.RData")
plot_data=proptab

order_sediment=rownames(plot_data[order(-plot_data$`Sediment (saline)`),])
plot_data=melt(plot_data,id.vars=c("sampleID"),variable.name="env_type",value.name="proportion")
plot_data$env_type<-factor(plot_data$env_type,levels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human Fecal",
                                                       "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                                                       "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                                                       "Water (non-saline)","Water (saline)","WWTPs","Unknown"))
plot_data$sampleID=factor(plot_data$sampleID, levels=order_sediment)


#plot for proportional bar
ggplot(data = plot_data, aes(x = sampleID, y = proportion, fill = env_type))+
  geom_bar(stat = "identity", width = 0.7) +
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position="none",
        panel.background = element_rect(fill='gray96', colour='gray'),
        text=element_text(size=25,color = "black"))+
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_fill_manual(values=c("#696969","#ffff02","#8b4512","#c98c23","#8ac8ed",
                             "#7cfc02","#006400","#ffa07a","#ee0000","#8FBC8F",
                             "#e46bef","#8a1894","#4169e0","#0a0a9f","#5c6b06",
                             "#d2b48c"), name="Environmental types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","WWTPs","Unknown"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 45,hjust = 1))+
  ylab("Predicted source proportion")+
  xlab("")

############################
#Global freshwater samples#----
###########################
setwd("./global_freshwater_samples")
load("wwtps_aquatic_541_classification_regroup.RData")
#world map
library(maps)
library(ggplot2)

wwtps_aquatic_1<-wwtps_aquatic[wwtps_aquatic[,1] <=0.01,]
wwtps_aquatic_2<-wwtps_aquatic[wwtps_aquatic[,1] <=0.05 &  wwtps_aquatic[,1] > 0.01,]
wwtps_aquatic_3<-wwtps_aquatic[wwtps_aquatic[,1] <=0.1 &  wwtps_aquatic[,1] > 0.05,]
wwtps_aquatic_4<-wwtps_aquatic[wwtps_aquatic[,1] <=0.4 &  wwtps_aquatic[,1] > 0.1,]

mp<-NULL
worldmap <- borders("world",colour = "gray50")
mp <- ggplot() + worldmap + ylim(-60,90)
mp


map <- mp + 
  geom_point(aes(x = wwtps_aquatic$Longitude, y = wwtps_aquatic$Latitude, color=(wwtps_aquatic$quantile)),position=position_jitter(width=0.1),size=2)+
  ggtitle("")+scale_color_manual(values=c("#78B7C5", "#07689f", "#EBCC2A","#90303d"), name="WWTPs abundance",
                                 labels=c("~ 1%  (29.8%)", "1% ~ 5%  (52.1%)", "5% ~ 10%  (14.4%)", "10% ~ (3.7%)"))+
  geom_point(aes(x = wwtps_aquatic_1$Longitude, y = wwtps_aquatic_1$Latitude),color="#78B7C5",
             position=position_jitter(width=0.1),size=2)+
  geom_point(aes(x = wwtps_aquatic_2$Longitude, y = wwtps_aquatic_2$Latitude),color="#07689f",
             position=position_jitter(width=0.1),size=2)+
  geom_point(aes(x = wwtps_aquatic_3$Longitude, y = wwtps_aquatic_3$Latitude),color="#EBCC2A",
             position=position_jitter(width=0.1),size=2)+
  geom_point(aes(x = wwtps_aquatic_4$Longitude, y = wwtps_aquatic_4$Latitude),color="#90303d",
             position=position_jitter(width=0.1),size=1)+
  ggtitle("")+
  theme_classic(base_size = 16)+
  theme(axis.line = element_blank(),
        axis.title= element_blank(),
        axis.text= element_blank(),
        axis.ticks= element_blank(),
        plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        panel.background = element_rect(fill = NA, color = "black", size = 1),
        legend.position="right", legend.title = element_text(size=12), legend.text  = element_text(size=12))

map


#barchart
wwtps_top <- read.table("fresh_water_all541_wwtps.txt", sep="\t", header = T, quote ="",check.names = F)
plot_data=wwtps_top
order_sediment=(plot_data[order(-plot_data$WWTPs),])[,1]
lable_sediment=(plot_data[order(-plot_data$WWTPs),])[,18]
library(reshape2)
library(ggplot2)
plot_data=melt(plot_data[,1:18],id.vars=c("sampleID","geographic_location"),variable.name="env_type",value.name="proportion")
plot_data$env_type<-factor(plot_data$env_type,levels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                                                       "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                                                       "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                                                       "Water (non-saline)","Water (saline)","WWTPs","Unknown"))
plot_data$env_type<-factor(plot_data$env_type,levels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                                                       "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                                                       "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                                                       "Water (non-saline)","Water (saline)","Unknown","WWTPs"))
plot_data$sampleID=factor(plot_data$sampleID, levels=order_sediment)

ggplot(data = plot_data, aes(x = sampleID, y = proportion, fill = env_type))+
  geom_bar(stat = "identity",width = 0.7) +
  scale_y_continuous(breaks=seq(0,1,by=0.1))+
  scale_fill_manual(values=c("#9ddfd3", "#17706e", "#edc988", "#c19065", "#19456b", "#62760c", "#9dad7f", "#d6b0b1", 
                             "#8b5e83","#f08a5d", "#a6b1e1", "#c05555",  "#59A5D8", "#07689f", "#d8d3cd", "#90303d"), 
                    name="Environmental types",
                    labels=c("Aerosol (non-saline)","Animal corpus","Animal fecal","Human fecal",
                             "Hypersaline (saline)","Plant corpus","Plant rhizosphere","Sediment (non-saline)",
                             "Sediment (saline)","Soil (non-saline)","Surface (non-saline)","Surface (saline)",
                             "Water (non-saline)","Water (saline)","Unknown","WWTPs"))+
  theme(axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background = element_rect(fill='gray96', colour='gray'),
        text=element_text(size=1,color = "black"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),
        plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))+
  ylab("Predicted source proportion")+
  xlab("")
graphics.off()
