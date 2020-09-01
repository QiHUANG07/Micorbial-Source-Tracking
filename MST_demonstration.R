##Microbial source tracking

#R codes for data analyses and plotting
#including cross-validation of the 3654-sample source dataset and marine sediment source tracking
#NOTE: the R script below has been simplified for illurstration

###################
#library preparation#
###################
install.packages("ggplot2")
install.packages("reshape")
install.packages("yardstick")
install.packages("dplyr")
library(ggplot2)
library(reshape)
library(yardstick)
library(dplyr)

source('./src/SourceTracker.r') #load SourceTracker

################################
#Leave-one-out cross-validation#
################################

inte_otu_abu = read.table(file="3654_source_dataset.txt", header=T, row.names=1, sep="\t",check.names=FALSE, quote="")
inte_otu_abu<-as.matrix(inte_otu_abu) 
class(inte_otu_abu)<-"numeric"

design = read.table("Training_design.csv",  header=T, sep=",",check.names=FALSE)

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

#pack the trained classifier
st_empo3<-sourcetracker(inte_otu_abu, design$empo_3)

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


#################################
#Marine sediment source tracking#
#################################
#import the source dataset
load(file="/media/HQ/EMP_data/model_new/Train_dataset/otu_3654.RData")

#import real dataset
real_dataset <- "sediment-cr-97.txt"

sink_otu <- read.table(real_dataset, row.names = 1,  header=T, sep="\t",check.names=FALSE, fill=TRUE)
sink_otu$OTU_ID <- rownames(sink_otu)

#merge the real source dataset and real dataset
library(dplyr)
index <- left_join(otu_3654, sink_otu, by="OTU_ID")
index[is.na(index)] <- 0 #convert the NA into 0
rownames(index) <- index$OTU_ID

sink_otu_index <- data.frame(t(index[,3656:3690]), check.names = F) #choose the sink sample
sink_otu_index_nor= sink_otu_index/(rowSums(sink_otu_index))

sink_otu_index_nor <-ceiling(sink_otu_index_nor*10000)

source("./src/SourceTracker.r")
load(file="/media/HQ/EMP_data/model_new/Train_dataset/st_empo3.RData")

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

plot_data=proptab[-5,]
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
