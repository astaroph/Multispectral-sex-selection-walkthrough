###First, taking the condensed output from the sensitivity analysis and reshaping it into
## the format needed for the sensitivity plots
library(plyr)
library(tidyr)
library(reshape2)
library(ggplot2)

# sens<-read.csv('Butterflies_Final_BS_sensitivity_analysis_results_error_range_mean_10000_calculated_v3.csv')
# sens<-read.csv('Moth_nocturnalONLY_Final_BS_sensitivity_analysis_results_error_range_mean_10000_calculated_v3.csv')
sens<-read.csv('Geometrids_Final_BS_sensitivity_analysis_results_error_range_mean_10000_calculated.csv')

sens[is.na(sens)] = 0
sens_med<-subset(ddply(sens,.(Value,Metrics),colwise(median)),select=-c(Dataset))
colnames(sens_med)[3:5]<-paste('Median',colnames(sens_med)[3:5],sep='_')
sens_max<-subset(ddply(sens,.(Value,Metrics),colwise(max)),select=-c(Dataset))
colnames(sens_max)[3:5]<-paste('Max',colnames(sens_max)[3:5],sep='_')
trial<-join(sens_med,sens_max,by=c('Value','Metrics'))
trial2<-reshape(trial, idvar = c('Metrics'), timevar = c('Value'),direction = "wide",sep='_')


# write.csv(trial2,'Butterflies_Final_BS_sensitivity_analysis_results_error_range_mean_10000_calculated_v3_reshaped.csv')
# write.csv(trial2,'Moth_nocturnalONLY_Final_BS_sensitivity_analysis_results_error_range_mean_10000_calculated_v3_reshaped.csv')
# write.csv(trial2,'Geometrids_Final_BS_sensitivity_analysis_results_error_range_mean_10000_calculated_reshaped.csv')

# sens<-subset(sens,Max_MDI_stat_error>0)

#butterflies
x_cutoff=170
y_cutoff=25
#nocturnal moths
x_cutoff=185
y_cutoff=40
# geometrids
x_cutoff=194
y_cutoff=86


####Strip-plot showing the % Disparity through time error across all bootstraps for each metric, ordered by decreasing error
ggplot(trial2, aes(y=Max_X95_percent_CI_percentofmean_DTT, x=reorder(Metrics, -Max_X95_percent_CI_percentofmean_DTT))) + 
  # geom_violin(alpha=0.4,draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_point(size=1,alpha=0.5)+
  xlab('Metrics (decreasing sensitivity)') + ylab('Maximum DTT mean CI %of mean') +
  geom_hline(yintercept=x_cutoff,linetype='dashed',color='red')+
  geom_hline(yintercept=0,linetype='dashed',color='black')+
  theme(text = element_text(size=14))+
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))+
  theme(axis.text.x = element_blank())+
  facet_grid(.~.)


####Stripplot showing the % mean value error across all bootstraps for each metric, ordered by decreasing error
ggplot(trial2, aes(y=Max_X95_percent_CI_percentofmean_Mean, x=reorder(Metrics, -Max_X95_percent_CI_percentofmean_Mean))) + 
  # geom_violin(alpha=0.4,draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_point(size=1,alpha=0.5)+
  xlab('Metrics (decreasing sensitivity)') + ylab('Maximum Mean value sensitivity') +
  geom_hline(yintercept=y_cutoff,linetype='dashed',color='red')+
  geom_hline(yintercept=0,linetype='dashed',color='black')+
  theme(text = element_text(size=14))+
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))+
  theme(axis.text.x = element_blank())+
  facet_grid(.~.)


####Trying them combined
pdf("Figures/Butterflies Final scatterplot and histogram DTT and mean value sensitivity threshold fixed maxima.pdf",width=11.042,height=7.396)
pdf("Figures/Nocturnal moths Final scatterplot and histogram DTT and mean value sensitivity threshold fixed maxima.pdf",width=11.042,height=7.396)
pdf("Figures/Geometridae Final scatterplot and histogram DTT and mean value sensitivity threshold fixed maxima.pdf",width=11.042,height=7.396)
ggplot(trial2, aes(y=Max_X95_percent_CI_percentofmean_Mean, x=(Max_X95_percent_CI_percentofmean_DTT))) + 
  # geom_violin(alpha=0.4,draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_rect(aes(xmin=0, xmax=x_cutoff, ymin=0,
                ymax=y_cutoff),fill='#BEC0E4',alpha=0.02)+
  geom_rect(aes(xmin=0, xmax=x_cutoff, ymin=y_cutoff,
                # ymax=max(Max_X95_percent_CI_percentofmean_Mean)),fill='#FEDADA',alpha=0.02)+
                ymax=300),fill='#FEDADA',alpha=0.02)+
  
  # geom_rect(aes(xmin=x_cutoff, xmax=max(Max_X95_percent_CI_percentofmean_DTT), ymin=0,
  geom_rect(aes(xmin=x_cutoff, xmax=2000, ymin=0,
                ymax=y_cutoff),fill='#FEDADA',alpha=0.02)+
  # geom_rect(aes(xmin=x_cutoff, xmax=max(Max_X95_percent_CI_percentofmean_DTT), ymin=y_cutoff,
  geom_rect(aes(xmin=x_cutoff, xmax=2000, ymin=y_cutoff,
                # ymax=max(Max_X95_percent_CI_percentofmean_Mean)),fill='#FEDADA',alpha=0.02)+
                ymax=300),fill='#FEDADA',alpha=0.02)+
  
  geom_point(size=1.5,alpha=0.5)+
  xlab('Max. 95% CI range (% of mean DTT)') + ylab('Max. 95% CI range (% of mean value)') +
  geom_hline(yintercept=y_cutoff,linetype='dashed',color='red',alpha=0.5)+
  geom_hline(yintercept=0.002,linetype='dashed',color='black',alpha=0.5)+
  geom_vline(xintercept=x_cutoff,linetype='dashed',color='red',alpha=0.5)+
  geom_vline(xintercept=0.002,linetype='dashed',color='black',alpha=0.5)+
  theme(text = element_text(size=30,family = "serif"))+
  theme(axis.text = element_text(size=26,family = "serif",color='black'))+
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_y_continuous(limits=c(0,2000))+
  scale_y_continuous(limits=c(0,300))+
  # theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5))+
  facet_grid(.~.)
dev.off()

###Bar histograms for subplots
ggplot(trial2, aes(x=(Max_X95_percent_CI_percentofmean_Mean))) + 
  # geom_violin(alpha=0.4,draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_bar()+
  # scale_x_binned(n.breaks=64,trans="sqrt")+
  scale_x_binned(n.breaks=30,limits=c(0,300))+
    xlab('Metrics (decreasing sensitivity)') +
  geom_vline(xintercept=y_cutoff,linetype='dashed',color='red')+
  theme(text = element_text(size=14))+
  theme(text = element_blank())+
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))+
  theme(axis.text.x = element_blank())+
  coord_flip()

ggplot(trial2, aes(x=(Max_X95_percent_CI_percentofmean_DTT))) + 
  # geom_violin(alpha=0.4,draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_bar()+
  # scale_x_binned(n.breaks=64,trans="sqrt")+
  scale_x_binned(n.breaks=60,limits=c(0,2000))+
  xlab('Metrics (decreasing sensitivity)') +
  geom_vline(xintercept=x_cutoff,linetype='dashed',color='red')+
  # theme(text = element_text(size=14))+
  theme(text = element_blank())+
  theme(axis.line = element_line(colour = "black"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  # scale_x_sqrt(breaks = c(0,25,50,100,250,500,1000,1500))+
  # theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust=1))+
  theme(axis.text.x = element_blank())
