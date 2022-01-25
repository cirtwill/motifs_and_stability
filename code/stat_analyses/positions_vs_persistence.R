library(lmerTest)
library(MuMIn)
library(GMCM) # colSds is oddly not in standard R

# Read in networks, collect into one data frame
fulldata=matrix(nrow=0,ncol=31)
allTLs=matrix(nrow=0,ncol=4) # Storage for TL, degree
for(S in seq(50,100,10)){
  print(S)
  for(C in seq(0.02,0.2,0.02)){
    print(C)
    exttime=matrix(nrow=0,ncol=5)
    allroles=matrix(nrow=0,ncol=26)
    for(i in seq(0,99)){
      datafile=paste0('../../data/positions/matched_to_extorder_condensed/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
      data=read.table(datafile,sep='\t',header=TRUE)

      roles=data[,(S+4):ncol(data)] # Species IDs, then extinction order for each removal, then roles
      allroles=rbind(allroles,roles)

      # Mean extinction order vs. roles, should also store strength of correlation.
      extcorr=mean(cor(data[,2:(S+1)]))
      mean_ext=rowSums(data[,2:(S+1)])/S
      extdat=cbind(rep(S,S),rep(C,S),rep(i,S),rep(extcorr,S),mean_ext)
      colnames(extdat)=c("S","C","Network","Extinction_correlation","Mean_order")
      exttime=rbind(exttime,extdat)
      allTLs=rbind(allTLs,cbind(rep(S,nrow(data)),rep(C,nrow(data)),data$STL,data$Degree))
    }
    alldata=cbind(exttime,allroles)
    fulldata=rbind(fulldata,alldata)
}}
colnames(allTLs)=c("S","C","STL","Deg")
allTLs=as.data.frame(allTLs)

fulldata=as.data.frame(fulldata)
fulldata$rando=paste0(fulldata$S,':',fulldata$C,':',fulldata$Network)
fulldata$rando=as.factor(fulldata$rando)
fulldata$globo=paste0(fulldata$S,':',fulldata$C)
fulldata$globo=as.factor(fulldata$globo) # Global-scale random effect

allTLs$globo=fulldata$globo

# # Section 1.2 : Does persistence relate to particular positions?
# # Still needs fig.

# Models with counts (one per motif)
all_count=with(fulldata,lmer(Mean_order~
  c12.0.1+c12.1.0+c12.1.1+
  c36.0.1+c36.2.0+
  c38.0.2+c38.1.1+c38.2.0+
  c6.0.2+c6.1.0
  +cother+(1|globo)+(1|rando)))
write.table(summary(all_count)$coefficients,file='../../data/summaries/persistence_countpositions.tsv')

# Models with frequencies
freq=as.data.frame(fulldata[,6:16]/rowSums(fulldata[,6:16]))
#non-independent, so can only include 4/5 motifs...
all_freq=with(freq,lmer(fulldata$Mean_order~c12.0.1+c12.1.0+c12.1.1+
  c36.0.1+c36.2.0+c38.0.2+c38.1.1+c38.2.0+
  c6.0.2+c6.1.0+(1|fulldata$globo)+(1|fulldata$rando)))
write.table(summary(all_freq)$coefficients,file='../../data/summaries/persistence_freqpositions.tsv')

# Models with Z-scores
all_Z=with(fulldata,lmer(Mean_order~  
  z12.0.1+z12.1.0+z12.1.1+
  z36.0.1+z36.2.0+
  z38.0.2+z38.1.1+z38.2.0+
  z6.0.2+z6.1.0
+zother+(1|globo)+(1|rando)))
write.table(summary(all_Z)$coefficients,file='../../data/summaries/persistence_Zpositions.tsv')


# # Section 2.1 : Does participation relate to degree, TL?
# Again dropping cother from frequency models
# Not converging with random effect of network....
degmotifs_count=with(fulldata,lmer(allTLs$Deg~
    c12.0.1+c12.1.0+c12.1.1+
  c36.0.1+c36.2.0+
  c38.0.2+c38.1.1+c38.2.0+
  c6.0.2+c6.1.0
  +cother+(1|globo)+(1|rando),control=lmerControl(optimizer="bobyqa")))
write.table(summary(degmotifs_count)$coefficients,file='../../data/summaries/positions_Count_Deg.tsv')
degmotifs_freq=with(freq,lmer(allTLs$Deg~
    c12.0.1+c12.1.0+c12.1.1+
  c36.0.1+c36.2.0+
  c38.0.2+c38.1.1+c38.2.0+
  c6.0.2+c6.1.0+(1|fulldata$globo)+(1|fulldata$rando),control=lmerControl(optimizer="bobyqa")))
write.table(summary(degmotifs_freq)$coefficients,file='../../data/summaries/positions_Freq_Deg.tsv')
degmotifs_Z=with(fulldata,lmer(allTLs$Deg~
  z12.0.1+z12.1.0+z12.1.1+
  z36.0.1+z36.2.0+
  z38.0.2+z38.1.1+z38.2.0+
  z6.0.2+z6.1.0+zother+(1|globo)+(1|rando),control=lmerControl(optimizer="bobyqa")))
write.table(summary(degmotifs_Z)$coefficients,file='../../data/summaries/positions_Z_Deg.tsv')

# Convergence OK
TLmotifs_count=with(fulldata,lmer(allTLs$STL~c12.0.1+c12.1.0+c12.1.1+
  c36.0.1+c36.2.0+
  c38.0.2+c38.1.1+c38.2.0+
  c6.0.2+c6.1.0+cother+(1|globo)+(1|rando)))
write.table(summary(TLmotifs_count)$coefficients,file='../../data/summaries/positions_Count_TL.tsv')
TLmotifs_freq=with(freq,lmer(allTLs$STL~c12.0.1+c12.1.0+c12.1.1+
  c36.0.1+c36.2.0+
  c38.0.2+c38.1.1+c38.2.0+
  c6.0.2+c6.1.0+(1|fulldata$globo)+(1|fulldata$rando)))
write.table(summary(TLmotifs_freq)$coefficients,file='../../data/summaries/positions_Freq_TL.tsv')
TLmotifs_Z=with(fulldata,lmer(allTLs$STL~z12.0.1+z12.1.0+z12.1.1+
  z36.0.1+z36.2.0+z38.0.2+z38.1.1+z38.2.0+
  z6.0.2+z6.1.0+zother+(1|globo)+(1|rando)))
write.table(summary(TLmotifs_Z)$coefficients,file='../../data/summaries/positions_Z_TL.tsv')



# Table of means for plotting
cmeans=colMeans(fulldata[,6:16])
write.table(as.matrix(cmeans),file='../../data/summaries/mean_positions_count.tsv',sep='\t')

fmeans=colMeans(freq)
write.table(as.matrix(fmeans),file='../../data/summaries/mean_positions_freq.tsv',sep='\t')

zmeans=colMeans(fulldata[,17:27])

save.image('position_tests.Rdata')

