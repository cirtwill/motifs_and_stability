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
			datafile=paste0('../../data/roles/matched_to_extorder_condensed/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
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

# # Section 1.2 : Does persistence relate to particular motifs?


# Models with counts (one per motif)
all_count=with(fulldata,lmer(Mean_order~cS5+cS4+cS2+cS1+cother+(1|globo)))

# Models with frequencies
freq=as.data.frame(fulldata[,6:10]/rowSums(fulldata[,6:10]))
#non-independent, so can only include 4/5 motifs...
all_freq=with(freq,lmer(fulldata$Mean_order~cS5+cS4+cS2+cS1+(1|fulldata$globo)))

# Models with Z-scores
all_Z=with(fulldata,lmer(Mean_order~zS5+zS4+zS2+zS1+zother+(1|globo)))


# # Section 2.1 : Does participation relate to degree, TL?
# Again dropping cother from frequency models
degmotifs_count=with(fulldata,lmer(allTLs$Deg~cS5+cS4+cS2+cS1+cother+(1|globo)))
write.table(summary(degmotifs_count)$coefficients,file='../../data/summaries/motifs_Count_Deg.tsv')
degmotifs_freq=with(freq,lmer(allTLs$Deg~cS5+cS4+cS2+cS1+(1|fulldata$globo)))
write.table(summary(degmotifs_freq)$coefficients,file='../../data/summaries/motifs_Freq_Deg.tsv')
degmotifs_Z=with(fulldata,lmer(allTLs$Deg~zS5+zS4+zS2+zS1+zother+(1|globo)))
write.table(summary(degmotifs_Z)$coefficients,file='../../data/summaries/motifs_Z_Deg.tsv')

TLmotifs_count=with(fulldata,lmer(allTLs$STL~cS5+cS4+cS2+cS1+cother+(1|globo)))
write.table(summary(TLmotifs_count)$coefficients,file='../../data/summaries/motifs_Count_TL.tsv')
TLmotifs_freq=with(freq,lmer(allTLs$STL~cS5+cS4+cS2+cS1+(1|fulldata$globo)))
write.table(summary(TLmotifs_freq)$coefficients,file='../../data/summaries/motifs_Freq_TL.tsv')
TLmotifs_Z=with(fulldata,lmer(allTLs$STL~zS5+zS4+zS2+zS1+zother+(1|globo)))
write.table(summary(TLmotifs_Z)$coefficients,file='../../data/summaries/motifs_Z_TL.tsv')


# # Section 2.2: Does persistence relate to degree, TL?

per_degTL=with(fulldata,lmer(Mean_order~allTLs$Deg*allTLs$STL+(1|globo)))



