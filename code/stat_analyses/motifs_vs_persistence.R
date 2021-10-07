library(lmerTest)
library(MuMIn)
library(GMCM) # colSds is oddly not in standard R

# Read in networks, collect into one data frame
fulldata=matrix(nrow=0,ncol=31)
allTLs=matrix(nrow=0,ncol=3) # Storage for TL, degree
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
			allTLs=rbind(allTLs,cbind(rep(S,nrow(data)),rep(C,nrow(data)),data$STL))
		}
		alldata=cbind(exttime,allroles)
		fulldata=rbind(fulldata,alldata)
}}

fulldata=as.data.frame(fulldata)
fulldata$rando=paste0(fulldata$S,':',fulldata$C,':',fulldata$Network)
fulldata$rando=as.factor(fulldata$rando)
fulldata$globo=paste0(fulldata$S,':',fulldata$C)
fulldata$globo=as.factor(fulldata$globo) # Global-scale random effect

# Models with counts (one per motif)
app_count=with(fulldata,lmer(Mean_order~cS5+(1|globo)))
dir_count=with(fulldata,lmer(Mean_order~cS4+(1|globo)))
omni_count=with(fulldata,lmer(Mean_order~cS2+(1|globo)))
chain_count=with(fulldata,lmer(Mean_order~cS1+(1|globo)))
other_count=with(fulldata,lmer(Mean_order~cother+(1|globo)))

# Models with frequencies
freq=as.data.frame(fulldata[,6:10]/rowSums(fulldata[,6:10]))
app_freq=with(freq,lmer(fulldata$Mean_order~cS5+(1|fulldata$globo)))
dir_freq=with(freq,lmer(fulldata$Mean_order~cS4+(1|fulldata$globo)))
omni_freq=with(freq,lmer(fulldata$Mean_order~cS2+(1|fulldata$globo)))
chain_freq=with(freq,lmer(fulldata$Mean_order~cS1+(1|fulldata$globo)))
other_freq=with(freq,lmer(fulldata$Mean_order~cother+(1|fulldata$globo)))


# Models with Z-scores
app_Z=with(fulldata,lmer(Mean_order~zS5+(1|globo)))
dir_Z=with(fulldata,lmer(Mean_order~zS4+(1|globo)))
omni_Z=with(fulldata,lmer(Mean_order~zS2+(1|globo)))
chain_Z=with(fulldata,lmer(Mean_order~zS1+(1|globo)))
other_Z=with(fulldata,lmer(Mean_order~zother+(1|globo)))


