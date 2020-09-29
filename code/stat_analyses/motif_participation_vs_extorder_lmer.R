library(lmerTest)
library(MuMIn)


# Subset species like with PCA axes. Then fit one big model, simplify.
# Testing: Do species with similar mean extinction orders have similar roles?
# No convergence??
bins=matrix(nrow=0,ncol=29)
fulldata=matrix(nrow=0,ncol=31)
# extorders=matrix(nrow=0,ncol=5)
# colnames(extorders)=c("S","C","Network","Extinction_correlation","Mean_order")
# allroles=matrix(nrow=0,ncol=30)
allTLs=matrix(nrow=0,ncol=3)
for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)
		exttime=matrix(nrow=0,ncol=5)
		allroles=matrix(nrow=0,ncol=26)
		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'_motifs.tsv',sep='')
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
		r=1
		for(i in 1:(nrow(alldata)/100)){
			rmax=r+99
			subset=alldata[r:rmax,]
			binres=c(colMeans(subset[,1:2]),colMeans(subset[,5:31]))
			bins=rbind(bins,binres)
			r=r+100
		}
}}

bins=as.data.frame(bins)
fulldata=as.data.frame(fulldata)
fulldata$rando=paste0(fulldata$S,':',fulldata$C,':',fulldata$Network)
fulldata$rando=as.factor(fulldata$rando)

fullmod=lmer(Mean_order~cS1+cS2+cS3+cS4+cS5+cD1+cD2+cD3+cD4+cD5+cD6+cD7+cD8+(1|rando),data=fulldata,na.action=na.fail)
write.table(fullmod$coefficients, file='../../data/summaries/full_lm_coefficients.tsv',sep='\t')
fulldredge=dredge(fullmod)
redmod=
write.table(redmod$coefficients,file='../../data/summaries/reduced_lm_coefficients.tsv',sep='\t')
