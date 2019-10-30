library(lme4)
library(lmerTest)

extorders=matrix(nrow=0,ncol=5)
colnames(extorders)=c("S","C","Network","Extinction_correlation","Mean_order")
for(S in seq(50,100,10)){
	# dir.create(paste0('../../data/summaries/extorder_perms/',as.character(S)))
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		# dir.create(paste0('../../data/summaries/extorder_perms/',as.character(S),'/',as.character(C)))
		print(C)
		# allroles=matrix(nrow=0,ncol=30)

		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)

			# roles=data[,(S+2):ncol(data)] # Species IDs, then extinction order for each removal, then roles
			# allroles=rbind(allroles,roles)

			# Mean extinction order vs. roles, should also store strength of correlation.
			extcorr=mean(cor(data[,2:(S+1)]))
			mean_ext=rowSums(data[,2:(S+1)])/S
			extdat=cbind(rep(S,S),rep(C,S),rep(i,S),rep(extcorr,S),mean_ext)
			colnames(extdat)=c("S","C","Network","Extinction_correlation","Mean_order")
			extorders=rbind(extorders,extdat)
		}
}}
		extorders=as.data.frame(extorders)
# 
extorders=as.data.frame(extorders)
extorders$rando=paste0(as.character(extorders$S),as.character(extorders$C),as.character(extorders$Network))
mod=lmer(Extinction_correlation~S*C+(1|rando),data=extorders)



