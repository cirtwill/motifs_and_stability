library(lmerTest)
fulldata=matrix(nrow=0,ncol=20)
for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)
		exttime=matrix(nrow=0,ncol=5)
		allroles=matrix(nrow=0,ncol=15)
		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'_motifs.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)

			roles=data[,(S+2):(S+16)] # Species IDs, then extinction order for each removal, then roles
			allroles=rbind(allroles,roles)

			# Mean extinction order vs. roles, should also store strength of correlation.
			extcorr=mean(cor(data[,2:(S+1)]))
			mean_ext=rowSums(data[,2:(S+1)])/S
			extdat=cbind(rep(S,S),rep(C,S),rep(i,S),rep(extcorr,S),mean_ext)
			colnames(extdat)=c("S","C","Network","Extinction_correlation","Mean_order")
			exttime=rbind(exttime,extdat)
		}
		alldata=cbind(exttime,allroles)
		fulldata=rbind(fulldata,alldata)
}}

fulldata=as.data.frame(fulldata)

results=matrix(nrow=13,ncol=5)
for(i in 1:13){
	motif=colnames(fulldata)[i+7]
	raw_lm=lm(fulldata[,i+7]~fulldata$STL)
	rawco=summary(raw_lm)$coefficients
	prop=fulldata[,i+7]/rowSums(fulldata[,8:20])
	prop_lm=lm(prop~fulldata$STL)
	propco=summary(prop_lm)$coefficients
	results[i,]=c(motif,rawco[2,1],rawco[2,4],propco[2,1],propco[2,4])
}
colnames(results)=c('Motif','Raw_beta','Raw_p','Prop_beta','Prop_p')
write.table(results,file='../../data/summaries/motifs_vs_degree.tsv',sep='\t')


