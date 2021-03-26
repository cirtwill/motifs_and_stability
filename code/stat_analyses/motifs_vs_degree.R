library(lmerTest)
fulldata=matrix(nrow=0,ncol=20)
Zdata=matrix(nrow=0,ncol=20)
for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)
		exttime=matrix(nrow=0,ncol=5)
		allroles=matrix(nrow=0,ncol=15)
		allZs=matrix(nrow=0,ncol=15)
		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'_motifs.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)

			roles=data[,(S+2):(S+16)] # Species IDs, then extinction order for each removal, then roles
			Zs=cbind(data[,(S+2):(S+3)],data[,(S+17):(S+29)])
			allroles=rbind(allroles,roles)
			allZs=rbind(allZs,Zs)

			# Mean extinction order vs. roles, should also store strength of correlation.
			extcorr=mean(cor(data[,2:(S+1)]))
			mean_ext=rowSums(data[,2:(S+1)])/S
			extdat=cbind(rep(S,S),rep(C,S),rep(i,S),rep(extcorr,S),mean_ext)
			colnames(extdat)=c("S","C","Network","Extinction_correlation","Mean_order")
			exttime=rbind(exttime,extdat)
		}
		alldata=cbind(exttime,allroles)
		fulldata=rbind(fulldata,alldata)

		Zex=cbind(exttime,allZs)
		Zdata=rbind(Zdata,Zex)
}}

fulldata=as.data.frame(fulldata)
Zdata=as.data.frame(Zdata)

Tresults=matrix(nrow=13,ncol=13)
for(i in 1:13){
	motif=colnames(fulldata)[i+7]
	raw_lm=lm(fulldata[,i+7]~fulldata$STL)
	rawco=summary(raw_lm)$coefficients
	prop=fulldata[,i+7]/rowSums(fulldata[,8:20])
	prop_lm=lm(prop~fulldata$STL)
	propco=summary(prop_lm)$coefficients
	zlm=lm(Zdata[,i+7]~fulldata$STL)
	zco=summary(zlm)$coefficients
	Tresults[i,]=c(motif,rawco[1,1],rawco[2,1],rawco[2,2],rawco[2,4],propco[1,1],propco[2,1],propco[2,2],propco[2,4],zco[1,1],zco[2,1],zco[2,2],zco[2,4])
}
colnames(Tresults)=c('Motif','Raw_intercept','Raw_beta','Raw_var','Raw_p','Prop_intercept','Prop_beta','Prop_var','Prop_p','Z_intercept','Z_beta','Z_var','Z_p')
write.table(Tresults,file='../../data/summaries/motifs_vs_STL.tsv',sep='\t')


Dresults=matrix(nrow=13,ncol=13)
for(i in 1:13){
	motif=colnames(fulldata)[i+7]
	raw_lm=lm(fulldata[,i+7]~fulldata$Degree)
	rawco=summary(raw_lm)$coefficients
	prop=fulldata[,i+7]/rowSums(fulldata[,8:20])
	prop_lm=lm(prop~fulldata$Degree)
	propco=summary(prop_lm)$coefficients
	zlm=lm(Zdata[,i+7]~fulldata$Degree)
	zco=summary(zlm)$coefficients
	Dresults[i,]=c(motif,rawco[1,1],rawco[2,1],rawco[2,2],rawco[2,4],propco[1,1],propco[2,1],propco[2,2],propco[2,4],zco[1,1],zco[2,1],zco[2,2],zco[2,4])
}
colnames(Dresults)=c('Motif','Raw_intercept','Raw_beta','Raw_var','Raw_p','Prop_intercept','Prop_beta','Prop_var','Prop_p','Z_intercept','Z_beta','Z_var','Z_p')
write.table(Dresults,file='../../data/summaries/motifs_vs_degree.tsv',sep='\t')


Deg_pers=with(fulldata,lmer(Mean_order~Degree+(1|S:C)))
STL_pers=with(fulldata,lmer(Mean_order~STL+(1|S:C)))

