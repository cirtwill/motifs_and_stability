library(lmerTest)
# Testing: Do species with similar mean extinction orders have similar roles?
# 90 and 100-species ones don't run here, need to be on the mac

extorders=matrix(nrow=0,ncol=5)
colnames(extorders)=c("S","C","Network","Extinction_correlation","Mean_order")
allroles=matrix(nrow=0,ncol=30)
for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)

		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'_motifs.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)

			roles=data[,(S+2):ncol(data)] # Species IDs, then extinction order for each removal, then roles
			allroles=rbind(allroles,roles)

			# Mean extinction order vs. roles, should also store strength of correlation.
			extcorr=mean(cor(data[,2:(S+1)]))
			mean_ext=rowSums(data[,2:(S+1)])/S
			extdat=cbind(rep(S,S),rep(C,S),rep(i,S),rep(extcorr,S),mean_ext)
			colnames(extdat)=c("S","C","Network","Extinction_correlation","Mean_order")
			extorders=rbind(extorders,extdat)
		}
}}
extorders=as.data.frame(cbind(extorders,allroles))
extorders$rando=paste0(as.character(extorders$S),as.character(extorders$C),as.character(extorders$Network))
# Most interested in motifs S1, S2, S4, S5

# Still not sure about how the random effects should work

# Stable motifs
S1=lmer(Mean_order~zS1+(1|S:C)+(1|rando),data=extorders)
S2=lmer(Mean_order~zS2+(1|S:C)+(1|rando),data=extorders)
S4=lmer(Mean_order~zS4+(1|S:C)+(1|rando),data=extorders)
S5=lmer(Mean_order~zS5+(1|S:C)+(1|rando),data=extorders)

# Unstable motifs
S3=lmer(Mean_order~zS3+(1|S:C)+(1|rando),data=extorders)
# Effect of loop should depend on web conditions
D1=lmer(Mean_order~zD1+(1|S:C)+(1|rando),data=extorders,
	control=lmerControl(optimizer="bobyqa"))
# Failure to converge.

D2=lmer(Mean_order~zD2+(1|S:C)+(1|rando),data=extorders)
D3=lmer(Mean_order~zD3+(1|S:C)+(1|rando),data=extorders)
D4=lmer(Mean_order~zD4+(1|S:C)+(1|rando),data=extorders)
D5=lmer(Mean_order~zD5+(1|S:C)+(1|rando),data=extorders)
D6=lmer(Mean_order~zD6+(1|S:C)+(1|rando),data=extorders)
D7=lmer(Mean_order~zD7+(1|S:C)+(1|rando),data=extorders)
# Effects of all two-ways should depend on web conditions
# Maybe I should add random effects somehow? Or plot the actual data?

deg=lmer(Mean_order~Degree+(1|S:C)+(1|rando),data=extorders)
TL=lmer(Mean_order~STL+(1|S:C)+(1|rando),data=extorders)
# TL there's still quite a few infinite trophic levels kicking around.  

write.table(summary(S1)$coef,file='tables/S1_lmer_table.tsv',sep='\t')
write.table(summary(S2)$coef,file='tables/S2_lmer_table.tsv',sep='\t')
write.table(summary(S4)$coef,file='tables/S4_lmer_table.tsv',sep='\t')
write.table(summary(S5)$coef,file='tables/S5_lmer_table.tsv',sep='\t')
write.table(summary(S3)$coef,file='tables/S3_lmer_table.tsv',sep='\t')
write.table(summary(D1)$coef,file='tables/D1_lmer_table.tsv',sep='\t')
write.table(summary(D2)$coef,file='tables/D2_lmer_table.tsv',sep='\t')
write.table(summary(D3)$coef,file='tables/D3_lmer_table.tsv',sep='\t')
write.table(summary(D4)$coef,file='tables/D4_lmer_table.tsv',sep='\t')
write.table(summary(D5)$coef,file='tables/D5_lmer_table.tsv',sep='\t')
write.table(summary(D6)$coef,file='tables/D6_lmer_table.tsv',sep='\t')
write.table(summary(D7)$coef,file='tables/D7_lmer_table.tsv',sep='\t')

