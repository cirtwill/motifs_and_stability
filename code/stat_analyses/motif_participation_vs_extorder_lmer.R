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

extorders$Total_stable=extorders$S1+extorders$S2+extorders$S4+extorders$S5
extorders$Prop_stable=extorders$Total_stable/rowSums(extorders[,6:18])

# Perhaps \% stable vs. unstable would be more efficient?
propmod=lmer(Mean_order~scale(Prop_stable)*scale(S)*scale(C)+(1|rando),data=extorders)
propmod_2=lmer(Mean_order~scale(Prop_stable)*scale(S)+scale(Prop_stable)*scale(C)+scale(S):scale(C)+(1|rando),data=extorders)
# So, higher proportion is longer time to extinction, weaker in higher-C webs.

# Stable motifs
S1=lmer(Mean_order~scale(S1)*scale(S)*scale(C)+(1|rando),data=extorders)
S2=lmer(Mean_order~scale(S2)*scale(S)*scale(C)+(1|rando),data=extorders)
S2_2=lmer(Mean_order~scale(S2)*scale(S)+scale(S2)*scale(C)+scale(S):scale(C)+(1|rando),data=extorders)
S4=lmer(Mean_order~scale(S4)*scale(S)*scale(C)+(1|rando),data=extorders)
S5=lmer(Mean_order~scale(S5)*scale(S)*scale(C)+(1|rando),data=extorders)

# Unstable motifs
S3=lmer(Mean_order~scale(S3)*scale(S)*scale(C)+(1|rando),data=extorders)
S3_2=lmer(Mean_order~scale(S3)*scale(S)+scale(S3)*scale(C)+scale(S):scale(C)+(1|rando),data=extorders)
S3_3=lmer(Mean_order~scale(S3)*scale(C)+scale(S)*scale(C)+(1|rando),data=extorders)
# Effect of loop should depend on web conditions
D1=lmer(Mean_order~scale(D1)*scale(S)*scale(C)+(1|rando),data=extorders)
D2=lmer(Mean_order~scale(D2)*scale(S)*scale(C)+(1|rando),data=extorders)
D3=lmer(Mean_order~scale(D3)*scale(S)*scale(C)+(1|rando),data=extorders)
D4=lmer(Mean_order~scale(D4)*scale(S)*scale(C)+(1|rando),data=extorders)
D4_2=lmer(Mean_order~scale(D4)*scale(S)+scale(D4)*scale(C)+scale(S):scale(C)+(1|rando),data=extorders)
D5=lmer(Mean_order~scale(D5)*scale(S)*scale(C)+(1|rando),data=extorders)
D6=lmer(Mean_order~scale(D6)*scale(S)*scale(C)+(1|rando),data=extorders)
D7=lmer(Mean_order~scale(D7)*scale(S)*scale(C)+(1|rando),data=extorders)
# Effects of all two-ways should depend on web conditions
# Maybe I should add random effects somehow? Or plot the actual data?

scaleS=scale(extorders$S)
scaleC=scale(extorders$C)
scaleS1=scale(extorders$S1)
scaleS2=scale(extorders$S2)
scaleS3=scale(extorders$S3)
scaleS4=scale(extorders$S4)
scaleS5=scale(extorders$S5)
scaleD1=scale(extorders$D1)
scaleD2=scale(extorders$D2)
scaleD3=scale(extorders$D3)
scaleD4=scale(extorders$D4)
scaleD5=scale(extorders$D5)
scaleD6=scale(extorders$D6)
scaleD7=scale(extorders$D7)
scaleP=scale(extorders$Prop_stable)

Sscale=attributes(scaleS)$'scaled:scale'
Cscale=attributes(scaleC)$'scaled:scale'
S1scale=attributes(scaleS1)$'scaled:scale'
S2scale=attributes(scaleS2)$'scaled:scale'
S3scale=attributes(scaleS3)$'scaled:scale'
S4scale=attributes(scaleS4)$'scaled:scale'
S5scale=attributes(scaleS5)$'scaled:scale'
D1scale=attributes(scaleD1)$'scaled:scale'
D2scale=attributes(scaleD2)$'scaled:scale'
D3scale=attributes(scaleD3)$'scaled:scale'
D4scale=attributes(scaleD4)$'scaled:scale'
D5scale=attributes(scaleD5)$'scaled:scale'
D6scale=attributes(scaleD6)$'scaled:scale'
D7scale=attributes(scaleD7)$'scaled:scale'
Pscale=attributes(scaleP)$'scaled:scale'

write.table(summary(S1)$coef,file='tables/S1_lmer_table.tsv',sep='\t')
write.table(summary(S2_2)$coef,file='tables/S2_lmer_table.tsv',sep='\t')
write.table(summary(S4)$coef,file='tables/S4_lmer_table.tsv',sep='\t')
write.table(summary(S5)$coef,file='tables/S5_lmer_table.tsv',sep='\t')
write.table(summary(S3_3)$coef,file='tables/S3_lmer_table.tsv',sep='\t')
write.table(summary(D1)$coef,file='tables/D1_lmer_table.tsv',sep='\t')
write.table(summary(D2)$coef,file='tables/D2_lmer_table.tsv',sep='\t')
write.table(summary(D3)$coef,file='tables/D3_lmer_table.tsv',sep='\t')
write.table(summary(D4_2)$coef,file='tables/D4_lmer_table.tsv',sep='\t')
write.table(summary(D5)$coef,file='tables/D5_lmer_table.tsv',sep='\t')
write.table(summary(D6)$coef,file='tables/D6_lmer_table.tsv',sep='\t')
write.table(summary(D7)$coef,file='tables/D7_lmer_table.tsv',sep='\t')
write.table(summary(propmod_2)$coef,file='tables/Proportion_stable_lmer_table.tsv',sep='\t')

scaletab=rbind(c("S",attributes(scaleS)$'scaled:center',attributes(scaleS)$'scaled:scale'),
	c("C",attributes(scaleC)$'scaled:center',attributes(scaleC)$'scaled:scale'),
	c("S1",attributes(scaleS1)$'scaled:center',attributes(scaleS1)$'scaled:scale'),
	c("S2",attributes(scaleS2)$'scaled:center',attributes(scaleS2)$'scaled:scale'),
	c("S3",attributes(scaleS3)$'scaled:center',attributes(scaleS3)$'scaled:scale'),
	c("S4",attributes(scaleS4)$'scaled:center',attributes(scaleS4)$'scaled:scale'),
	c("S5",attributes(scaleS5)$'scaled:center',attributes(scaleS5)$'scaled:scale'),
	c("D1",attributes(scaleD1)$'scaled:center',attributes(scaleD1)$'scaled:scale'),
	c("D2",attributes(scaleD2)$'scaled:center',attributes(scaleD2)$'scaled:scale'),
	c("D3",attributes(scaleD3)$'scaled:center',attributes(scaleD3)$'scaled:scale'),
	c("D4",attributes(scaleD4)$'scaled:center',attributes(scaleD4)$'scaled:scale'),
	c("D5",attributes(scaleD5)$'scaled:center',attributes(scaleD5)$'scaled:scale'),
	c("D6",attributes(scaleD6)$'scaled:center',attributes(scaleD6)$'scaled:scale'),
	c("D7",attributes(scaleD7)$'scaled:center',attributes(scaleD7)$'scaled:scale'),
	c("prop",attributes(scaleP)$'scaled:center',attributes(scaleP)$'scaled:scale')
	)
colnames(scaletab)=c("Val","Center","Scale")
write.table(scaletab,file='tables/scales_for_motif_lmer.tsv',sep='\t')

