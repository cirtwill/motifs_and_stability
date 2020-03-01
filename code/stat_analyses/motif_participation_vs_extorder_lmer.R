library(lmerTest)
library(MuMIn)


# Subset species like with PCA axes. Then fit one big model, simplify.




# Testing: Do species with similar mean extinction orders have similar roles?
# No convergence??
bins=matrix(nrow=0,ncol=29)
# extorders=matrix(nrow=0,ncol=5)
# colnames(extorders)=c("S","C","Network","Extinction_correlation","Mean_order")
# allroles=matrix(nrow=0,ncol=30)
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
		}
		alldata=cbind(exttime,allroles)
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

# extorders=as.data.frame(cbind(extorders,allroles))
# extorders$rando=paste0(as.character(extorders$S),as.character(extorders$C),as.character(extorders$Network))
# # Most interested in motifs S1, S2, S4, S5

# # Still not sure about how the random effects should work

# # Stable motifs
# S1=lmer(Mean_order~zS1+(1|rando),data=extorders,family='poisson')
# S2=lmer(Mean_order~zS2+(1|rando),data=extorders,family='poisson')
# S4=lmer(Mean_order~zS4+(1|rando),data=extorders,family='poisson')
# S5=lmer(Mean_order~zS5+(1|rando),data=extorders,family='poisson')

# # Unstable motifs
# S3=lmer(Mean_order~zS3+(1|rando),data=extorders,family='poisson')
# # Effect of loop should depend on web conditions
# D1=lmer(Mean_order~zD1+(1|rando),data=extorders,
# 	control=lmerControl(optimizer="bobyqa"),family='poisson')
# # Failure to converge.

# D2=lmer(Mean_order~zD2+(1|rando),data=extorders,family='poisson')
# D3=lmer(Mean_order~zD3+(1|rando),data=extorders,family='poisson')
# D4=lmer(Mean_order~zD4+(1|rando),data=extorders,family='poisson')
# D5=lmer(Mean_order~zD5+(1|rando),data=extorders,family='poisson')
# D6=lmer(Mean_order~zD6+(1|rando),data=extorders,family='poisson')
# D7=lmer(Mean_order~zD7+(1|rando),data=extorders,family='poisson')
# # Effects of all two-ways should depend on web conditions
# # Maybe I should add random effects somehow? Or plot the actual data?

# deg=lmer(Mean_order~Degree+(1|rando),data=extorders,family='poisson')
# TL=lmer(Mean_order~STL+(1|rando),data=extorders,family='poisson')
# # TL there's still quite a few infinite trophic levels kicking around.  

# write.table(summary(S1)$coef,file='tables/S1_lmer_table.tsv',sep='\t')
# write.table(summary(S2)$coef,file='tables/S2_lmer_table.tsv',sep='\t')
# write.table(summary(S4)$coef,file='tables/S4_lmer_table.tsv',sep='\t')
# write.table(summary(S5)$coef,file='tables/S5_lmer_table.tsv',sep='\t')
# write.table(summary(S3)$coef,file='tables/S3_lmer_table.tsv',sep='\t')
# write.table(summary(D1)$coef,file='tables/D1_lmer_table.tsv',sep='\t')
# write.table(summary(D2)$coef,file='tables/D2_lmer_table.tsv',sep='\t')
# write.table(summary(D3)$coef,file='tables/D3_lmer_table.tsv',sep='\t')
# write.table(summary(D4)$coef,file='tables/D4_lmer_table.tsv',sep='\t')
# write.table(summary(D5)$coef,file='tables/D5_lmer_table.tsv',sep='\t')
# write.table(summary(D6)$coef,file='tables/D6_lmer_table.tsv',sep='\t')
# write.table(summary(D7)$coef,file='tables/D7_lmer_table.tsv',sep='\t')


biglm=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zS4)+scale(zS5)+scale(zD1)+scale(zD2)+scale(zD3)+scale(zD4)+scale(zD5)+scale(zD6)+scale(zD7)+scale(S)+scale(C),data=bins,na.action='na.fail')
lm2=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zS4)+scale(zS5)+scale(zD1)+scale(zD2)+scale(zD3)+scale(zD5)+scale(zD6)+scale(zD7)+scale(S)+scale(C),data=bins)
lm3=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zS4)+scale(zS5)+scale(zD2)+scale(zD3)+scale(zD5)+scale(zD6)+scale(zD7)+scale(S)+scale(C),data=bins)
lm4=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zS4)+scale(zD2)+scale(zD3)+scale(zD5)+scale(zD6)+scale(zD7)+scale(S)+scale(C),data=bins)
lm5=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zD2)+scale(zD3)+scale(zD5)+scale(zD6)+scale(zD7)+scale(S)+scale(C),data=bins)
lm6=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zD2)+scale(zD3)+scale(zD5)+scale(zD6)+scale(S)+scale(C),data=bins)
lm7=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zD2)+scale(zD3)+scale(zD6)+scale(S)+scale(C),data=bins)
lm8=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zD2)+scale(zD3)+scale(S)+scale(C),data=bins)
lm9=lm(Mean_order~scale(zS1)+scale(zS2)+scale(zS3)+scale(zD2)+scale(S)+scale(C),data=bins)
lm10=lm(Mean_order~scale(zS1)+scale(zS3)+scale(zD2)+scale(S)+scale(C),data=bins)
lm11=lm(Mean_order~scale(zS1)+scale(zS3)+scale(S)+scale(C),data=bins)
lm12=lm(Mean_order~scale(zS1)+scale(S)+scale(C),data=bins)
# S3 was significant up to lm10, but now it's not.
# Weird.
# Also not significant in the average model though, so I guess it's right.

avmod=model.avg(dredge(biglm),subset=delta<2)
# S1 is the only significant one after S and C.

