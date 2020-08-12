library(lmerTest)
library(vegan)
library(MuMIn)

# The glmer's are EXTREMELY slow... data may not be distributed right?
# Or group species based on extinction order, 100>1.

all_extorders=matrix(nrow=0,ncol=5)
all_positions=read.table(file=paste0('../../data/summaries/role_PCA/role_PCA_positions_global.tsv',sep=''),sep='\t',header=TRUE)		
bins=matrix(nrow=0,ncol=6)
for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		# Using positions based on PCA of webs within an S:C block
		# all_positions=rbind(all_positions,positions)
		exttime=matrix(nrow=0,ncol=5)
		print(C)
		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)
			extorders=data[,2:(S+1)]
			meanorders=rowMeans(extorders)

			orders=as.data.frame(cbind(rep(S,length(meanorders)),rep(C,length(meanorders)),rep(i,length(meanorders)),as.character(data$Species),meanorders))
			colnames(orders)=c("S","C","Network","Species","Extorder")
			for(r in 1:nrow(orders)){
				rownames(orders)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',orders$Species[r])
			}
			all_extorders=rbind(all_extorders,orders)
			exttime=rbind(exttime,orders)
		}
		positions=all_positions[rownames(exttime),]
		if(length(which(rownames(exttime)==rownames(positions)))==nrow(exttime)){
			alldata=cbind(exttime,positions)
		}
		# alldata$rando=paste0(as.character(alldata$S),as.character(alldata$C),as.character(alldata$Network))
		alldata=as.data.frame(alldata)
		alldata$S=as.numeric(as.character(alldata$S))
		alldata$C=as.numeric(as.character(alldata$C))
		alldata$Extorder=as.numeric(as.character(alldata$Extorder))

		# Group species based on extorder, take means
		alldata=alldata[order(alldata$Extorder),]
		r=1
		for(i in 1:(nrow(alldata)/100)){
			rmax=r+99
			subset=alldata[r:rmax,]
			binres=c(colMeans(subset[,1:2]),colMeans(subset[,5:8]))
			bins=rbind(bins,binres)
			r=r+100
		}
		# That is fast :)

		}
	}
bins=as.data.frame(bins)
colnames(bins)=c('S','C','Extorder','PC1','PC2','PC3')
# Works quickly for all C within 50 S
# These glms work, so now I need to figure out what I want to do. Dredge isn't working super well (can't calculate Loglik).

# Extinction orders non-integer, Poisson really doesn't work. Back to LM.
PCs=lm(Extorder~scale(PC1)+scale(PC2)+scale(PC3)+scale(S)+scale(C)+scale(S):scale(C),data=bins,na.action='na.fail')
# S:C not significant
# write.table(summary(PCs)$coef,file=paste0('tables/PC_lmer_table.tsv'),sep='\t')

all_extorders$S=as.numeric(as.character(all_extorders$S))
all_extorders$C=as.numeric(as.character(all_extorders$C))
all_extorders$Extorder=as.numeric(as.character(all_extorders$Extorder))
allPCs=lm(all_extorders$Extorder~scale(PC1)+scale(PC2)+scale(PC3)+scale(all_extorders$S)+scale(all_extorders$C)+scale(all_extorders$S):scale(all_extorders$C),data=all_positions,na.action='na.fail')
write.table(summary(allPCs)$coef,file=paste0('tables/PC_lmer_table.tsv'),sep='\t')

scaleS=scale(all_extorders$S)
scaleC=scale(all_extorders$C)
scalePC1=scale(all_positions$PC1)
scalePC2=scale(all_positions$PC2)
scalePC3=scale(all_positions$PC3)

scaletab=rbind(c("S",attributes(scaleS)$'scaled:center',attributes(scaleS)$'scaled:scale'),
	c("C",attributes(scaleC)$'scaled:center',attributes(scaleC)$'scaled:scale'),
	c("PC1",attributes(scalePC1)$'scaled:center',attributes(scalePC1)$'scaled:scale'),
	c("PC2",attributes(scalePC2)$'scaled:center',attributes(scalePC2)$'scaled:scale'),
	c("PC3",attributes(scalePC3)$'scaled:center',attributes(scalePC3)$'scaled:scale')
	)

colnames(scaletab)=c("Val","Center","Scale")
write.table(scaletab,file='tables/scales_for_motif_lmer.tsv',sep='\t')