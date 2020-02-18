library(lmerTest)

all_positions=matrix(nrow=0,ncol=3)
all_extorders=matrix(nrow=0,ncol=5)

for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		positions=read.table(file=paste0('../../data/summaries/role_PCA/',as.character(S),'/role_PCA_positions_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t',header=TRUE)		
		all_positions=rbind(all_positions,positions)

		print(C)
		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)
			extorders=data[,2:(S+1)]
			meanorders=rowMeans(extorders)

			orders=as.data.frame(cbind(rep(S,length(meanorders)),rep(C,length(meanorders)),rep(i,length(meanorders)),as.character(data$Species),meanorders))
			colnames(orders)=c("S","C","Network","Species","Extorder")
			for(r in 1:nrow(orders)){
				rownames(orders)[r]=paste0(as.character(i),orders[,4][r])
			}
			all_extorders=rbind(all_extorders,orders)
		}
	}
}

if(length(which(rownames(all_extorders)==rownames(all_positions)))==nrow(all_extorders)){
	alldata=cbind(all_extorders,all_positions)
}
alldata$rando=paste0(as.character(alldata$S),as.character(alldata$C),as.character(alldata$Network))
alldata=as.data.frame(alldata)
alldata$S=as.numeric(as.character(alldata$S))
alldata$C=as.numeric(as.character(alldata$C))
alldata$Extorder=as.numeric(as.character(alldata$Extorder))

PC1=lmer(Extorder~scale(PC1)*scale(S)*scale(C)+(1|rando),data=alldata,family='poisson')
PC2=lmer(Extorder~scale(PC2)*scale(S)*scale(C)+(1|rando),data=alldata,family='poisson')
PC3_full=lmer(Extorder~scale(PC3)*scale(S)*scale(C)+(1|rando),data=alldata,family='poisson')
PC3_2=lmer(Extorder~scale(PC3)*scale(S)+scale(S)*scale(C)+scale(PC3)*scale(C)+(1|rando),data=alldata,family='poisson')
PC3_3=lmer(Extorder~scale(PC3)*scale(S)+scale(S)*scale(C)+(1|rando),data=alldata,family='poisson')


scaleS=scale(alldata$S)
scaleC=scale(alldata$C)
scalePC1=scale(alldata$PC1)
scalePC2=scale(alldata$PC2)
scalePC3=scale(alldata$PC3)

Sscale=attributes(scaleS)$'scaled:scale'
Cscale=attributes(scaleC)$'scaled:scale'
scale1=attributes(scalePC1)$'scaled:scale'
scale2=attributes(scalePC2)$'scaled:scale'
scale3=attributes(scalePC3)$'scaled:scale'

write.table(summary(PC1)$coef,file='tables/PC1_lmer_table.tsv',sep='\t')
write.table(summary(PC2)$coef,file='tables/PC2_lmer_table.tsv',sep='\t')
write.table(summary(PC3_3)$coef,file='tables/PC3_lmer_table.tsv',sep='\t')

scaletab=rbind(c("S",attributes(scaleS)$'scaled:center',attributes(scaleS)$'scaled:scale'),
	c("C",attributes(scaleC)$'scaled:center',attributes(scaleC)$'scaled:scale'),
	c("PC1",attributes(scalePC1)$'scaled:center',attributes(scalePC1)$'scaled:scale'),
	c("PC2",attributes(scalePC2)$'scaled:center',attributes(scalePC2)$'scaled:scale'),
	c("PC3",attributes(scalePC3)$'scaled:center',attributes(scalePC3)$'scaled:scale')
	)

colnames(scaletab)=c("Val","Center","Scale")
write.table(scaletab,file='tables/scales_for_motif_lmer.tsv',sep='\t')