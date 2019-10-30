library(vegan)
library(lme4)
library(lmerTest)

# Would be good to get the predicted extinction orders and observed, 

# Test mean extinction order vs. dissimilarit to the removed species.
for(S in seq(50,100,10)){
	meta_extorders=matrix(nrow=0,ncol=8)
	colnames(meta_extorders)=c("S","C","Network","Dissim","Extorder","Path","Removed_sp","Focal_sp")
	print(S)
	for(C in seq(0.02,0.2,0.02)){

		extorders=matrix(nrow=0,ncol=6)
		colnames(extorders)=c("Network","Dissim","Extorder","Path","Removed_sp","Focal_sp")

		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)
			pathfile=paste0('../../data/paths/pre_disturbance/',as.character(S),'/',as.character(C),'/initial_paths_',as.character(i),'.tsv',sep='')
			paths=read.table(pathfile,sep='\t',header=TRUE)

			data=data[order(data$Species),]
			paths=paths[order(rownames(paths)),]

			roles=data[,(S+2):ncol(data)] # Species IDs, then extinction order for each removal, then roles
			roledist=vegdist(roles,method="bray")

			for(j in seq(1,S)){
				dists=as.matrix(roledist)[,j]
				exts=data[,1+j]
				colname=paste0("sp",j,sep='')
				minipaths=paths[,colname]				
				smalldat=cbind(rep(i,S),dists,exts,minipaths,rep(j,S),data[,1])
				extorders=rbind(extorders,smalldat)
			}
		}
		mets=cbind(rep(S,nrow(extorders)),rep(C,nrow(extorders)),extorders)
		meta_extorders=rbind(meta_extorders,mets)

		extorders=as.data.frame(extorders[which(!extorders[,5]==extorders[,6]),])

		mod=lmer(Extorder~scale(Dissim)*scale(Path)+(1|Network),data=extorders)
		preds=predict(mod)
		extorders=cbind(extorders,preds)
		colnames(extorders)[7]="Predicted"

		write.table(summary(mod)$coefficients,file=paste0('../../data/summaries/rolesim_extorder/',as.character(S),'/role_exttime_lm_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
		write.table(extorders,file=paste0('../../data/summaries/rolesim_extorder/',as.character(S),'/role_exttime_preds_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
	}
	meta_extorders=as.data.frame(meta_extorders[which(!meta_extorders[,5]==meta_extorders[,6]),])
	meta_mod=lmer(Extorder~scale(Dissim)*scale(Path)*scale(C)+(1|Network),data=meta_extorders)
	write.table(summary(meta_mod)$coefficients,file=paste0('../../data/summaries/rolesim_extorder/',as.character(S),'/role_exttime_lm_with_C.tsv'),sep='\t')
	scales=rbind(c(attributes(scale(meta_extorders$Dissim))$'scaled:center',attributes(scale(meta_extorders$Path))$'scaled:center',attributes(scale(meta_extorders$C))$'scaled:center'),
				c(attributes(scale(meta_extorders$Dissim))$'scaled:scale',attributes(scale(meta_extorders$Path))$'scaled:scale',attributes(scale(meta_extorders$C))$'scaled:scale'))
	scales=as.data.frame(scales)
	rownames(scales)=c("center","scale")
	colnames(scales)=c("Dissim","Path","C")
	write.table(scales,file=paste0('../../data/summaries/rolesim_extorder/scales_',as.character(S),'.tsv'),sep='\t')
}


# Test whether first extinction more similar to removed spp. than others

