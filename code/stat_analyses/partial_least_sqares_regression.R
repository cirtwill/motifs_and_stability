library(vegan)
library(pls)

# Normalizing based on network prevalence of motifs (Z-scores)

allC=matrix(nrow=0,ncol=13)
allZ=matrix(nrow=0,ncol=13)
metadata=matrix(nrow=0,ncol=5)
for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)
		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'_motifs.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE) # Species IDs, then extinction order for each removal, then roles
			persist=data[,2:(S+1)]
			mets=cbind(rowMeans(persist),rep(S,S),rep(C,S),data[,(S+2):(S+3)])
			colnames(mets)[1:3]=c("Persistence","S","C")
			for(r in 1:nrow(mets)){
				rownames(mets)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',data$Species[r])
			}
			counts=data[,(S+4):(S+16)]
			for(r in 1:nrow(counts)){
				rownames(counts)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',data$Species[r])
			}
			Zs=data[,(S+17):(S+29)]
			for(r in 1:nrow(Zs)){
				rownames(Zs)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',data$Species[r])
			}
			metadata=rbind(metadata,mets)
			allC=rbind(allC,counts)
			allZ=rbind(allZ,Zs)
		}
	}
}

# Centering and scaling all variables
# Raw roles
	# Scaling raw motif roles also
	raw=plsr(metadata$Persistence~metadata$S*metadata$C+metadata$STL+metadata$Degree+as.matrix(allC),validation="CV",scale=TRUE,center=TRUE)
	write.table(as.data.frame(RMSEP(raw)$val),sep='\t',file='../../data/PLS_regression/raw_errors.tsv')
	ncomp_raw=selectNcomp(raw,method="onesigma",plot=TRUE)
	raw_opt=plsr(metadata$Persistence~metadata$S*metadata$C+metadata$STL+metadata$Degree+as.matrix(allC),validation="CV",,scale=TRUE,center=TRUE,ncomp=5)
	write.table(as.data.frame(raw_opt$coefficients),file='../../data/PLS_regression/raw_coefficients.tsv',sep='\t')


# Network-normalized roles
	netnorm=plsr(metadata$Persistence~metadata$S*metadata$C+metadata$STL+metadata$Degree+as.matrix(allZ),validation="CV",scale=TRUE,center=TRUE)
	write.table(as.data.frame(RMSEP(netnorm)$val),sep='\t',file='../../data/PLS_regression/netnorm_errors.tsv')
	ncomp_netnorm=selectNcomp(netnorm,method="onesigma",plot=TRUE)
	netnorm_opt=plsr(metadata$Persistence~metadata$S*metadata$C+metadata$STL+metadata$Degree+as.matrix(allZ),validation="CV",scale=TRUE,center=TRUE,ncomp=5)
	write.table(as.data.frame(netnorm_opt$coefficients),file='../../data/PLS_regression/netnorm_coefficients.tsv',sep='\t')


# Degree-normalized roles
	allC_norm=allC/rowSums(allC)
	degnorm=plsr(metadata$Persistence~metadata$S*metadata$C+metadata$STL+metadata$Degree+as.matrix(allC_norm),validation="CV",scale=TRUE,center=TRUE)
	write.table(as.data.frame(RMSEP(degnorm)$val),sep='\t',file='../../data/PLS_regression/degnorm_errors.tsv')
	ncomp_degnorm=selectNcomp(degnorm,method="onesigma",plot=TRUE)
	degnorm_opt=plsr(metadata$Persistence~metadata$S*metadata$C+metadata$STL+metadata$Degree+as.matrix(allC_norm),validation="CV",scale=TRUE,center=TRUE,ncomp=5)
	write.table(as.data.frame(degnorm_opt$coefficients),file='../../data/PLS_regression/degnorm_coefficients.tsv',sep='\t')

scales=c(attributes(scale(metadata))$'scaled:scale',attributes(scale(allC))$'scaled:scale',attributes(scale(allZ))$'scaled:scale')
write.table(scales,file='../../data/PLS_regression/scales.tsv',sep='\t')

# # Test datasets
# TallC=matrix(nrow=0,ncol=13)
# TallZ=matrix(nrow=0,ncol=13)
# Tmetadata=matrix(nrow=0,ncol=5)
# for(S in seq(50,100,10)){
# 	print(S)
# 	for(C in seq(0.02,0.2,0.02)){
# 		print(C)
# 		for(i in seq(50,99)){
# 			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'_motifs.tsv',sep='')
# 			data=read.table(datafile,sep='\t',header=TRUE) # Species IDs, then extinction order for each removal, then roles
# 			persist=data[,2:(S+1)]
# 			mets=cbind(rowMeans(persist),rep(S,S),rep(C,S),data[,(S+2):(S+3)])
# 			colnames(mets)[1:3]=c("Persistence","S","C")
# 			for(r in 1:nrow(mets)){
# 				rownames(mets)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',data$Species[r])
# 			}
# 			counts=data[,(S+4):(S+16)]
# 			for(r in 1:nrow(counts)){
# 				rownames(counts)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',data$Species[r])
# 			}
# 			Zs=data[,(S+17):(S+29)]
# 			for(r in 1:nrow(Zs)){
# 				rownames(Zs)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',data$Species[r])
# 			}
# 			Tmetadata=rbind(Tmetadata,mets)
# 			TallC=rbind(TallC,counts)
# 			TallZ=rbind(TallZ,Zs)
# 		}
# 	}
# }

# # Raw motif predictions
# 	Tmets=scale(Tmetadata)
# 	newdat=cbind(Tmets,scale(TallC))
# 	raw_predictions=predict(raw_opt,newdata=newdat)
# 	predout=cbind(Tmets[,1],raw_predictions)
# 	colnames(predout)=c("Obs","Pred")
# 	write.table(predout,file='../../data/PLS_regression/raw_predictions.tsv',sep='\t')

# # Degree-normalized predictions
# 	TallC_norm=scale(TallC/rowSums(TallC))

# 	newdat2=cbind(Tmets,TallC_norm)
# 	degnorm_predictions=predict(degnorm_opt,newdata=newdat2)
# 	predout2=cbind(Tmets[,1],degnorm_predictions)
# 	colnames(predout2)=c("Obs","Pred")
# 	write.table(predout2,file='../../data/PLS_regression/degnorm_predictions.tsv',sep='\t')

# # Network normalized predictions
# 	Tmets=scale(Tmetadata)
# 	newdat=cbind(Tmets,scale(TallZ))
# 	netnorm_predictions=predict(netnorm_opt,newdata=newdat)
# 	predout3=cbind(Tmets[,1],netnorm_predictions)
# 	colnames(predout3)=c("Obs","Pred")
# 	write.table(predout3,file='../../data/PLS_regression/netnorm_predictions.tsv',sep='\t')
