library(vegan)

# PCA of roles across all webs in an S:C combination
for(S in seq(50,100,10)){
	dir.create(paste0('../../data/summaries/role_PCA/',as.character(S)))
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)
		allroles=matrix(nrow=0,ncol=30)

		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE) # Species IDs, then extinction order for each removal, then roles
			roles=data[,(S+2):ncol(data)]
			for(r in 1:nrow(roles)){
				rownames(roles)[r]=paste0(as.character(i),data$Species[r])
			}
			allroles=rbind(allroles,roles)
		}
		PCA=prcomp(allroles)

		write.table(PCA$rotation[,1:3],file=paste0('../../data/summaries/role_PCA/',as.character(S),'/role_PCA_axes_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
		write.table(PCA$sdev,file=paste0('../../data/summaries/role_PCA/',as.character(S),'/PCA_axis_var_explained_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
		write.table(PCA$x[,1:3],file=paste0('../../data/summaries/role_PCA/',as.character(S),'/role_PCA_positions_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')		
	}
}

# Mega PCA ocross all webs, with S and C.
allroles=matrix(nrow=0,ncol=30)
for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)

		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE) # Species IDs, then extinction order for each removal, then roles
			roles=cbind(rep(S,S),rep(C,S),data[,(S+2):ncol(data)])
			colnames(roles)[1:2]=c("S","C")
			for(r in 1:nrow(roles)){
				rownames(roles)[r]=paste0(as.character(S),'_',as.character(C),'_',as.character(i),'_',data$Species[r])
			}
			allroles=rbind(allroles,roles)
		}
	}
}

PCA=prcomp(allroles)

write.table(PCA$rotation[,1:3],file=paste0('../../data/summaries/role_PCA/role_PCA_axes_global.tsv',sep=''),sep='\t')
write.table(PCA$sdev,file=paste0('../../data/summaries/role_PCA/PCA_axis_var_explained_global.tsv',sep=''),sep='\t')
write.table(PCA$x[,1:3],file=paste0('../../data/summaries/role_PCA/role_PCA_positions_global.tsv',sep=''),sep='\t')		
