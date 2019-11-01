library(vegan)

# I think we want to do one permanova across all 100 networks for a size/connectance
# Testing: Do species with similar mean extinction orders have similar roles?
# 90 and 100-species ones don't run here, need to be on the mac

for(S in seq(50,100,10)){
	dir.create(paste0('../../data/summaries/role_PCA/',as.character(S)))
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		print(C)
		allroles=matrix(nrow=0,ncol=30)

		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)

			roles=data[,(S+2):ncol(data)] # Species IDs, then extinction order for each removal, then roles
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

# Test whether first extinction more similar to removed spp. than others
