library(vegan)

# Not sure how this is going, really.
# Testing: Does role similarity affect extinction order?

# Test mean extinction order vs. roles
for(S in seq(50,100,10)){
	dir.create(paste0('../../data/summaries/rolesim_extorder/',as.character(S)))
	print(S)
	for(C in seq(0.02,0.2,0.02)){
		if(C==0.1){
			print(C)
		}
		extorders=matrix(nrow=0,ncol=3)
		colnames(extorders)=c("Network","Mantel","p")

		for(i in seq(0,99)){
			datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
			data=read.table(datafile,sep='\t',header=TRUE)

			roles=data[,(S+2):ncol(data)] # Species IDs, then extinction order for each removal, then roles
			roledist=vegdist(roles,method="bray")

			exts=data[,2:(S+1)]
			rownames(exts)=data$Species
			extdist=vegdist(exts,method="euclid")

			Mtest=mantel(roledist,exts,permutations=9999)
			extorders=rbind(extorders,c(i,Mtest$statistic,Mtest$signif))
		}
		write.table(extorders,file=paste0('../../data/summaries/rolesim_extorder/',as.character(S),'/role_exttime_Mantel_',as.character(S),'_',as.character(C),'_',as.character(i),'.tsv',sep=''),sep='\t')
	}
}


# Test whether first extinction more similar to removed spp. than others

