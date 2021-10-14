library(vegan)

# Testing: Do species with similar mean extinction orders have similar roles?
# For each size-connectance combination, running a 9999-perm PERMANOVA stratified by network
# Also calculating an ANOVA of dispersion vs. quantile of extinction
# Also printing dispersion values vs. quantile for later plotting.
# Will need to be made human-readable with Python

for(method in c("count","freq","Z")){
	for(S in seq(50,100,10)){
	# for(S in seq(50,100,10)){
		dir.create(paste0('../../data/permanova/',method,'/perms/',as.character(S)))
		dir.create(paste0('../../data/permanova/',method,'/disps/',as.character(S)))
		print(S)
		for(C in seq(0.02,0.2,0.02)){
			dir.create(paste0('../../data/permanova/',method,'/perms/',as.character(S),'/',as.character(C)))
			dir.create(paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C)))
			print(C)
			extorders=matrix(nrow=0,ncol=5)
			colnames(extorders)=c("Network","Extinction_correlation","Mean_order","Rank","Quantile")
			allroles=matrix(nrow=0,ncol=13) # Count, freq, or Z depending on method.

			for(i in seq(0,99)){
				datafile=paste0('../../data/roles/matched_to_extorder/',as.character(S),'/',as.character(C),'/network_',as.character(i),'.tsv',sep='')
				data=read.table(datafile,sep='\t',header=TRUE)

				roles=data[,(S+4):(S+16)]
				Zs=data[,(S+17):ncol(data)] 
				if(method=='count'){
					allroles=rbind(allroles,roles)
				} else {
					if(method=='freq'){
						allroles=rbind(allroles,roles/rowSums(roles))
					} else {
						allroles=rbind(allroles,Zs)
					}
				}

				# Mean extinction order vs. roles, should also store strength of correlation.
				extcorr=mean(cor(data[,2:(S+1)]))
				mean_ext=rowSums(data[,2:(S+1)])/S # Higher rank = longer mean time to extinction.
				ranked=rank(mean_ext,ties.method=c("average"))/S # Rank as a fraction of the maximum rank, can be rounded.
				quant=round(ranked,1) # grouping in blocks of 10% - blocks may not contain equal numbers of species since there were tied ranks
				extdat=cbind(rep(i,S),rep(extcorr,S),mean_ext,ranked,quant)

				colnames(extdat)=c("Network","Extinction_correlation","Mean_order","Rank","Quantile")
				extorders=rbind(extorders,extdat)

			}

			extorders=as.data.frame(extorders)
			roledist=vegdist(allroles,method="bray")
			# For speed, especially in the biggest webs, doing fewer permutation at a time. Will collect them and get p-values later.
			if(S<100){
				for(j in seq(1,10)){
					meanperm=adonis(roledist~extorders$Mean_order,strata=as.factor(extorders$Network),permutations=999)			
					summa=meanperm$aov.tab[1,]
					results=c(S,C,mean(extorders$Extinction_correlation),summa,meanperm$f.perms[,1])
					names(results)=c("S","C","ext_corr","DF","SS","MS","F.model","R2","pval","perms")
					write.table(results,file=paste0('../../data/permanova/',method,'/perms/',as.character(S),'/',as.character(C),'/mean_extorder_vs_roles_',as.character(S),'_',as.character(C),'_',as.character(j),'.tsv',sep=''),sep='\t')
				}			
				meanperm=adonis(roledist~extorders$Mean_order,strata=as.factor(extorders$Network),permutations=9)			
				summa=meanperm$aov.tab[1,]
				results=c(S,C,mean(extorders$Extinction_correlation),summa,meanperm$f.perms[,1])
				names(results)=c("S","C","ext_corr","DF","SS","MS","F.model","R2","pval","perms")
				write.table(results,file=paste0('../../data/permanova/',method,'/perms/',as.character(S),'/',as.character(C),'/mean_extorder_vs_roles_',as.character(S),'_',as.character(C),'_x.tsv',sep=''),sep='\t')
				disper=betadisper(roledist,extorders$Quantile,type="centroid")
				dists=cbind(disper$distances,disper$group)
				colnames(dists)=c("Dispersion","Quantile")
				anny=anova(disper)
				write.table(dists,file=paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C),'/dispersion_vs_quantile_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
				write.table(anny,file=paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C),'/dispersion_anova_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
			} else {
				for(j in seq(1,999)){
					meanperm=adonis(roledist~extorders$Mean_order,strata=as.factor(extorders$Network),permutations=10)			
					summa=meanperm$aov.tab[1,]
					results=c(S,C,mean(extorders$Extinction_correlation),summa,meanperm$f.perms)
					names(results)=c("S","C","ext_corr","DF","SS","MS","F.model","R2","pval","perm")
					write.table(results,file=paste0('../../data/permanova/',method,'/perms/',as.character(S),'/',as.character(C),'/mean_extorder_vs_roles_',as.character(S),'_',as.character(C),'_',as.character(j),'.tsv',sep=''),sep='\t')
				}						
				meanperm=adonis(roledist~extorders$Mean_order,strata=as.factor(extorders$Network),permutations=9)			
				summa=meanperm$aov.tab[1,]
				results=c(S,C,mean(extorders$Extinction_correlation),summa,meanperm$f.perms[,1])
				names(results)=c("S","C","ext_corr","DF","SS","MS","F.model","R2","pval","perms")
				write.table(results,file=paste0('../../data/permanova/',method,'/perms/',as.character(S),'/',as.character(C),'/mean_extorder_vs_roles_',as.character(S),'_',as.character(C),'_x.tsv',sep=''),sep='\t')		
				disper=betadisper(roledist,extorders$Quantile,type="centroid")
				dists=cbind(disper$distances,disper$group)
				colnames(dists)=c("Dispersion","Quantile")
				anny=anova(disper)
				write.table(dists,file=paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C),'/dispersion_vs_quantile_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
				write.table(anny,file=paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C),'/dispersion_anova_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
			}
		}
	}
}
# Test whether first extinction more similar to removed spp. than others
