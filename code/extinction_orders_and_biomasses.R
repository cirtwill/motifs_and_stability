library(igraph)
library(rjson)
# Convert JSON object to more R/python-readible files
# Want post-equilibrium, pre-removal structure and biomasses.

badnets=matrix(nrow=0,ncol=3)
colnames(badnets)=c("S","C","i")
for(s in seq(50,100,10)){
	print(paste0(as.character(s),' species'))
	dir.create(paste0('../data/networks/extorder/',as.character(s)))
	for(c in seq(0.02,0.2,0.02)){
		print(paste0('connectance=',as.character(c)))
		dir.create(paste0('../data/networks/extorder/',as.character(s),'/',as.character(c)))

		# For some reason R can't make this work with list.files.
		for(i in seq(0,99)){
			meta_extorder=matrix(nrow=s,ncol=0)
			meta_biomasses=matrix(nrow=0,ncol=s+1)
			for(rem in seq(1,s)){
				# Select the json file to work with
				infile=paste0('../data/networks/removals/',as.character(s),'/',as.character(c),'/network_',as.character(i),':removal_',as.character(rem),'.json',sep='')

				# Open the file
				# Is a list of biomasses after each round of simulations
				# [[i]] gives the biomasses for species i
				ext_time=matrix(nrow=s,ncol=2)
				ext_biomes=matrix(nrow=51,ncol=s+1)
				ext_biomes[,1]=c(rep(rem,51))
				json_data=fromJSON(file=infile)
				for(j in 1:length(json_data)){
					zeros=which(json_data[[j]][1:51]==0)
					ext_time[j,]=c(j,50-length(zeros))
					ext_biomes[,j+1]=json_data[[j]]
				}
				colnames(ext_time)=c(paste0("Species_",as.character(rem)),paste0("Ext_step_",as.character(rem)))
				ext_time=ext_time[order(ext_time[,2]),]
				meta_extorder=cbind(meta_extorder,ext_time)
				meta_biomasses=rbind(meta_biomasses,ext_biomes)
			}
			write.table(t(meta_extorder),file=paste0('../data/networks/extorder/',as.character(s),'/',as.character(c),'/network_',as.character(i),'_extinction_orders.tsv',sep=''),sep='\t',row.names=TRUE,col.names=FALSE)
			meta_biomasses=as.data.frame(meta_biomasses)
			colnames(meta_biomasses)=c("Removal","Biomasses")
			write.table(meta_biomasses,file=paste0('../data/networks/extorder/',as.character(s),'/',as.character(c),'/network_',as.character(i),'_biomasses.tsv',sep=''),sep='\t',row.names=TRUE,col.names=FALSE)
		}
	}
}

