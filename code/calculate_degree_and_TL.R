library(igraph)
library(rjson)
# Convert JSON object to more R/python-readible files
# Want post-equilibrium, pre-removal structure and biomasses.

badnets=matrix(nrow=0,ncol=3)
colnames(badnets)=c("S","C","i")
for(s in seq(50,100,10)){
	print(paste0(as.character(s),' species'))
	# dir.create(paste0('../data/TLs/',as.character(s)))
	# dir.create(paste0('../data/degrees/',as.character(s)))
	for(c in seq(0.02,0.2,0.02)){
		print(paste0('connectance=',as.character(c)))
		# dir.create(paste0('../data/TLs/',as.character(s),'/',as.character(c)))
		# dir.create(paste0('../data/degrees/',as.character(s),'/',as.character(c)))
		# For some reason R can't make this work with list.files.
		for(i in seq(0,99)){
			# Select the json file to work with
			infile=paste0('../data/networks/pre_disturbance/',as.character(s),'/',as.character(c),'/initial_net_',as.character(i),'.json',sep='')
			json_data=fromJSON(file=infile)
			A=as.data.frame(json_data)
			# Preds in rows, prey in columns.
			colnames(A)=c(paste0('sp',seq(1,ncol(A))))
			rownames(A)=colnames(A)
			A <- as.matrix(A)
			# Do not count cannibalism for degree or TL
			# There's the occasional species with a cannibalistic link but no other prey
			# This is actually kind of an analogue for "bacteria" nodes
			diag(A)<-0

			# convert to edgelist and save
			edges=get.edgelist(graph.adjacency(as.matrix(A)))

			# Prey-averaged trophic level not working well since there are redundant species
			# Just going to use shortest trophic level instead - no reason why not.
			# Calculate paths between species 
			g=graph_from_edgelist(edges,directed=TRUE)
			degrees=degree(g,mode="total",loops=TRUE)
			write.table(degrees,file=paste0('../data/degrees/',as.character(s),'/',as.character(c),'/degrees_',as.character(i),'.csv',sep=''),sep=',')
			basals=rownames(A[rowSums(A)==0,])
			paths=distances(g,v=basals,mode='in')
			STL <- 1 + as.numeric(apply(paths,2,min))
			if(length(which(STL==Inf)>0)){
				badnets=rbind(badnets,c(s,c,i))
			}
			write.table(STL,file=paste0('../data/TLs/',as.character(s),'/',as.character(c),'/SWTL_',as.character(i),'.csv',sep=''),sep=',')

		}
	}
}

write.table(badnets,file='multicomponent_networks_to_rerun.tsv',sep='\t')



