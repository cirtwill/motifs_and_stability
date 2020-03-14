library(igraph)
# Want post-equilibrium, pre-removal structure and biomasses.

badnets=matrix(nrow=0,ncol=3)
colnames(badnets)=c("S","C","i")
for(s in seq(50,100,10)){
	print(paste0(as.character(s),' species'))
	dir.create(paste0('../data/matrices/pre_disturbance/',as.character(s)))
	dir.create(paste0('../data/edgelists/pre_disturbance/',as.character(s)))
	dir.create(paste0('../data/paths/pre_disturbance/',as.character(s)))
	for(c in seq(0.02,0.2,0.02)){
		print(paste0('connectance=',as.character(c)))
		dir.create(paste0('../data/matrices/pre_disturbance/',as.character(s),'/',as.character(c)))
		dir.create(paste0('../data/edgelists/pre_disturbance/',as.character(s),'/',as.character(c)))
		dir.create(paste0('../data/paths/pre_disturbance/',as.character(s),'/',as.character(c)))
		# For some reason R can't make this work with list.files.
		for(i in seq(0,99)){
			# Select the json file to work with
			infile=paste0('../data/networks/pre_disturbance/',as.character(s),'/',as.character(c),'/initial_net_',as.character(i),'.csv',sep='')

			# Open the network file
			A=as.data.frame(read.csv(infile))
			# Preds in rows, prey in columns.
			colnames(A)=c(paste0('sp',seq(1,ncol(A))))
			rownames(A)=colnames(A)

			# convert to edgelist and save
			edges=get.edgelist(graph.adjacency(as.matrix(A)))

			# Calculate paths between species 
			g=graph_from_edgelist(edges,directed=TRUE)
			write.table(A,file=paste0('../data/matrices/pre_disturbance/',as.character(s),'/',as.character(c),'/initial_matrix_',as.character(i),'.csv',sep=''),sep=',')
			num_components=components(g)$no
			# If there is only one component, save the paths.
			if(num_components==1){
				# Obtain a matrix of shorest path lengths between all species
				write.table(edges,file=paste0('../data/edgelists/pre_disturbance/',as.character(s),'/',as.character(c),'/initial_edges_',as.character(i),'.tsv',sep=''),sep='\t',row.names=FALSE,col.names=FALSE)
				paths=distances(g)
				write.table(paths,file=paste0('../data/paths/pre_disturbance/',as.character(s),'/',as.character(c),'/initial_paths_',as.character(i),'.tsv',sep=''),sep='\t',row.names=TRUE,col.names=TRUE)
			# If there are multiple components, remove the networks and append the number for another round.
			} else {
				badnets=rbind(badnets,c(s,c,i))
			}
		}
	}
}

write.table(badnets,file='multicomponent_networks_to_rerun.tsv',sep='\t')



