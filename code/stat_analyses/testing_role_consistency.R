library(vegan)
library(nsprcomp)

# I think we want to do one permanova across all 100 networks for a size/connectance
# Testing: Do species with similar mean extinction orders have similar roles?
# 90 and 100-species ones don't run here, need to be on the mac

coordinates=matrix(nrow=0,ncol=6)
colnames(coordinates)=c("S","C","Position","Ax1","Ax2","Ax3")

stdevs=matrix(nrow=0,ncol=32)

for(S in seq(50,100,10)){
	print(S)
	for(C in seq(0.02,0.2,0.02)){

		axisfile=paste0('../../data/summaries/role_PCA/',as.character(S),'/role_PCA_axes_',as.character(S),'_',as.character(C),'.tsv',sep='')
		axes=read.table(axisfile,header=TRUE,sep='\t')

		devfile=paste0('../../data/summaries/role_PCA/',as.character(S),'/PCA_axis_var_explained_',as.character(S),'_',as.character(C),'.tsv',sep='')
		sdevs=read.table(devfile)
		devdat=c(S,C,sdevs$x)
		stdevs=rbind(stdevs,devdat)

		for(r in 1:nrow(axes)){
			rowdat=c(S,C,rownames(axes)[r],axes[r,1],axes[r,2],axes[r,3])
			coordinates=rbind(coordinates,rowdat)
		}
	}
}
print('Collecting done')
stdevs=as.data.frame(stdevs)
colnames(stdevs)[1:2]=c("S","C")
stdevs$Var1=stdevs[,3]/rowSums(stdevs[,3:32])
stdevs$Var2=stdevs[,4]/rowSums(stdevs[,3:32])
stdevs$Var3=stdevs[,5]/rowSums(stdevs[,3:32])


coordinates=as.data.frame(coordinates)
coordinates$S=as.numeric(as.character(coordinates$S))
coordinates$C=as.numeric(as.character(coordinates$C))
coordinates$Ax1=as.numeric(as.character(coordinates$Ax1))
coordinates$Ax2=as.numeric(as.character(coordinates$Ax2))
coordinates$Ax3=as.numeric(as.character(coordinates$Ax3))
meanload1=tapply(coordinates$Ax1,coordinates$Position,mean)
meanload2=tapply(coordinates$Ax2,coordinates$Position,mean)
meanload3=tapply(coordinates$Ax3,coordinates$Position,mean)

coords=cbind(meanload1,meanload2,meanload3)
write.table(coords,file='../figure_creation/mean_loadings_PCAaxes.tsv',sep='\t')
print('Permanova time')
# Want to see if the axes are consistent across webs.
coordist=vegdist(coordinates[,4:6],method="euclid")
# Do positions maintain consistent loadings on the three major axes?
distperm=adonis(coordist~coordinates$Position,permutations=9999)

print(mean(stdevs$Var1))
print(mean(stdevs$Var2))
print(mean(stdevs$Var3))
# How much variance do the first 3 axes explain across networks?

print(distperm)


