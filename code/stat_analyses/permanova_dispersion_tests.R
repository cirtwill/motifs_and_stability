library(vegan)

# Testing: Do species with similar mean extinction orders have similar roles?
# For each size-connectance combination, running a 9999-perm PERMANOVA stratified by network
# Also calculating an ANOVA of dispersion vs. quantile of extinction
# Also printing dispersion values vs. quantile for later plotting.
# Will need to be made human-readable with Python

for(method in c("count","freq","Z")){
  for(S in seq(50,100,10)){
    print(S)
    for(C in seq(0.02,0.2,0.02)){
      print(C)

      dists=read.table(file=paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C),'/dispersion_vs_quantile_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
      distmod=with(dists,lm(Dispersion~Quantile))
      distan=with(dists,aov(Dispersion~as.factor(Quantile)))
      distuk=TukeyHSD(distan,conf.level=.95)

      write.table(summary(distmod)$coefficients,file=paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C),'/LM_dispersion_vs_quantile_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')
      write.table(distuk$'as.factor(Quantile)',file=paste0('../../data/permanova/',method,'/disps/',as.character(S),'/',as.character(C),'/Tukey_dispersion_vs_quantile_',as.character(S),'_',as.character(C),'.tsv',sep=''),sep='\t')

    }
  }
}
# Test whether first extinction more similar to removed spp. than others
