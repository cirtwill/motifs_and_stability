import sys
import os

def permfile_reader(permfile,request):
	temp={}
	f=open(permfile,'r')
	for line in f:
		if line.split()[0]!='"S"':
			S=int(line.split()[1])
			temp['S']=S
			C=float(line.split()[2])
			temp['C']=C
			ext_corr=float(line.split()[3])
			temp['ext_corr']=ext_corr
			df=int(line.split()[4])
			temp['df']=df
			SS=float(line.split()[5])
			temp['SS']=SS
			MS=float(line.split()[6])
			temp['MS']=MS
			model_F=float(line.split()[7])
			temp['model_F']=model_F
			R2=float(line.split()[8])
			temp['R2']=R2
			pval=line.split()[9]
			perms=line.split()[10:]
	f.close()

	permuts=[float(x) for x in perms]


	if request=="observed":
		return temp
	elif request=="random":
		return permuts

def pval_calculator(model_F,perms):
	greater=[x for x in perms if x>=model_F]

	p=float(len(greater))/float(len(perms))

	return p

def recreate_aovtab(datadict,outfile):

	f=open(outfile,'w')
	f.write('S\tC\text_corr\tdf\tSS\tMS\tModel_F\tR2\tpval\n')
	for C in sorted(datadict):
		dats=datadict[C]
		f.write(str(dats['S'])+'\t'+str(dats['C'])+'\t'+str(dats['ext_corr'])+'\t')
		f.write(str(dats['df'])+'\t'+str(dats['SS'])+'\t'+str(dats['MS'])+'\t')
		f.write(str(dats['model_F'])+'\t'+str(dats['R2'])+'\t'+str(dats['pval'])+'\n')
	f.close()




def main():

  motifdir='../../data/summaries/extorder_perms/'

  for s in sorted(os.listdir(motifdir)):  	
    outfile='../../data/summaries/extorder_perms/'+s+'/extorder_roles_permanova_summary_'+s+'.tsv'
    datadict={}
    for c in [x for x in sorted(os.listdir(motifdir+s)) if 'summary_' not in x]:
      filelist=sorted(os.listdir(motifdir+s+'/'+c))
      permfiles=[x for x in filelist if 'mean_extorder_vs_roles' in x]

      datadict[c]=permfile_reader(motifdir+s+'/'+c+'/'+permfiles[0],"observed")
      all_perms=[]

      for permfile in permfiles:
		perms=permfile_reader(motifdir+s+'/'+c+'/'+permfile,"random")      	
		for perm in perms:
			all_perms.append(perm)

      if len(all_perms)==9999:
      	datadict[c]['pval']=pval_calculator(datadict[c]['model_F'],all_perms)
      	print 'All present and correct with ',s,c
      else:
      	print 'Missing ',str(9999-len(all_perms)),' perms for run: ',s,' ',c

    recreate_aovtab(datadict,outfile)

if __name__ == '__main__':
  main()
