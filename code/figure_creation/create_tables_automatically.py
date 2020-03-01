import os
import sys
import numpy as np

def read_permfile(infile):
  subdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='S':
      s=int(line.split()[0])
      c=float(line.split()[1])
      ext_corr=np.round(float(line.split()[2]),3)
      F=str(line.split()[6])    
      if len(F.split('.')[0])==3:
        F=int(np.round(float(F)))
      elif len(F.split('.')[0])==2:
        F=np.round(float(F),1)
      elif len(F.split('.')[0])==1:
        F=np.round(float(F),2)
      else:
        print F
      p=np.round(float(line.split()[7]),3)
      if p<0.00083:
        print s,c,p,'Significant'
      subdict[c]=(ext_corr,F,p)
  # print ticker, ' non-significant Permanovas'
  f.close()

  return subdict

def read_lmfile(infile):
  subdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      name=line.split()[0][1:-1]
      beta=line.split()[1]
      if beta[0]=='-' and float(beta)>-1:
        beta=np.round(float(beta),3)
      elif beta[0]=='-' and float(beta)<-1:
        if len(beta.split('.')[0])==3:
          beta=int(np.round(float(beta),1))
        elif len(beta.split('.')[0])==2:
          beta=np.round(float(beta),2)
        elif len(F.split('.')[0])==1:
          beta=np.round(float(beta),3)
        else:
          print beta
      else:
        if len(beta.split('.')[0])==3:
          beta=int(np.round(float(beta)))
        elif len(beta.split('.')[0])==2:
          beta=np.round(float(beta),1)
        elif len(beta.split('.')[0])==1 and float(beta)>1:
          beta=np.round(float(beta),2)
        elif len(beta.split('.')[0])==1 and float(beta)<1:
          beta=np.round(float(beta),3)
        else:
          print beta
      p=np.round(float(line.split()[-1]),3)
      subdict[name]=(beta,p)
  f.close()

  return subdict

def write_perm_table(datadict):
  f=open('../../manuscript/permanova_summary_table.txt','w')
  f.write('S&C&Extinction correlation&pseudo-$F$&$p$-value\\\\ \n')
  for s in sorted(datadict):
    for c in sorted(datadict[s]):
      dat=datadict[s][c]
      f.write(str(s)+'&'+str(c)+'&')
      if len(str(dat[0]))==5:
        f.write(str(dat[0])+'&')
      elif len(str(dat[0]))==4:
        f.write(str(dat[0])+'0&')
      elif len(str(dat[0]))==3:
        f.write(str(dat[0])+'00&')
      else:
        print len(str(dat[0])), str(dat[0])
      f.write(str(dat[1])+'&')
      if len(str(dat[2]))==5:
        f.write(str(dat[2])+'\\\\ \n')
      elif len(str(dat[2]))==4:
        f.write(str(dat[2])+'0\\\\ \n')
      else:
        print len(str(dat[2]))
  f.close()

def write_lm_table(datadict):
  f=open('../../manuscript/lmer_summary_table.txt','w')
  f.write('Term&S&Beta&$p$-value\\\\ \n')
  for s in sorted(datadict):
    for term in sorted(datadict[s]):   
      dat=datadict[s][term]
      f.write(term+'&'+str(s)+'&'+str(dat[0])+'&'+str(dat[1])+'\\\\ \n') 
  f.close()


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

datadir='../../data/summaries/extorder_perms/'
datadict={}
for s in sorted(os.listdir(datadir)):
  datadict[int(s)]=read_permfile(datadir+s+'/extorder_roles_permanova_summary_'+s+'.tsv')

write_perm_table(datadict)

sys.exit()

lmdict={}
lmdir='../../data/summaries/rolesim_extorder/'
for s in ['50','60','70','80','90','100']:
  lmdict[int(s)]=read_lmfile(lmdir+s+'/role_exttime_lm_with_C.tsv')

write_lm_table(lmdict)
