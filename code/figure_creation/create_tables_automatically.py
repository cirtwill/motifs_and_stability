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

def read_dispfile(infile):
  subdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"x"':
      key=line.split()[0][1:-1]
      val=float(line.split()[1])
      if key in ['F','P','lm_slope','lm_pval']:
        if key in ['P','lm_pval']:
          val=np.round(val,3)
          if val<0.0001:
            val='\\textless0.001'
        if key=='F':
          if val<10 and val>1:
            val=np.round(val,2)
          else:  
            print val
        subdict[key]=val
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

def write_perm_table(datadict,outfile):
  f=open(outfile,'w')
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

def write_disp_table(datadict):
  f=open('../../manuscript/betadisper_summary_table.txt','w')
  f.write('S&C& ANOVA & & Regression & \\\\ \n')
  f.write('S&C&$F$&$p$-value&\\beta&$p$-value\\\\ \n')
  for s in sorted(datadict):
    for c in sorted(datadict[s]):
      dat=datadict[s][c]
      f.write(str(s)+'&'+str(c)+'&')
      f.write(str(dat['F'])+'&'+str(dat['P'])+'&')
      f.write(str(dat['lm_slope'])+'&'+str(dat['lm_pval'])+'\n')
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

datadir='../../data/summaries/permanova/'
for flavour in ['count','freq','Z']:
  outfile='../../manuscript/tables/permanova_summary_table_'+flavour+'.txt'
  datadict={}
  for sfil in sorted(os.listdir(datadir+flavour)):
    s=sfil.split('_')[-2]
    print s
    datadict[int(s)]=read_permfile(datadir+flavour+'/'+sfil)

  write_perm_table(datadict,outfile)

sys.exit()

dispdir='../../data/summaries/extorder_disps/'
dispdict={}
for s in sorted(os.listdir(dispdir)):
  dispdict[int(s)]={}
  for c in sorted(os.listdir(dispdir+'/'+s)):
    dispdict[int(s)][float(c)]=read_dispfile(dispdir+s+'/'+c+'/mean_extorder_vs_roles_'+s+'_'+c+'.tsv')

write_disp_table(dispdict)

lmdict={}
lmdir='../../data/summaries/rolesim_extorder/'
for s in ['50','60','70','80','90','100']:
  lmdict[int(s)]=read_lmfile(lmdir+s+'/role_exttime_lm_with_C.tsv')

write_lm_table(lmdict)
