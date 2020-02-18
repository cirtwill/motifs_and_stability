import sys
import os
import numpy as np

def read_motiffile(motiffile):
  motifdict={}
  f=open(motiffile,'r')
  for line in f:
    if len(line.split())>1:
      if line.split()[0]!='node':
        sp=line.split()[0][3:-1]
        role=line.split()[1:]
        motifdict[sp]=role
  f.close()
  return motifdict

def read_particfile(particfile):
  header=['S5', 'S1', 'D3', 'S4', 'S2', 'D1', 'D4', 'D8', 'S3', 'D5', 'D2', 'D7', 'D6']
  particdict={}
  f=open(particfile,'r')
  for line in f:
    if len(line.split())==14:
      sp=line.split()[0][3:-1]
      positions=line.split()[1:]
      particdict[sp]={}
      for i in range(0,len(header)):
        particdict[sp][header[i]]=float(positions[i])
    elif line.split() in [['S5', 'S1', 'D3', 'S4', 'S2', 'D1', 'D4', 'D8', 'S3', 'D5', 'D2', 'D7', 'D6'],[]]:
      pass
    else:
      print particfile
      print line.split()
      print "A wild header has appeared"
      sys.exit()
  f.close()
  return particdict

def read_degrees(degfile):
  degdict={}
  f=open(degfile,'r')
  for line in f:
    if len(line.split(','))==2:
      sp=line.split(',')[0][1:-1]
      val=line.split(',')[1].split()[0] # Remove pesky line breaks.
      degdict[sp]=val
  f.close()
  return degdict

def normalise_participation(particdict):
  normdict={}
  motiflists={}
  for motif in particdict['24']:
    motiflists[motif]=[particdict[species][motif] for species in particdict]

  for sp in particdict:
    normdict[sp]={}
    for motif in motiflists:
      mean=np.mean(motiflists[motif])
      sd=np.std(motiflists[motif])
      if mean!=0:
        dev=float(particdict[sp][motif]-mean)/sd
      else:
        dev=0
      normdict[sp][motif]=((particdict[sp][motif],dev))

  return normdict

def read_extorders(extfile):
  extdict={}
  # Dict of [removal][species][step]
  f=open(extfile,'r')
  for line in f:
    if 'Species' in line.split()[0]:
      remmy=line.split()[0].split('_')[1][:-1]
      extdict[remmy]={}
      i=1
      for sp in line.split()[1:]:
        extdict[remmy]['sp'+sp]=i
        i=i+1
    # Match to the extinction order for this removal
    elif 'Ext_step' in line.split()[0] and line.split()[0].split('_')[-1][:-1]==remmy: 
      for sp in extdict[remmy]:
        extdict[remmy][sp]='step_'+line.split()[extdict[remmy][sp]]
  f.close()
  return extdict

def write_matched_file(motifdict,extdict,outfile,s):
  f=open(outfile,'w')
  f.write('Species')
  for i in range(1,int(s)+1):
    f.write('\tRemoval_'+str(i))
  for j in range(1,31):
    f.write('\tPosition_'+str(j))
  f.write('\n')
  for species in sorted(motifdict):
    role=motifdict[species]
    f.write('sp'+str(species))
    for removal in sorted(extdict):
      step=extdict[str(removal)]['sp'+str(species)].split('_')[1]
      f.write('\t'+step)
    f.write('\t'+'\t'.join(role))
    f.write('\n')

  f.close()

def write_extorder_participation_file(particdict,degdict,TLdict,extdict,outfile,s):
  motifs=['S5', 'S1', 'D3', 'S4', 'S2', 'D1', 'D4', 'D8', 'S3', 'D5', 'D2', 'D7', 'D6']
  f=open(outfile,'w')
  f.write('Species')
  for i in range(1,int(s)+1):
    f.write('\tRemoval_'+str(i))
  f.write('\tSTL\tDegree')
  for motif in sorted(motifs):
    f.write('\tc'+motif+'\tz'+motif)
  f.write('\n')
  for species in sorted(particdict):
    f.write('sp'+str(species))
    for removal in sorted(extdict):
      step=extdict[str(removal)]['sp'+str(species)].split('_')[1]
      f.write('\t'+step) 
    f.write('\t'+str(TLdict[species])+'\t'+str(degdict['sp'+species]))
    for motif in sorted(motifs):
      (count,z)=particdict[species][motif]
      f.write('\t'+str(count)+'\t'+str(z))
    f.write('\n')
  f.close()

def main():

  motifdir='../data/roles/pre_disturbance/'
  particdir='../data/participation/pre_disturbance/'
  TLdir='../data/TLs/'
  degdir='../data/degrees/'

  for s in sorted(os.listdir(motifdir)):
    try:
      os.mkdir('../data/roles/matched_to_extorder/'+s)
    except OSError:
      pass
    for c in sorted(os.listdir(motifdir+s)):
    # for c in ['0.06','0.08','0.1']:
      try:
        os.mkdir('../data/roles/matched_to_extorder/'+s+'/'+c)
      except OSError:
        pass
      print s,c
      for j in range(0,100):
        motifdict=read_motiffile(motifdir+s+'/'+c+'/initial_edges_'+str(j)+'.roles')
        particdict=read_particfile(particdir+s+'/'+c+'/initial_edges_'+str(j)+'.participation')
        normdict=normalise_participation(particdict)
        degdict=read_degrees(degdir+s+'/'+c+'/degrees_'+str(j)+'.csv')
        TLdict=read_degrees(TLdir+s+'/'+c+'/SWTL_'+str(j)+'.csv')

        network=str(j)

        extfile='../data/networks/extorder/'+s+'/'+c+'/network_'+network+'_extinction_orders.tsv'
        extdict=read_extorders(extfile)

        outfile='../data/roles/matched_to_extorder/'+s+'/'+c+'/network_'+network+'.tsv'
        outfile2='../data/roles/matched_to_extorder/'+s+'/'+c+'/network_'+network+'_motifs.tsv'

        write_matched_file(motifdict,extdict,outfile,s)
        write_extorder_participation_file(normdict,degdict,TLdict,extdict,outfile2,s)


if __name__ == '__main__':
  main()
