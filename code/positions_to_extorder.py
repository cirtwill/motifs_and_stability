import sys
import os
import numpy as np

def read_particfile(particfile):
  particdict={}
  header=[]
  f=open(particfile,'r')
  for line in f:
    if len(line.split())>0:
      if line.split()[0]=='node':
        for pos in line.split('(')[1:]:
          (mot,pred,prey)=pos.split(')')[0].split(',')
          header.append((mot,pred,prey))
      else:
        sp=line.split()[0][3:-1]
        positions=line.split()[1:]
        particdict[sp]={}
        for i in range(0,len(header)):
          particdict[sp][header[i]]=float(positions[i])
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

def normalise_participation(particdict,normstyle):
  normdict={}
  motiflists={}
  grouped_motifs={}
  ungroup=['6','36','38','12'] # All other positions are in loopies
  for (motif,preds,prey) in particdict['24']:
    motiflists[(motif,preds,prey)]=[particdict[species][(motif,preds,prey)] for species in particdict]
  for species in particdict:
    groupies=[]
    for (motif,preds,prey) in particdict[species]:
      if motif not in ungroup:
        groupies.append(particdict[species][(motif,preds,prey)])
    grouped_motifs[species]=sum(groupies)

  for sp in particdict:
    normdict[sp]={}
    if normstyle=='full':
      for motif in motiflists:
        mean=np.mean(motiflists[motif])
        sd=np.std(motiflists[motif])
        if mean!=0:
          dev=float(particdict[sp][motif]-mean)/sd
        else:
          dev=0
        normdict[sp][motif]=((particdict[sp][motif],dev))
    elif normstyle=='condensed':
      for motif in [(m,p1,p2) for (m,p1,p2) in motiflists if m in ungroup]:
        mean=np.mean(motiflists[motif])
        sd=np.std(motiflists[motif])
        if mean!=0:
          dev=float(particdict[sp][motif]-mean)/sd
        else:
          dev=0
        normdict[sp][motif]=((particdict[sp][motif],dev))
      # Now do the `other' motifs (loop-containing)
      netmean=np.mean(grouped_motifs.values())
      netsd=np.std(grouped_motifs.values())
      obs=grouped_motifs[sp]
      if netmean!=0:
        normdict[sp]['other']=((obs,float(obs-netmean)/netsd))
      else:
        normdict[sp]['other']=((obs,0))

  return normdict
# You are calculating Z-scores for summed 'other' motifs
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

# Printing counts and Z-scores
def write_extorder_participation_file(particdict,degdict,TLdict,extdict,outfile,s):
  motifs=[m+':'+p1.split()[0]+':'+p2.split()[0] for (m,p1,p2) in particdict['24'].keys()]

  f=open(outfile,'w')
  f.write('Species')
  for i in range(1,int(s)+1):
    f.write('\tRemoval_'+str(i))
  f.write('\tSTL\tDegree')
  for motif in sorted(motifs):
    f.write('\tc'+motif)
  for motif in sorted(motifs):
    f.write('\tz'+motif)    
  f.write('\n')
  for species in sorted(particdict):
    f.write('sp'+str(species))
    for removal in sorted(extdict):
      step=extdict[str(removal)]['sp'+str(species)].split('_')[1]
      f.write('\t'+step) 
    f.write('\t'+str(TLdict[species])+'\t'+str(degdict['sp'+species]))
    for m in sorted(particdict[species].keys()):
      (count,z)=particdict[species][m]
      f.write('\t'+str(count))
    for m in sorted(particdict[species].keys()):
      (count,z)=particdict[species][m]
      f.write('\t'+str(z))
    f.write('\n')
  f.close()

def write_extorder_condensed_file(particdict,degdict,TLdict,extdict,outfile,s):
  pos=['other']
  for (m,p1,p2) in [k for k in particdict['24'].keys() if k!='other']:
    pos.append(m+':'+p1.split()[0]+':'+p2.split()[0])

  f=open(outfile,'w')
  f.write('Species')
  for i in range(1,int(s)+1):
    f.write('\tRemoval_'+str(i))
  f.write('\tSTL\tDegree')
  for p in sorted(pos):
    f.write('\tc'+p)
  for p in sorted(pos):
    f.write('\tz'+p)    
  f.write('\n')
  for species in sorted(particdict):
    f.write('sp'+str(species))
    for removal in sorted(extdict):
      step=extdict[str(removal)]['sp'+str(species)].split('_')[1]
      f.write('\t'+step) 
    f.write('\t'+str(TLdict[species])+'\t'+str(degdict['sp'+species]))
    for motif in sorted(particdict[species]):
      (count,z)=particdict[species][motif]
      f.write('\t'+str(count))
    for motif in sorted(particdict[species]):
      (count,z)=particdict[species][motif]
      f.write('\t'+str(z))
    f.write('\n')
  f.close()


def main():

  particdir='../data/positions/pre_disturbance/'
  TLdir='../data/TLs/'
  degdir='../data/degrees/'

  for s in sorted(os.listdir(particdir)):
    try:
      os.mkdir('../data/positions/matched_to_extorder/'+s)
      os.mkdir('../data/positions/matched_to_extorder_condensed/'+s)
    except OSError:
      pass
    for c in sorted(os.listdir(particdir+s)):
    # for c in ['0.06','0.08','0.1']:
      try:
        os.mkdir('../data/positions/matched_to_extorder/'+s+'/'+c)
        os.mkdir('../data/positions/matched_to_extorder_condensed/'+s+'/'+c)
      except OSError:
        pass
      print s,c
      for j in range(0,100):
        particdict=read_particfile(particdir+s+'/'+c+'/initial_edges_'+str(j)+'.roles')
        normdict=normalise_participation(particdict,'full')
        condensed_normdict=normalise_participation(particdict,'condensed')
        degdict=read_degrees(degdir+s+'/'+c+'/degrees_'+str(j)+'.csv')
        TLdict=read_degrees(TLdir+s+'/'+c+'/SWTL_'+str(j)+'.csv')

        network=str(j)

        extfile='../data/networks/extorder/'+s+'/'+c+'/network_'+network+'_extinction_orders.tsv'
        extdict=read_extorders(extfile)

        outfile='../data/positions/matched_to_extorder/'+s+'/'+c+'/network_'+network+'.tsv'
        outfile2='../data/positions/matched_to_extorder_condensed/'+s+'/'+c+'/network_'+network+'.tsv'

        write_extorder_participation_file(normdict,degdict,TLdict,extdict,outfile,s)

        write_extorder_condensed_file(condensed_normdict,degdict,TLdict,extdict,outfile2,s)

if __name__ == '__main__':
  main()
