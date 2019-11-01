import sys
import os
import pymfinder as py

def main():

  size=3

  webdir='../../data/edgelists/pre_disturbance/'
  motifdir='../../data/motifs/pre_disturbance/'
  roledir='../../data/roles/pre_disturbance/'
  pardir='../../data/participation/pre_disturbance/'

  for s in sorted(os.listdir(webdir)):
    try:
      os.mkdir(roledir+s)
    except:
      pass
    try:
      os.mkdir(motifdir+s)
    except:
      pass
    for c in sorted(os.listdir(webdir+s)):
      print s,c
      try:
        os.mkdir(roledir+s+'/'+c)
      except:
        pass
      try:
        os.mkdir(motifdir+s+'/'+c)
      except:
        pass
      for web in sorted(os.listdir(webdir+s+'/'+c)):
        infile=webdir+s+'/'+c+'/'+web

        rolefile=roledir+s+'/'+c+'/'+web.split('.tsv')[0]+'.roles'
        motfile=motifdir+s+'/'+c+'/'+web.split('.tsv')[0]+'.motifs'
        parfile=pardir+s+'/'+c+'/'+web.split('.tsv')[0]+'.participation'

        result=py.motif_roles(infile,motifsize=size,stoufferIDs=True,allroles=True)
        motifs=str(result).split('node')[0]
        roles='node'+str(result).split('node')[1]


        presult=py.motif_participation(infile,motifsize=3,stoufferIDs=True,allmotifs=True)
        partic=str(presult).split('node')[1]

        e=open(parfile,'w')
        e.write(partic)
        e.close()


        f=open(rolefile,'w')
        f.write(roles)
        f.close()
        g=open(motfile,'w')
        g.write(motifs)
        g.close()


if __name__ == '__main__':
  main()
