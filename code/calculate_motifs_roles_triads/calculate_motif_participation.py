import sys
import os
import pymfinder as py

def main():

  size=3

  webdir='../../data/edgelists/pre_disturbance/'
  motifdir='../../data/participation/pre_disturbance/'

  for s in sorted(os.listdir(webdir)):
    try:
      os.mkdir(motifdir+s)
    except:
      pass
    for c in sorted(os.listdir(webdir+s)):
      print s,c
      try:
        os.mkdir(motifdir+s+'/'+c)
      except:
        pass
      for web in sorted(os.listdir(webdir+s+'/'+c)):
        infile=webdir+s+'/'+c+'/'+web

        motfile=motifdir+s+'/'+c+'/'+web.split('.tsv')[0]+'.participation'

        result=py.motif_participation(infile,motifsize=3,stoufferIDs=True,allmotifs=True)
        motifs=str(result).split('node')[1]

        g=open(motfile,'w')
        g.write(motifs)
        g.close()

if __name__ == '__main__':
  main()
