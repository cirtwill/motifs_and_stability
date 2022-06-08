import sys
import os
import math
import random
from decimal import *
import numpy as np

#Pygrace libraries
from PyGrace.grace import Grace
from PyGrace.colors import ColorBrewerScheme
from PyGrace.dataset import SYMBOLS
from PyGrace.Extensions.panel import Panel,MultiPanelGrace
from PyGrace.Extensions.distribution import CDFGraph, PDFGraph
from PyGrace.Extensions.latex_string import LatexString, CONVERT
from PyGrace.drawing_objects import DrawText, DrawLine,DrawBox
from PyGrace.Extensions.colorbar import SolidRectangle, ColorBar
from PyGrace.axis import LINEAR_SCALE, LOGARITHMIC_SCALE
from PyGrace.Styles.el import ElGraph, ElLogColorBar


colors=ColorBrewerScheme('RdYlBu',n=253,reverse=True)  # The blue is very beautiful but maybe harder to see.
# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

def read_file(infile):
  f=open(infile,'r')
  for line in f:
    if line.split()[0]=='"1"':
      rho=float(line.split()[1])
    elif line.split()[0]=='"2"':
      p=float(line.split()[1])
  # print ticker, ' non-significant Permanovas'
  f.close()

  return ((rho,p))

def format_graph(graph,yaxis,flavour):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=-1
  graph.world.ymax=1
  graph.yaxis.tick.major=.5
  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  if flavour=='freq':
    graph.yaxis.label.configure(text='Spearman rank correlation',char_size=1.25,just=2)

  graph.world.xmin=45
  graph.world.xmax=105
  graph.xaxis.tick.major=10
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  if flavour=='Z':
    graph.xaxis.label.configure(text='Species richness',char_size=1.25,just=2)
    graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  else:
    graph.xaxis.ticklabel.char_size=0
  return graph

def populate_graph(graph,datadict):
  for S in datadict:
    for C in datadict[S]:
      dat=graph.add_dataset([(int(S)+10*float(C),datadict[S][C][0])])
      dat.symbol.configure(fill_color=colorbar.z2color(float(C)),shape=1,size=0.5,color=1)
      if datadict[S][C][1]<0.05:
        dat.symbol.configure(shape=2,size=0.5,color=0)


  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################


# Dispersions have gone awry somehow.


datadir='../../data/mantel/'
datadict={}
for flavour in ['count','freq','Z']:
  datadict[flavour]={}
  # for s in ['50','60','70','80','90','100']:
  for s in ['50','100']:
    datadict[flavour][s]={}
    for c in ['0.02','0.04','0.06','0.08','0.1','0.12','0.14','0.16','0.18','0.2']:
    # for c in ['0.02','0.04','0.06','0.08','0.1','0.12','0.14']:
      datadict[flavour][s][c]=read_file(datadir+flavour+'/'+s+'/mantel_'+s+'_'+c+'.tsv')

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('smarty',['','A: count','B: frequency','C: Z-score'])
grace.set_label_scheme('smarty')
colorbar = grace.add_graph(ElLogColorBar,domain=(0,0.2),
                           scale=LINEAR_SCALE,autoscale=False)
colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
for flavour in datadict:
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,'F',flavour)
  graph=populate_graph(graph,datadict[flavour])

  graph.panel_label.configure(char_size=1,placement='iul',dy=0.02,dx=0.02)

grace.graphs[1].set_view(0.4,0.65,0.9,0.9)
grace.graphs[2].set_view(0.4,0.35,0.9,0.6)
grace.graphs[3].set_view(0.4,0.05,0.9,0.3)
# graph.xaxis.label.text=''
colorbar.set_view(0.95,0.05,1.0,0.9)

grace.write_file('../../manuscript/figures/mantel.eps')
