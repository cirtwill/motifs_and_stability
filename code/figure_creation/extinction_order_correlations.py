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

# For correlations
modeldict={'S':0.001349,'C':-0.05893,'S:C':0.001963,'intercept':0.7858}

# For means
def read_coeffile(coeffile):
  dicto={}
  f=open(coeffile)
  for line in f:
    if line.split()[0]!='"Estimate"':
      pred=line.split()[0]
      if '$' in pred and len(pred.split('$'))==2:
        pred=pred.split('$')[1][:-1]
      elif '$' in pred and len(pred.split('$'))>2:
        pred='S:C'
      else:
        pred='Intercept'
      coef=float(line.split()[1])

      dicto[pred]=coef
  f.close()

  return dicto

def read_permfile(infile):
  subdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='S':
      s=int(line.split()[0])
      c=float(line.split()[1])
      ext_corr=float(line.split()[2])
      F=float(line.split()[6])    
      p=float(line.split()[7])
      subdict[(s,ext_corr)]=c
  f.close()

  return subdict

def format_graph(graph,yaxis):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.panel_label.configure(char_size=0.75,placement='iul',dx=0.03,dy=0.03)
  graph.world.ymin=0.75
  graph.world.xmin=45
  graph.world.xmax=105
  graph.xaxis.tick.major=10
  if yaxis=='dots':
    graph.world.ymax=1
    graph.yaxis.tick.major=0.1
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.configure(text=\
      LatexString(r'R\S2\N time to extinction'),char_size=.75,just=2)
    graph.xaxis.ticklabel.char_size=0
  elif yaxis=='means': 
    graph.world.ymin=0
    graph.world.ymax=50
    graph.yaxis.tick.major=10
    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='Mean time to extinction (x10)',char_size=.75,just=2)   
    graph.xaxis.label.configure(text='Species richness',char_size=1,just=2)
    graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  elif yaxis=='Cmeans': 
    graph.world.ymin=0
    graph.world.ymax=50
    graph.yaxis.tick.major=10
    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='Mean time to extinction (x10)',char_size=.75,just=2)   
    graph.xaxis.label.configure(text='Connectance',char_size=1,just=2)
    graph.xaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
    graph.world.xmin=0.0001
    graph.world.xmax=0.201
    graph.xaxis.tick.major=.04
  else:
    graph.world.ymax=0.06
    graph.yaxis.tick.major=0.01
    graph.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
    graph.yaxis.label.configure(text='p-value',char_size=1,just=2)   


  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  return graph

def populate_graph(graph,level,datadict):
  if level=='point':
    pointy=sorted(datadict['50'].keys())[0]
    dat=graph.add_dataset([pointy])
    dat.symbol.configure(fill_color=colorbar.z2color(datadict['50'][pointy]),shape=1,size=0.75,color=1)

  else:
    for S in datadict:
      for pointy in sorted(datadict[S]):
        dat=graph.add_dataset([pointy])
        dat.symbol.configure(fill_color=colorbar.z2color(datadict[S][pointy]),shape=1,size=0.75,color=1)

  return graph

def fill_predictions(graph,level):
  if level=='point':
    Cs=[0.02]
  else:
    Cs=[0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2]
  for C in Cs:
    datapoints=[]
    for S in range(50,101):
      pred=modeldict['intercept']+S*modeldict['S']+C*modeldict['C']+S*C*modeldict['S:C']
      datapoints.append((S,pred))
    dats=graph.add_dataset(datapoints)
    dats.symbol.shape=0
    dats.line.configure(color=colorbar.z2color(C),linewidth=1)

  return graph

def fill_Sgraph(graph,modeldict):
  Cs=[0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2]
  for C in Cs:
    datapoints=[]
    for S in range(50,101):
      pred=modeldict['Intercept']+S*modeldict['S']+C*modeldict['C']+S*C*modeldict['S:C']
      datapoints.append((S,pred))
    dats=graph.add_dataset(datapoints)
    dats.symbol.shape=0
    dats.line.configure(color=colorbar.z2color(C),linewidth=1)

  return graph

def fill_Cgraph(graph,modeldict):
  Ss=[50,60,70,80,90,100]
  for S in Ss:
    print S, colorbar2.z2color(S)
    datapoints=[]
    for C in [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2]:
      pred=modeldict['Intercept']+S*modeldict['S']+C*modeldict['C']+S*C*modeldict['S:C']
      datapoints.append((C,pred))
    dats=graph.add_dataset(datapoints)
    dats.symbol.shape=0
    dats.line.configure(color=colorbar2.z2color(S),linewidth=1)

  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

datadir='../../data/summaries/extorder_perms/'
datadict={}
for s in ['50','60','70','80','90','100']:
  datadict[s]=read_permfile(datadir+s+'/extorder_roles_permanova_summary_'+s+'.tsv')
form='paper'
level='full'

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['','A','B','','C'])
grace.set_label_scheme('dummy')
colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
                           scale=LINEAR_SCALE,autoscale=False)
colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)

# Correlation graph
graph=grace.add_graph(Panel)
graph=format_graph(graph,'dots')

graph=fill_predictions(graph,level)
graph=populate_graph(graph,level,datadict)

# for graph in grace.graphs:
#   print graph.get_view()

lmcoefs=read_coeffile('../stat_analyses/mean_exttime_lm.tsv')

# Mean time to extinction vs. S
Sgraph=grace.add_graph(Panel)
Sgraph=format_graph(Sgraph,'means')
Sgraph=fill_Sgraph(Sgraph,lmcoefs)

colorbar2 = grace.add_graph(ElLogColorBar,domain=(50,100),
                           scale=LINEAR_SCALE,autoscale=False)
colorbar2.yaxis.label.configure(text="Species Richness",char_size=1,just=2)
colorbar2.yaxis.tick.configure(major=10,major_size=.5,minor_ticks=0)
colorbar2.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)

# Mean time to extinction vs. C
Cgraph=grace.add_graph(Panel)
Cgraph=format_graph(Cgraph,'Cmeans')
Cgraph=fill_Cgraph(Cgraph,lmcoefs)



graph.set_view(0.15,0.7,0.6,0.95)
colorbar.set_view(0.625,0.43,0.675,0.95)
Sgraph.set_view(0.15,0.43,0.6,0.68)
Cgraph.set_view(0.15,0.10,0.6,0.35)
colorbar2.set_view(0.625,0.10,0.675,0.35)

grace.write_file('../../manuscript/figures/extinction_order/extorder_correlations.eps')
