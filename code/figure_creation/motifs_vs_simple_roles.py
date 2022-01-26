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


# Do we also want a plot of extinction order vs. betas?

colors=ColorBrewerScheme('Spectral')  # The blue is very beautiful but maybe harder to see.
colors.add_color(136,204,238,'Tol1')
colors.add_color(122,170,153,'Tol2')
# colors.add_color(17,119,51,'Tol3')
colors.add_color(51,34,136,'Tol4')
# colors.add_color(221,204,119,'Tol5')
# colors.add_color(153,153,51,'Tol6')
colors.add_color(204,142,159,'Tol7')
colors.add_color(136,34,85,'Tol8')
colors.add_color(170,68,153,'Tol9')
# colors.add_color(221,221,221,'Tol10')

motif_means={'Count':{'Chain':152.89997,'Omni':151.80711,'DC':91.24644,'AC':274.19061,'Other':56.46391},
'Freq':{'Chain':0.23425806,'Omni':0.15057811,'DC':0.14015832,'AC':0.42344680,'Other':0.05155872},
'Z':{'Chain':3.317912e-16,'Omni':2.107330e-15,'DC':5.396885e-15,'AC':-2.344344e-16,'Other':-1.111925e-15}}

# Adding overall R2(fixed) to each panel
R2={'Deg':{'Count':0.926,'Freq':0.629,'Z':0.896},
'TL':{'Count':0.044,'Freq':0.030,'Z':0.021}}

def read_file(infile):
  lmdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      print line.split()
      name=line.split()[0][1:-1]
      if 'Intercept' in name:
        name='intercept'
      elif 'S5' in name:
        name='AC'
      elif 'S4' in name:
        name='DC'
      elif 'S2' in name:
        name='Omni'
      elif 'S1' in name:
        name='Chain'
      else:
        name='Other'
      effect=float(line.split()[1])
      pval=float(line.split()[-1])

      lmdict[name]=((effect,pval))

  return lmdict

def format_linegraph(graph,normtype,simple):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=1)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  if normtype=='Count':
    graph.world.xmax=4000
    graph.xaxis.tick.major=1000
    graph.xaxis.label.configure(text='Count',char_size=1.5,just=2)
  elif normtype=='Freq':
    graph.world.xmax=1
    graph.xaxis.tick.major=.2
    graph.xaxis.label.configure(text='Frequency',char_size=1.5,just=2)
    graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=1)
  else:
    graph.world.xmin=-3
    graph.world.xmax=10
    graph.xaxis.tick.major=2
    graph.xaxis.label.configure(text='Z-score',char_size=1.5,just=2)


  if simple=='Deg':
    graph.world.ymin=0
    graph.world.ymax=125
    graph.yaxis.tick.major=25
    graph.yaxis.label.configure(text='Degree',char_size=1.5,just=2)
  else:
    graph.world.ymin=0
    graph.world.ymax=3
    graph.yaxis.tick.major=1
    graph.yaxis.label.configure(text='STL',char_size=1.5,just=2)
  graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='iul',dy=0.02,dx=0.03)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def add_R2(graph,normtype,simple):
  if simple=='Deg':
    if normtype=='Count':
      x=2500
      y=25
    elif normtype=='Freq':
      x=0.6
      y=100
    else:
      x=5
      y=100
  else:
    if normtype=='Count':
      x=250
      y=.5
    elif normtype=='Freq':
      x=0.1
      y=.5
    else:
      x=-2
      y=.5

  graph.add_drawing_object(DrawText,text='R2='+str(R2[simple][normtype]),loctype='world',x=x,y=y,char_size=.75)

  return graph

def populate_graph(graph,datadict,normtype,simple):
  if normtype!='Freq':
    motifs=['AC','DC','Omni','Chain','Other']
  else:
    motifs=['AC','DC','Omni','Chain']
  # Actually I want to plot the marginal effects, not the whole thing...
  j=13
  if normtype=='Count':
    xs=range(1,4101,100)
  elif normtype=='Freq':
    bix=range(0,101,5)
    xs=[float(x)/100 for x in bix]
  else:
    bix=range(-5,21)
    xs=[float(x)/2 for x in bix]
  for motif in motifs:
    dats=[]
    for x in xs:
      motpart=datadict['intercept'][0]+x*datadict[motif][0]
      otherpart=sum([motif_means[normtype][mot]*datadict[mot][0] for mot in motifs if mot!=motif])
      dats.append((x,motpart+otherpart))

    pointy=graph.add_dataset(dats)
    pointy.symbol.shape=0
    pointy.line.configure(linestyle=1,linewidth=3,color=j)
    if motif=='Other':
      pointy.line.linestyle=2
    elif motif=='Omni':
      pointy.line.linestyle=3

    if normtype=='Z' and simple=='Deg':
      pointy.legend=motif
    j+=1

  # graph.add_drawing_object(DrawText,text='STL',x=90,y=16,loctype='world',just=0,char_size=1)
  if normtype=='Z' and simple=='Deg':
    graph.legend.configure(box_linestyle=0,char_size=.5,loc=(0,120),loctype='world')

  # # Add an arrow for S
  # graph.add_drawing_object(DrawLine,end=(14,0),start=(14,-1.45),arrow=1,arrow_type=1,linewidth=2,linestyle=1,loctype='world')

  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################


grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A','B','C','D','E','F'])
grace.set_label_scheme('dummy')
# colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
#                            scale=LINEAR_SCALE,autoscale=False)
# colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
# colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
# colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)

for simple in ['Deg','TL']:
  for normtype in ['Count','Freq','Z']:
    datfile='../../data/summaries/motifs_'+normtype+'_'+simple+'.tsv'
    datdict=read_file(datfile)
    # Need to update populate_graph
    graph=grace.add_graph(Panel)
    graph=format_linegraph(graph,normtype,simple)
    graph=add_R2(graph,normtype,simple)
    graph=populate_graph(graph,datdict,normtype,simple)


# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=2,cols=3,vgap=.03,hgap=.06)
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/motifs_vs_degTL.eps')

