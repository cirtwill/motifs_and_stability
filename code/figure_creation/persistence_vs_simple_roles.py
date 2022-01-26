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

colors=ColorBrewerScheme('RdBu',n=9)  # The blue is very beautiful but maybe harder to see.
colors.add_color(136,204,238,'Tol1')
colors.add_color(122,170,153,'Tol2')
colors.add_color(17,119,51,'Tol3')
colors.add_color(51,34,136,'Tol4')
colors.add_color(221,204,119,'Tol5')
colors.add_color(153,153,51,'Tol6')
colors.add_color(204,142,159,'Tol7')
colors.add_color(136,34,85,'Tol8')
colors.add_color(170,68,153,'Tol9')
colors.add_color(221,221,221,'Tol10')
colors.add_color(220,220,220,'lightgrey')

def read_file(infile):
  lmdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      name=line.split()[0][1:-1]
      if 'Intercept' in name:
        name='intercept'
      else:
        nam=name.split('allTLs$')
        name=''.join(nam)
      effect=float(line.split()[1])
      pval=float(line.split()[-1])

      lmdict[name]=((effect,pval))

  f.close()
  return lmdict

def format_linegraph(graph):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  # graph.frame.configure(background_color='lightgrey',background_pattern=1)

  graph.world.xmin=0
  graph.world.xmax=125
  graph.xaxis.tick.major=25
  graph.world.ymin=0
  graph.world.ymax=50
  graph.yaxis.tick.major=10

  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=1)
  graph.xaxis.label.configure(text='Degree',char_size=1.5,just=2)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=1)
  graph.yaxis.label.configure(text='Mean persistence time',char_size=1.5,just=2)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  # graph.panel_label.configure(char_size=0.75,placement='iul',dy=0.02,dx=0.03)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def populate_graph(graph,datadict):

  j=2
  for TL in [1,2,3,4,5,6]:
    dats=[]
    for deg in range(1,125):
      x=deg
      y=datadict['intercept'][0]+TL*datadict['STL'][0]+deg*datadict['Deg'][0]+deg*TL*datadict['Deg:STL'][0]
      dats.append((x,y))

    pointy=graph.add_dataset(dats)
    pointy.symbol.shape=0
    if j==6:
      pointy.line.configure(linestyle=1,linewidth=1,color=1)
    elif j==12:
      pointy.line.configure(linestyle=2,linewidth=1,color=1)      
    else:
      pointy.line.configure(linestyle=1,linewidth=1,color=j)

    pointy.legend=str(TL)
    j+=2

  graph.add_drawing_object(DrawText,text='STL',x=90,y=16,loctype='world',just=0,char_size=1)
  graph.legend.configure(box_linestyle=0,char_size=1,loc=(90,15),loctype='world')

  # # Add an arrow for S
  # graph.add_drawing_object(DrawLine,end=(14,0),start=(14,-1.45),arrow=1,arrow_type=1,linewidth=2,linestyle=1,loctype='world')
  graph.add_drawing_object(DrawText,text='R2=0.211',x=90,y=40,loctype='world',just=0,char_size=1)

  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

datfile='../../data/summaries/persistence_degTL.tsv'
datdict=read_file(datfile)

grace=Grace(colors=colors)
# grace.add_label_scheme('dummy',['A','B','C'])
# grace.set_label_scheme('dummy')
# colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
#                            scale=LINEAR_SCALE,autoscale=False)
# colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
# colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
# colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)

graph=grace.add_graph()
graph=format_linegraph(graph)
graph=populate_graph(graph,datdict)


# graph.set_view(0.1,0.45,0.9,0.95)
# grace.multi(rows=3,cols=1,vgap=.03)
# grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/roles/persistence_vs_degTL.eps')

