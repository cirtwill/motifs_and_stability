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

colors=ColorBrewerScheme('Reds',n=15,reverse=True)  # The blue is very beautiful but maybe harder to see.
colors.add_color(0,12,54,'blue2') 
colors.add_color(0,50,117,'blue3') 
colors.add_color(0,245,181,'blue4') 
colors.add_color(1,53,245,'blue5') 
colors.add_color(1,109,245,'blue6')  
# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

def read_file(infile):
  lmdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"Motif"':
      name=line.split()[1][2:-1]
      raw_beta=float(line.split()[2][1:-1])
      raw_p=float(line.split()[3][1:-1])
      prop_beta=float(line.split()[4][1:-1])
      prop_p=float(line.split()[5][1:-1])
      lmdict[name]={'raw':(raw_beta,raw_p),'prop':(prop_beta,prop_p)}
  f.close()
  return lmdict

def format_linegraph(graph,form):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.xmax=99
  graph.xaxis.tick.major=10

  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.label.configure(text='Degree',char_size=1,just=2)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  if form=='raw':
    graph.world.ymin=0
    graph.world.ymax=1000
    graph.yaxis.tick.major=200

    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='Count of motif',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  else:
    graph.world.ymin=0
    graph.world.ymax=1
    graph.yaxis.tick.major=.2

    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='Proportion of role',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='iul',dy=0.02,dx=0.03)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def populate_graph(graph,datadict,form):

  col=2
  sty=1
  for motif in ['S1','S2','S3','S4','S5','D1','D2','D3','D4','D5','D6','D7']:
    points=[]
    for x in range(0,100):
      points.append((x,x*datadict[motif][form][0]))
    pointy=graph.add_dataset(points)
    pointy.symbol.shape=0
    if motif[0]=='S':
      nucol='blue'+str(col)
    else:
      nucol=col
    pointy.line.configure(linestyle=sty,linewidth=1,color=nucol)
    col+=1
    sty+=1
    if sty>8:
      sty=1
    pointy.legend=motif

  # null=graph.add_dataset([(12,0)])
  # null.symbol.configure(shape=9,size=1,fill_color=4)

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

# scalefile='motif_lm_scales.tsv'
degfile='../../data/summaries/motifs_vs_degree.tsv'
degdict=read_file(degfile)

TLfile='../../data/summaries/motifs_vs_STL.tsv'
TLdict=read_file(TLfile)

grace=MultiPanelGrace(colors=colors)
# grace.add_label_scheme('dummy',['','S1: food chain','S2: omnivory','S4: direct competition','S5: apparent competition',''])
# grace.set_label_scheme('dummy')
# colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
#                            scale=LINEAR_SCALE,autoscale=False)
# colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
# colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
# colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
raw=grace.add_graph(Panel)
raw=format_linegraph(raw,'raw')
raw=populate_graph(raw,degdict,'raw')

prop=grace.add_graph(Panel)
prop=format_linegraph(prop,'prop')
raw=populate_graph(prop,degdict,'prop')


# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=2,cols=1,vgap=.03)
grace.hide_redundant_labels()
# # for graph in grace.graphs:
# #   print graph.get_view()
# colorbar.set_view(0.9,0.15,0.975,0.9)
# grace.graphs[1].set_view(0.125, 0.55, 0.475, 0.9)
# grace.graphs[2].set_view(0.525, 0.55, 0.875, 0.9)
# grace.graphs[3].set_view(0.125, 0.15, 0.475, 0.5)
# grace.graphs[4].set_view(0.525, 0.15, 0.875, 0.5)

grace.write_file('../../manuscript/figures/roles/motif_vs_degree.eps')
