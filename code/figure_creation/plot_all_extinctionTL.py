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
colors.add_color(17,119,51,'Tol3')
colors.add_color(51,34,136,'Tol4')
colors.add_color(221,204,119,'Tol5')
colors.add_color(153,153,51,'Tol6')
colors.add_color(204,142,159,'Tol7')
colors.add_color(136,34,85,'Tol8')
colors.add_color(170,68,153,'Tol9')
colors.add_color(221,221,221,'Tol10')

def read_data(meanfile):
  points={}
  f=open(meanfile,'r')
  for line in f:
    if line.split()[1]!='"S"':
      extorder=float(line.split()[1])
      S=int(line.split()[2])
      if S not in points:
        points[S]={}
      C=float(line.split()[3])
      if C not in points[S]:
        points[S][C]=[]
      STL=int(line.split()[4])+random.random()/2 -.25
      deg=int(line.split()[5])
      points[S][C].append((STL,10*extorder))
  return points


def format_graph(graph):
  # motdict={'6':'App. Comp.','36':'Dir. Comp.','38':'Omnivory','12':'Chain'}
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.xmax=8
  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.tick.configure(major=1,onoff='on',minor_ticks=0,major_size=.5,minor_size=.2,place='both',major_linewidth=1,minor_linewidth=1)
  graph.xaxis.label.configure(text='Shortest Trophic Level',char_size=1,just=2)

  graph.world.ymin=0
  graph.world.ymax=500
  graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.yaxis.tick.configure(major=100,onoff='on',minor_ticks=0,major_size=.5,minor_size=.2,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.label.configure(text='Mean time to extinction',char_size=1,just=2)

  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def populate_graph(graph,dat):
  for S in dat:
    for C in dat[S]:
      dots=graph.add_dataset(dat[S][C])
      dots.line.linestyle=0
      dots.symbol.configure(size=.1)

  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

dat=read_data('../stat_analyses/extinction_TLs.tsv')

grace=Grace(colors=colors)
graph=grace.add_graph()
graph=format_graph(graph)
graph=populate_graph(graph,dat)

# graph.set_view(0.1,0.45,0.9,0.95)
# grace.multi(rows=2,cols=2,vgap=.03,hgap=.05)
# grace.set_col_yaxislabel(col=0,rowspan=(None,None),label="Mean persistence time",char_size=1.5,just=2,position='normal')
# grace.set_row_xaxislabel(row=1,colspan=(None,None),label=altnames[normtype],char_size=1.5,just=2,position='normal')
# grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/extorder_vs_TL.eps')

