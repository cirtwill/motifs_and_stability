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
      # if p>0.00083:
      #   print s,c,p,'Non-significant'
      #   if p<0.05:
      #     print 'But close'
      if p<0.00083:
        print s,c,p,'Significant'
      subdict[(s,F)]=(c,p)
  # print ticker, ' non-significant Permanovas'
  f.close()

  return subdict

def format_graph(graph,form,yaxis):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  if yaxis=='F':
    graph.world.ymax=200
    graph.yaxis.tick.major=50
    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='pseudo-F statistic',char_size=1,just=2)
  else:
    graph.world.ymax=0.06
    graph.yaxis.tick.major=0.01
    graph.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
    graph.yaxis.label.configure(text='p-value',char_size=1,just=2)    

  graph.world.xmin=45
  graph.world.xmax=105
  graph.xaxis.tick.major=10
  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.xaxis.label.configure(text='Species richness',char_size=1,just=2)
  return graph

def populate_graph(graph,level,datadict):
  if level=='point':
    pointy=datadict['50'].keys()[0]
    dat=graph.add_dataset([pointy])
    dat.symbol.configure(fill_color=colorbar.z2color(datadict['50'][pointy][0]),shape=1,size=0.75,color=1)

  else:
    for S in datadict:
      for pointy in datadict[S]:
        dat=graph.add_dataset([pointy])
        dat.symbol.configure(fill_color=colorbar.z2color(datadict[S][pointy][0]),shape=1,size=0.75,color=1)

  return graph

def populate_pgraph(graph,level,datadict):
  if level=='full':
    for S in datadict:
      for pointy in datadict[S]:
        dat=graph.add_dataset([(S,datadict[S][pointy][1])])
        dat.symbol.configure(fill_color=colorbar.z2color(datadict[S][pointy][0]),shape=1,size=0.75,color=1)

  dashed=graph.add_dataset([(-1,0.05),(1000,0.05)])
  dashed.symbol.shape=0
  dashed.line.configure(linestyle=2,linewidth=1)
  return graph

        # color = colorbar.z2color(pdf)
        # # you can change the opacity percentage of a single color, as well
        # # color.change_opacity(60)
        # graph.add_dataset([(x0,y0), (x1,y1)], SolidRectangle, color)

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

for form in ['talk','paper']:
  for level in ['axis','point','full']:
    grace=MultiPanelGrace(colors=colors)
    grace.add_label_scheme('dummy',['','',''])
    grace.set_label_scheme('dummy')
    colorbar = grace.add_graph(ElLogColorBar,domain=(0,0.2),
                               scale=LINEAR_SCALE,autoscale=False)
    colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
    colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
    colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
    graph=grace.add_graph(Panel)
    graph=format_graph(graph,form,'F')
    pgraph=grace.add_graph(Panel)
    pgraph=format_graph(pgraph,form,'p')

    if level in ['point','full']:
      graph=populate_graph(graph,level,datadict)
      pgraph=populate_pgraph(pgraph,level,datadict)

    if form=='talk':
      graph.set_view(0.1,0.2,0.45,0.5)
      pgraph.set_view(0.55,0.2,0.9,0.5)
      colorbar.set_view(0.925,0.2,0.975,0.5)
    else:
      grace.add_label_scheme('smarty',['','A','B'])
      grace.set_label_scheme('smarty')
      for graphy in [graph,pgraph]:
        graphy.panel_label.configure(char_size=1.25,placement='iul',dy=0.02,dx=0.02)
      graph.set_view(0.4,0.55,0.9,0.9)
      graph.xaxis.label.text=''
      pgraph.set_view(0.4,0.15,0.9,0.5)
      colorbar.set_view(0.95,0.1,1.0,0.9)

    grace.write_file('../../manuscript/figures/extinction_order/permanova_summary_'+form+'_'+level+'.eps')
