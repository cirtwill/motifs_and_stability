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

def read_dispfile(infile):
  subdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"x"':
      key=line.split()[0][1:-1]
      val=float(line.split()[1])
      if key in ['F','P','lm_slope','lm_pval']:
        if key in ['P','lm_pval']:
          val=np.round(val,3)
          if val<0.0001:
            val='\\textless0.001'
        if key=='F':
          if val<10 and val>1:
            val=np.round(val,2)
          else:  
            print val
        subdict[key]=val
  # print ticker, ' non-significant Permanovas'
  f.close()

  return subdict

def format_graph(graph,form,yaxis,flavour):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  graph.world.ymax=200
  graph.yaxis.tick.major=50
  graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  if flavour=='freq':
    graph.yaxis.label.configure(text='pseudo-F statistic',char_size=1.25,just=2)

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

def populate_sgraph(graph,level,datadict):
  if level=='full':
    for S in datadict:
      for C in datadict[S]:
        slope=datadict[S][C]["lm_slope"]
        dat=graph.add_dataset([(S,slope)])
        dat.symbol.configure(fill_color=colorbar.z2color(C),shape=1,size=0.75,color=1)

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


# Dispersions have gone awry somehow.


datadir='../../data/summaries/permanova/'
datadict={}
# dispdict={}
for flavour in ['count','freq','Z']:
  datadict[flavour]={}
  # dispdict[flavour]={}
  for s in ['50','60','70','80','90','100']:
    datadict[flavour][s]=read_permfile(datadir+flavour+'/permanova_summary_'+s+'_'+flavour+'.tsv')

  # For each S in each flavour, [(s,F)]=(c,p)

#   dispdir='../../data/permanova/'+flavour+'/disps/'
#   for s in sorted(os.listdir(dispdir)):
#     dispdict[flavour][int(s)]={}
#     for c in sorted(os.listdir(dispdir+'/'+s)):
#       dispdict[flavour][int(s)][float(c)]=read_dispfile(dispdir+s+'/'+c+'/mean_extorder_vs_roles_'+s+'_'+c+'.tsv')

# print dispdict
# sys.exit()

level="full"
form=''

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
  graph=format_graph(graph,form,'F',flavour)
  # pgraph=grace.add_graph(Panel)
  # pgraph=format_graph(pgraph,form,'p')
  # sgraph=grace.add_graph(Panel)
  # sgraph=format_graph(sgraph,form,'slopes')

  graph=populate_graph(graph,level,datadict[flavour])
  # pgraph=populate_pgraph(pgraph,level,datadict)
  # sgraph=populate_sgraph(sgraph,level,dispdict)

  graph.panel_label.configure(char_size=1,placement='iul',dy=0.02,dx=0.02)

grace.graphs[1].set_view(0.4,0.65,0.9,0.9)
grace.graphs[2].set_view(0.4,0.35,0.9,0.6)
grace.graphs[3].set_view(0.4,0.05,0.9,0.3)
# graph.xaxis.label.text=''
colorbar.set_view(0.95,0.05,1.0,0.9)

grace.write_file('../../manuscript/figures/extinction_order/permanova_summary.eps')
