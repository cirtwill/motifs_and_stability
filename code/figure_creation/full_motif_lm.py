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

colors=ColorBrewerScheme('RdYlBu',n=253,reverse=True)  # The blue is very beautiful but maybe harder to see.
# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

def process_lms(lmfile):
  # scaledict={}
  # f=open(scalefile,'r')
  # for line in f:
  #   if line.split()!=['"x"']:
  #     pred=line.split()[0][1:-1]
  #     scale=float(line.split()[1])
  #     scaledict[pred]=scale
  # f.close()

  lmdict={}
  f=open(lmfile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      name=line.split()[0][2:-1]
      shortname=name.split(')')[0]
      beta=float(line.split()[1])
      # if shortname!='Intercept':
      #   adjbeta=beta*scaledict[shortname]      
      # else:
      #   adjbeta=beta
      p=float(line.split()[-1])
      lmdict[shortname]=(beta,p)
  f.close()
  return lmdict

def process_ranges(rangefile):
  rangedict={}

  f=open(rangefile,'r')
  for line in f:
    if line.split()[0]=='"S1"':
      header=line.split()
    else:
      keyname=line.split()[0][1:-1]
      vals=line.split()[1:]
      for i in range(0,len(header)):
        motif=header[i]
        if motif[1:-1] not in rangedict:
          rangedict[motif[1:-1]]={}
        rangedict[motif[1:-1]][keyname]=float(vals[i])
  f.close()
  return rangedict

def format_linegraph(graph,form):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.world.xmax=13.5
  graph.xaxis.tick.major=1

  graph.xaxis.tick.set_spec_ticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[],tick_labels=['S1','S2','S3','S4','S5','D1','D2','D3','D4','D5','D6','D7','D8','S','C'])
  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.label.configure(text='Motif',char_size=1,just=2)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  if form=='coefs':
    graph.world.ymin=-.5
    graph.world.ymax=1
    graph.yaxis.tick.major=0.5

    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.configure(text='Coefficient',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  else:
    graph.world.ymin=-27
    graph.world.ymax=20
    graph.yaxis.tick.major=5

    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='Change in time to extinction',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='iul',dy=0.02,dx=0.03)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def populate_graph(graph,datadict):
  # Adding color blocks for the stable motifs
  block1=graph.add_dataset([(0.5,3),(5.5,3)])
  block1.fill.configure(color=70,type=2)
  block2=graph.add_dataset([(0.5,-3),(5.5,-3)])
  block2.fill.configure(color=70,type=2)
  block3=graph.add_dataset([(2.5,3),(3.5,3)])
  block3.fill.configure(color=0,type=2)
  block4=graph.add_dataset([(2.5,-3),(3.5,-3)])
  block4.fill.configure(color=0,type=2)

  x=1
  points=[]
  for motif in ['S1','S2','S3','S4','S5','D1','D2','D3','D4','D5','D6','D7','D8']:
    y=datadict[motif][0]
    points.append((x,y))
    x=x+1

  baseline=graph.add_dataset([(0,0),(100,0)])
  baseline.symbol.shape=0
  baseline.line.configure(linestyle=1,color=1,linewidth=0.5)

  dat=graph.add_dataset(points)
  dat.line.linestyle=0
  dat.symbol.configure(shape=1,size=.5,color=3,fill_color=3)

  # null=graph.add_dataset([(12,0)])
  # null.symbol.configure(shape=9,size=1,fill_color=4)

  # # Add an arrow for S
  # graph.add_drawing_object(DrawLine,end=(14,0),start=(14,-1.45),arrow=1,arrow_type=1,linewidth=2,linestyle=1,loctype='world')

  return graph

def populate_contextgraph(graph,coefdict,rangedict):

  # Adding color blocks for the stable motifs
  block1=graph.add_dataset([(0.5,50),(5.5,50)])
  block1.fill.configure(color=70,type=2)
  block2=graph.add_dataset([(0.5,-50),(5.5,-50)])
  block2.fill.configure(color=70,type=2)
  block3=graph.add_dataset([(2.5,50),(3.5,50)])
  block3.fill.configure(color=0,type=2)
  block4=graph.add_dataset([(2.5,-50),(3.5,-50)])
  block4.fill.configure(color=0,type=2)

  x=1
  maxes=[]
  means=[]
  for motif in ['S1','S2','S3','S4','S5','D1','D2','D3','D4','D5','D6','D7','D8']:
    coef=coefdict[motif][0]
    ranges=rangedict[motif]
    maxes.append((x,coef*ranges['max']))
    means.append((x,coef*ranges['mean']))
    x=x+1


  baseline=graph.add_dataset([(0,0),(100,0)])
  baseline.symbol.shape=0
  baseline.line.configure(linestyle=1,color=1,linewidth=0.5)

  dat=graph.add_dataset(maxes)
  dat.line.linestyle=0
  dat.symbol.configure(shape=3,size=.5,color=3,fill_color=3)
  dat.legend='Maximum'

  meany=graph.add_dataset(means)
  meany.line.linestyle=0
  meany.symbol.configure(shape=2,size=.5,color=3,fill_color=0)
  meany.legend='Mean'

  graph.legend.configure(box_linestyle=0,char_size=0.75,loc=(9,10),loctype='world')

  return graph


###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

# scalefile='motif_lm_scales.tsv'
lmfile='../../data/summaries/full_lm_coefficients.tsv'
datadict=process_lms(lmfile)

rangefile='../../data/summaries/motif_ranges.tsv'
rangedict=process_ranges(rangefile)

grace=MultiPanelGrace(colors=colors)
# grace.add_label_scheme('dummy',['','S1: food chain','S2: omnivory','S4: direct competition','S5: apparent competition',''])
# grace.set_label_scheme('dummy')
# colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
#                            scale=LINEAR_SCALE,autoscale=False)
# colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
# colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
# colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
graph=grace.add_graph(Panel)
graph=format_linegraph(graph,'coefs')
graph=populate_graph(graph,datadict)

context=grace.add_graph(Panel)
context=format_linegraph(context,'context')
context=populate_contextgraph(context,datadict,rangedict)

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

grace.write_file('../../manuscript/figures/extinction_order/motif_lmer_summary.eps')
