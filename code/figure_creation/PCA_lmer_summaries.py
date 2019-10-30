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

colors=ColorBrewerScheme('RdBu',n=253,reverse=True)  # The blue is very beautiful but maybe harder to see.
# colors.add_color(120,120,120,'grey')
# colors.add_color(255,125,125,'lightish_red')
# colors.add_color(200,200,200,'lightgrey')

def read_scalefile(scalefile):
  scales={}
  f=open(scalefile,'r')
  for line in f:
    if line.split()[0]!='"Val"':
      pred=line.split()[1][1:-1]
      center=float(line.split()[2][1:-1])
      scale=float(line.split()[3][1:-1])
      scales[pred]=(center,scale)
  f.close()

  return scales

def read_lmfile(infile):
  lmdict={}
  f=open(infile,'r')
  for line in f:
    # print line.split()
    if line.split()[0]!='"Estimate"':
      name=line.split()[0][1:-1]
      if len(name.split(':'))==1:
        shortname=name.split('(')[1].split(')')[0]
      else:
        comps=[]
        for comp in name.split(':'):
          shorty=comp.split('(')[1].split(')')[0]
          comps.append(shorty)
        shortname=':'.join(comps)
      beta=float(line.split()[1])      
      p=float(line.split()[-1])
      lmdict[shortname]=(beta,p)
  f.close()

  return lmdict

def format_linegraph(graph,form,focus):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  graph.world.xmin=-1300
  graph.world.xmax=4000
  graph.world.ymax=50
  graph.yaxis.tick.major=10
  graph.xaxis.tick.major=1000
  graph.yaxis.ticklabel.char_size=0
  graph.xaxis.ticklabel.char_size=0

  if form=='paper':
    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    if focus in ['PC2']:
      graph.yaxis.label.configure(text='Time to extinction',char_size=1,just=2)
    if focus in ['PC3']:
      graph.xaxis.ticklabel.configure(format='decimal',prec=0,angle=30,char_size=.75)
      graph.xaxis.label.configure(text='Position on axis',char_size=1,just=2)
  else:
    graph.xaxis.ticklabel.configure(format='decimal',prec=0,angle=90,char_size=.75)
    if focus in ['PC1']:
      graph.yaxis.label.configure(text='Time to extinction',char_size=1,just=2)
      graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    if focus in ['PC2']:
      graph.xaxis.label.configure(text='Position on axis',char_size=1,just=2)

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.4,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.panel_label.configure(char_size=0.75,placement='ouc',dx=0.02,dy=0.005)

  if form=='paper':
    graph.legend.configure(char_size=0.5,box_linestyle=0,loctype='world',loc=(0,25))
  else:
    graph.legend.configure(char_size=0.5,box_linestyle=0,loctype='world',loc=(-400,22))

  return graph

def populate_graph(graph,level,datadict,scaledict,focus):
  if level=='point':
    Ss=[50,70,100]
    Cs=[0.02]
  else:
    Ss=[50,70,100]
    Cs=[0.02,0.1,0.2]

  coefs=datadict[focus]
  for S in Ss:
    scale_s=(S-scaledict['S'][0])/scaledict['S'][1]
    if S==50:
      typ=2
    elif S==70:
      typ=3
    else:
      typ=1
    for C in Cs:
      col=colorbar.z2color(C)
      scale_c=(C-scaledict['C'][0])/scaledict['C'][1]
      baseline=datadict[focus]['Intercept'][0]+datadict[focus]['S'][0]*scale_s+datadict[focus]['C'][0]*scale_c+datadict[focus]['S:C'][0]*scale_s*scale_c

      points=[]
      if level in ['point','single'] and focus=='PC1':
        for count in range(-2000,4000):
          scale_count=(count-scaledict[focus][0])/scaledict[focus][1]
          Spart=datadict[focus][focus][0]*scale_count+datadict[focus][focus+':S'][0]*scale_s*scale_count
          Cpart=datadict[focus][focus+':C'][0]*scale_c*scale_count+datadict[focus][focus+':S:C'][0]*scale_count*scale_s*scale_c
          y=baseline+Spart+Cpart
          points.append((count,y))
      elif level=='two' and focus in ['PC1','PC2']:
        for count in range(-2000,4000):
          scale_count=(count-scaledict[focus][0])/scaledict[focus][1]
          Spart=datadict[focus][focus][0]*scale_count+datadict[focus][focus+':S'][0]*scale_s*scale_count
          Cpart=datadict[focus][focus+':C'][0]*scale_c*scale_count+datadict[focus][focus+':S:C'][0]*scale_count*scale_s*scale_c
          y=baseline+Spart+Cpart
          points.append((count,y))
      elif level=='full':
        for count in range(-2000,4000):
          scale_count=(count-scaledict[focus][0])/scaledict[focus][1]
          Spart=datadict[focus][focus][0]*scale_count+datadict[focus][focus+':S'][0]*scale_s*scale_count
          if focus!='PC3':
            Cpart=datadict[focus][focus+':C'][0]*scale_c*scale_count+datadict[focus][focus+':S:C'][0]*scale_count*scale_s*scale_c
          else:
            Cpart=0
          y=baseline+Spart+Cpart
          points.append((count,y))

      dat=graph.add_dataset(points)
      dat.line.configure(color=col,linestyle=typ,linewidth=2)
      dat.symbol.shape=0

      if C==0.02 and focus=='PC1':
        dat.legend='Species richness: '+str(S)

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


datafiles=['PC1_lmer_table.tsv','PC2_lmer_table.tsv','PC3_lmer_table.tsv']
datadict={}
scalefile='../stat_analyses/scales_for_PCA_lmer.tsv'

scaledict=read_scalefile(scalefile)
for datfile in datafiles:
  motif=datfile.split('_')[0]
  datadict[motif]=read_lmfile('../stat_analyses/'+datfile)

for form in ['talk','paper']:
  if form=='talk':
    levelar=['axis','point','single','two','full']
  else: 
    levelar=['full']
  for level in levelar:
    grace=MultiPanelGrace(colors=colors)
    grace.add_label_scheme('dummy',['','A: Axis 1','B: Axis 2','C: Axis 3','D: S5',''])
    grace.set_label_scheme('dummy')
    colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
                               scale=LINEAR_SCALE,autoscale=False)
    colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
    colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
    colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
    for axis in ['PC1','PC2','PC3']:
      graph=grace.add_graph(Panel)
      graph=format_linegraph(graph,form,axis)

      if level!='axis':
        graph=populate_graph(graph,level,datadict,scaledict,axis)

    # grace.multi(rows=3,cols=1)
    # for graph in grace.graphs:
    #   print graph.get_view()

    if form=='talk':
      grace.graphs[1].set_view(0.07, 0.74, 0.32, 0.95)
      grace.graphs[2].set_view(0.33, 0.74, 0.58, 0.95)
      grace.graphs[3].set_view(0.59, 0.74, 0.85, 0.95)
      colorbar.set_view(0.87,0.74,0.91,0.95)
    else:
      colorbar.set_view(0.46,0.27,0.51,0.95)
      grace.graphs[1].set_view(0.15, 0.75, 0.436, 0.95)
      grace.graphs[2].set_view(0.15, 0.51, 0.436, 0.71)
      grace.graphs[3].set_view(0.15, 0.27, 0.436, 0.47)   

    print level
    grace.write_file('../../manuscript/figures/extinction_order/PCA_position_lmer_summary_'+form+'_'+level+'.eps')
