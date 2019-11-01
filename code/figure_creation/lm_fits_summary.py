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

def read_scalefile(scalefile):
  scales={}
  f=open(scalefile,'r')
  for line in f:
    if line.split()[0]=='"center"':
      d_center=float(line.split()[1])
      p_center=float(line.split()[2])
      c_center=float(line.split()[3])
    elif line.split()[0]=='"scale"':
      d_scale=float(line.split()[1])
      p_scale=float(line.split()[2])
      c_scale=float(line.split()[3])  
  f.close()

  scales['Dissim']=(d_center,d_scale)
  scales['Path']=(p_center,p_scale)
  scales['C']=(c_center,c_scale)  
  return scales

def read_lmfile(infile):
  lmdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"Estimate"':
      beta=float(line.split()[1])      
      p=float(line.split()[-1])
      if line.split()[0]=='"scale(Dissim)"':
        lmdict['Dissim']=(beta,p)
      elif line.split()[0]=='"(Intercept)"':
        lmdict['Intercept']=(beta,p)
      elif line.split()[0]=='"scale(Path)"':
        lmdict['Path']=(beta,p)
      elif line.split()[0]=='"scale(C)"':
        lmdict['C']=(beta,p)
      elif line.split()[0]=='"scale(Dissim):scale(Path)"':
        lmdict['Dissim:Path']=(beta,p)
      elif line.split()[0]=='"scale(Dissim):scale(C)"':
        lmdict['Dissim:C']=(beta,p)
      elif line.split()[0]=='"scale(Path):scale(C)"':
        lmdict['Path:C']=(beta,p)
      elif line.split()[0]=='"scale(Dissim):scale(Path):scale(C)"':
        lmdict['3way']=(beta,p)
      else:
        print line.split()
  f.close()

  return lmdict

def format_linegraph(graph,form,focus):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  graph.world.xmin=0
  graph.world.ymax=50
  graph.yaxis.tick.major=10

  if focus=='dissim':
    graph.world.xmax=1
    graph.xaxis.tick.major=0.2
    graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.xaxis.label.configure(text='Dissimiliarity',char_size=1,just=2)
    graph.yaxis.label.configure(text='Time to extinction',char_size=1,just=2)
    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  elif focus=='path':
    graph.world.xmax=8
    graph.xaxis.tick.major=2
    graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.xaxis.label.configure(text='Path length',char_size=1,just=2)
    if form=='talk': 
      graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=0)
    if form=='paper':
      graph.yaxis.label.configure(text='Time to extinction',char_size=1,just=2)
      graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)


  return graph


def populate_graph(graph,level,datadict,scaledict,focus):
  # if level=='point':
  #   c=datadict[50].keys()[0]
  #   inter=datadict[50][c][0]
  #   beta=datadict[50][c][1]
  #   points=[]
  #   for dissim in range(0,21):
  #     x=float(dissim)/float(20)
  #     y=inter+beta*x
  #     points.append((x,y))
  #   dat=graph.add_dataset(points)
  #   dat.line.configure(color=colorbar.z2color(c),linestyle=1,linewidth=1)
  #   dat.symbol.shape=0
  # else:
  if focus=='dissim':
    for S in datadict:
      for c in [0.02,0.06,0.1,0.14,0.18,0.2]:
        col=colorbar.z2color(c)
        scale_c=(c-scaledict[str(S)]['C'][0])/scaledict[str(S)]['C'][1]
        baseline=datadict[S]['Intercept'][0]+datadict[S]['C'][0]*scale_c
        for path in [1,2,4,8]:
          scale_p=(path-scaledict[str(S)]['Path'][0])/scaledict[str(S)]['Path'][1]
          pluspath=baseline+datadict[S]['Path'][0]*scale_p+datadict[S]['Path:C'][0]*scale_p*scale_c
          points=[]
          for dissim in range(0,21):
            scale_d=(float(dissim)/20-scaledict[str(S)]['Dissim'][0])/scaledict[str(S)]['Dissim'][1]
            d1=datadict[S]['Dissim'][0]*scale_d+datadict[S]['Dissim:Path'][0]*scale_d*scale_p
            d2=datadict[S]['Dissim:C'][0]*scale_d*scale_c+datadict[S]['3way'][0]*scale_d*scale_p*scale_c
            x=float(dissim)/float(20)
            y=pluspath+d1+d2
            points.append((x,y))
        dat=graph.add_dataset(points)
        dat.line.configure(color=col,linestyle=1,linewidth=1)
        dat.symbol.shape=0
  elif focus=='path':
    for S in datadict:
      for c in [0.02,0.06,0.1,0.16,0.2]:
        col=colorbar.z2color(c)
        scale_c=(c-scaledict[str(S)]['C'][0])/scaledict[str(S)]['C'][1]
        baseline=datadict[S]['Intercept'][0]+datadict[S]['C'][0]*scale_c
        for dissim in range(0,21):
          scale_d=(float(dissim)/20-scaledict[str(S)]['Dissim'][0])/scaledict[str(S)]['Dissim'][1]
          plusdis=baseline+datadict[S]['Dissim'][0]*scale_d+datadict[S]['Dissim:C'][0]*scale_d*scale_c
          points=[]
          for path in range(0,8):
            scale_p=(path-scaledict[str(S)]['Path'][0])/scaledict[str(S)]['Path'][1]
            p1=datadict[S]['Path'][0]*scale_p+datadict[S]['Path:C'][0]*scale_p*scale_c
            p2=datadict[S]['Dissim:Path'][0]*scale_d*scale_p+datadict[S]['3way'][0]*scale_d*scale_p*scale_c
            x=path
            y=plusdis+p1+p2
            points.append((x,y))
        dat=graph.add_dataset(points)
        dat.line.configure(color=col,linestyle=1,linewidth=1)
        dat.symbol.shape=0

  # graph.add_drawing_object(DrawText,text='p<0.001',x=50,y=0.5,loctype='world',char_size=0.75,just=2)

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

datadir='../../data/summaries/rolesim_extorder/'
datadict={}
scaledict={}
for s in ['50','60','70','80','90','100']:
  datadict[int(s)]={}
  # options=os.listdir(datadir+s)
  lmfile='role_exttime_lm_with_C.tsv'
  datadict[int(s)]=read_lmfile(datadir+s+'/'+lmfile)
  scalefile=datadir+'scales_'+s+'.tsv'
  scaledict[s]=read_scalefile(scalefile)

for form in ['talk','paper']:
  for level in ['axis','point','full']:
    grace=MultiPanelGrace(colors=colors)
    grace.add_label_scheme('dummy',['','','','',''])
    grace.set_label_scheme('dummy')
    colorbar = grace.add_graph(ElLogColorBar,domain=(0,0.2),
                               scale=LINEAR_SCALE,autoscale=False)
    colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
    colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
    colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
    dlinegraph=grace.add_graph(Panel)
    plinegraph=grace.add_graph(Panel)

    dlinegraph=format_linegraph(dlinegraph,form,'dissim')
    plinegraph=format_linegraph(plinegraph,form,'path')

    if level in ['point','full']:
      dlinegraph=populate_graph(dlinegraph,level,datadict,scaledict,'dissim')
      plinegraph=populate_graph(plinegraph,level,datadict,scaledict,'path')

    if form=='talk':
      dlinegraph.set_view(0.1,0.2,0.45,0.5)
      plinegraph.set_view(0.5,0.2,0.9,0.5)
      colorbar.set_view(0.95,0.2,1.0,0.5)
    elif form=='paper':
      dlinegraph.set_view(0.1,0.55,0.55,0.9)
      plinegraph.set_view(0.1,0.1,0.55,0.45)
      colorbar.set_view(0.6,0.1,0.65,0.9)

    print level
    grace.write_file('../../manuscript/figures/extinction_order/dissimilarity_fits_summary_'+form+'_'+level+'.eps')
