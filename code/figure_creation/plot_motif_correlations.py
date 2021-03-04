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

colors=ColorBrewerScheme('Greys',reverse=True)  # The blue is very beautiful but maybe harder to see.
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


def read_file(infile):
  lmdict={}
  f=open(infile,'r')
  for line in f:
    if line.split()[0]!='"Motif"':
      name=line.split()[1][2:-1]
      raw_int=float(line.split()[2][1:-1])
      raw_beta=float(line.split()[3][1:-1])
      raw_se=float(line.split()[4][1:-1])
      raw_p=float(line.split()[5][1:-1])
      prop_int=float(line.split()[6][1:-1])
      prop_beta=float(line.split()[7][1:-1])
      prop_se=float(line.split()[8][1:-1])
      prop_p=float(line.split()[9][1:-1])
      zed_int=float(line.split()[10][1:-1])
      zed_beta=float(line.split()[11][1:-1])
      zed_se=float(line.split()[12][1:-1])
      zed_p=float(line.split()[13][1:-1])

      lmdict[name]={'raw':(raw_int,raw_beta,raw_se,raw_p),
                    'prop':(prop_int,prop_beta,prop_se,prop_p),
                    'zed':(zed_int,zed_beta,zed_se,zed_p)}
  f.close()
  return lmdict

def format_linegraph(graph,form):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=1
  graph.world.xmax=100
  graph.xaxis.tick.major=20

  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.label.configure(text='Degree',char_size=1,just=2)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  if form=='raw':
    graph.world.ymin=-50
    graph.world.ymax=1300
    graph.yaxis.tick.major=300

    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='Count of motif',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  elif form=='prop':
    graph.world.ymin=-.01
    graph.world.ymax=.55
    graph.yaxis.tick.major=.1

    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.configure(text='Proportion of role',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  else:
    graph.world.ymin=-1
    graph.world.ymax=3
    graph.yaxis.tick.major=1

    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.configure(text='Z-score of motif',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='iul',dy=0.02,dx=0.03)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def format_Tgraph(graph,form):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=1
  graph.world.xmax=6
  graph.xaxis.tick.major=1

  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.label.configure(text='Trophic level',char_size=1,just=2)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  if form=='raw':
    graph.world.ymin=0
    graph.world.ymax=400
    graph.yaxis.tick.major=100

    graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
    graph.yaxis.label.configure(text='Count of motif',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  elif form=='prop':
    graph.world.ymin=0
    graph.world.ymax=.45
    graph.yaxis.tick.major=.1

    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.configure(text='Proportion of role',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  else:
    graph.world.ymin=-.5
    graph.world.ymax=0.30000000000000000001
    graph.yaxis.tick.major=.2

    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
    graph.yaxis.label.configure(text='Z-score of motif',char_size=1,just=2)
    graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='iul',dy=0.02,dx=0.03)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def populate_graph(graph,datadict,form):
  if form=='raw':
    order=['S5','S2','S1','S4','D1','D3','D2','D7','D4','D6','D8','S3']
  else:
    order=['S5','S1','S2','S4','D1','D3','D2','D4','D7','D6','D8','S3']

  sty=1
  j=30
  for motif in order:
    points=[]
    upper=[]
    lower=[]
    for x in range(1,101):
      points.append((x,datadict[motif][form][0]+x*datadict[motif][form][1]))
      upper.append((x,datadict[motif][form][0]+x*(datadict[motif][form][1]+1.75*datadict[motif][form][2])))
      lower.append((x,datadict[motif][form][0]+x*(datadict[motif][form][1]-1.75*datadict[motif][form][2])))
    if form=='raw':
      graph.add_drawing_object(DrawText,text=motif,x=j,y=(datadict[motif][form][0]+j*datadict[motif][form][1])+10,loctype='world',just=2,char_size=.5)
    elif form=='prop':
      if motif in ['S1','S4','S5']:
        graph.add_drawing_object(DrawText,text=motif,x=j,y=(datadict[motif][form][0]+j*datadict[motif][form][1])+.005,loctype='world',just=2,char_size=.5)
      elif motif=='S2':
        graph.add_drawing_object(DrawText,text=motif,x=45,y=(datadict[motif][form][0]+45*datadict[motif][form][1])+.005,loctype='world',just=2,char_size=.5)
      else:
        graph.add_drawing_object(DrawText,text=motif,x=j,y=(datadict[motif][form][0]+j*datadict[motif][form][1])+.001,loctype='world',just=2,char_size=.5)
    else:
      if motif in ['D1','D3']:
        graph.add_drawing_object(DrawText,text=motif,x=87,y=(datadict[motif][form][0]+87*datadict[motif][form][1])-.07,loctype='world',just=2,char_size=.5)
      elif motif in ['S3','D8','D4','S4','D6']:
        graph.add_drawing_object(DrawText,text=motif,x=90,y=(datadict[motif][form][0]+90*datadict[motif][form][1])+.02,loctype='world',just=2,char_size=.5)
      elif motif in ['D7']:
        graph.add_drawing_object(DrawText,text=motif,x=92,y=(datadict[motif][form][0]+92*datadict[motif][form][1])+.02,loctype='world',just=2,char_size=.5)
      else:
        graph.add_drawing_object(DrawText,text=motif,x=90,y=(datadict[motif][form][0]+90*datadict[motif][form][1])+.05,loctype='world',just=2,char_size=.5)

    j+=5

    pointy=graph.add_dataset(points)
    pointy.symbol.shape=0
    # uppy=graph.add_dataset(upper)
    # uppy.symbol.shape=0
    if motif in ['S1','S2','S4','S5']:
      nucol='Tol4'
      col2=5
      if motif=='S1' and form=='zed':
        pointy.legend='Stable'      
    else:
      nucol='Tol7'
      col2=8
      if motif=='D3' and form=='zed':
        pointy.legend='Unstable'

    pointy.line.configure(linestyle=sty,linewidth=1,color=nucol)

    # pointy.legend=motif

  # null=graph.add_dataset([(12,0)])
  # null.symbol.configure(shape=9,size=1,fill_color=4)
  if form=='zed':
    graph.legend.configure(box_linestyle=0,char_size=.5,loc=(10,2.5),loctype='world')

  # # Add an arrow for S
  # graph.add_drawing_object(DrawLine,end=(14,0),start=(14,-1.45),arrow=1,arrow_type=1,linewidth=2,linestyle=1,loctype='world')

  return graph

def populate_Tgraph(graph,datadict,form):
  if form=='raw':
    order=['S5','S2','S1','S4','D3','D1','D2','D4','D7','D6','D8','S3']
  elif form=='prop':
    order=['S5','S1','S2','S4','D3','D1','D2','D4','D7','D6','D8','S3']
  else:
    order=['S5','S1','S2','S4','D1','D3','D2','D4','D7','D6','D8','S3']

  sty=1
  j=2
  for motif in order:
    points=[]
    upper=[]
    lower=[]
    for w in range(1,61):
      x=float(w/10)
      points.append((x,datadict[motif][form][0]+x*datadict[motif][form][1]))
      upper.append((x,datadict[motif][form][0]+x*(datadict[motif][form][1]+1.75*datadict[motif][form][2])))
      lower.append((x,datadict[motif][form][0]+x*(datadict[motif][form][1]-1.75*datadict[motif][form][2])))
    if form=='raw':
      if motif not in ['S1','S2']:
        graph.add_drawing_object(DrawText,text=motif,x=j,y=(datadict[motif][form][0]+j*datadict[motif][form][1])+5,loctype='world',just=2,char_size=.5)
      elif motif=='S1':
        graph.add_drawing_object(DrawText,text=motif,x=4,y=(datadict[motif][form][0]+4*datadict[motif][form][1])+5,loctype='world',just=2,char_size=.5)
      else:
        graph.add_drawing_object(DrawText,text=motif,x=4,y=(datadict[motif][form][0]+4*datadict[motif][form][1])-20,loctype='world',just=2,char_size=.5)
    elif form=='prop':
      if motif in ['S4','S5']:
        graph.add_drawing_object(DrawText,text=motif,x=j,y=(datadict[motif][form][0]+j*datadict[motif][form][1])-.025,loctype='world',just=2,char_size=.5)
      elif motif in ['S1','S2']:
        graph.add_drawing_object(DrawText,text=motif,x=j,y=(datadict[motif][form][0]+j*datadict[motif][form][1])+.005,loctype='world',just=2,char_size=.5)
      else:
        graph.add_drawing_object(DrawText,text=motif,x=j,y=(datadict[motif][form][0]+j*datadict[motif][form][1])+.003,loctype='world',just=2,char_size=.5)
    else:
      if motif in ['D4','S3','D3','D1','S4']:
        graph.add_drawing_object(DrawText,text=motif,x=4.5,y=(datadict[motif][form][0]+4.5*datadict[motif][form][1])+.02,loctype='world',just=2,char_size=.5)
      elif motif in ['D8','D7','D6','D2','D5']:
        graph.add_drawing_object(DrawText,text=motif,x=5,y=(datadict[motif][form][0]+5*datadict[motif][form][1])+.02,loctype='world',just=2,char_size=.5)
      elif motif in ['S1','S2','S5']:
        graph.add_drawing_object(DrawText,text=motif,x=4,y=(datadict[motif][form][0]+4*datadict[motif][form][1])-.055,loctype='world',just=2,char_size=.5)

    j+=.3

    pointy=graph.add_dataset(points)
    pointy.symbol.shape=0
    # uppy=graph.add_dataset(upper)
    # uppy.symbol.shape=0
    if motif in ['S1','S2','S4','S5']:
      nucol='Tol3'
      col2=5
      if motif=='S1' and form=='zed':
        pointy.legend='Stable'
    else:
      nucol='Tol5'
      col2=8
      if motif=='D3' and form=='zed':
        pointy.legend='Unstable'
    pointy.line.configure(linestyle=sty,linewidth=1,color=nucol)

    # pointy.legend=motif

  # null=graph.add_dataset([(12,0)])
  # null.symbol.configure(shape=9,size=1,fill_color=4)
  if form=='zed':
    graph.legend.configure(box_linestyle=0,char_size=.5,loc=(1.25,-0.2),loctype='world')

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

# # Degree vs motifs
# scalefile='motif_lm_scales.tsv'
degfile='../../data/summaries/motifs_vs_degree.tsv'
degdict=read_file(degfile)

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A','B','C'])
grace.set_label_scheme('dummy')
# colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
#                            scale=LINEAR_SCALE,autoscale=False)
# colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
# colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
# colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)

for roletype in ['raw','prop','zed']:
# for roletype in ['zed']:
  graph=grace.add_graph(Panel)
  graph=format_linegraph(graph,roletype)
  graph=populate_graph(graph,degdict,roletype)


# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=3,cols=1,vgap=.03)
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/roles/motif_vs_degree.eps')


# # TL vs motifs
TLfile='../../data/summaries/motifs_vs_STL.tsv'
TLdict=read_file(TLfile)

Tgrace=MultiPanelGrace(colors=colors)
Tgrace.add_label_scheme('dummy',['A','B','C'])
Tgrace.set_label_scheme('dummy')
# colorbar = grace.add_graph(ElLogColorBar,domain=(0.02,0.2),
#                            scale=LINEAR_SCALE,autoscale=False)
# colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
# colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
# colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)

for roletype in ['raw','prop','zed']:
# for roletype in ['zed']:
  graph=Tgrace.add_graph(Panel)
  graph=format_Tgraph(graph,roletype)
  graph=populate_Tgraph(graph,TLdict,roletype)

# graph.set_view(0.1,0.45,0.9,0.95)
Tgrace.multi(rows=3,cols=1,vgap=.03)
Tgrace.hide_redundant_labels()

Tgrace.write_file('../../manuscript/figures/roles/motif_vs_TL.eps')


grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A','B','C','D','E','F'])
grace.set_label_scheme('dummy')

for roletype in ['raw','prop','zed']:
  for pred in ['Deg','TL']:
    graph=grace.add_graph(Panel)
    if pred=='Deg':
      graph=format_linegraph(graph,roletype)
      graph=populate_graph(graph,degdict,roletype)
    elif pred=='TL':
      graph=format_Tgraph(graph,roletype)
      graph=populate_Tgraph(graph,TLdict,roletype)

grace.multi(rows=3,cols=2,vgap=.03,hgap=.05)
grace.hide_redundant_labels()
for graph in grace.graphs:
  graph.panel_label.configure(dy=.015)

grace.write_file('../../manuscript/figures/roles/motif_vs_oneD.eps')
