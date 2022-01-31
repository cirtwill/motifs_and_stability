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

def read_means(meanfile):
  motif_means={}
  f=open(meanfile,'r')
  for line in f:
    if len(line.split())>1:
      name=line.split()[0][1:-1]
      if len(name.split('.'))==3:
        motif=str(int(name.split('.')[0][1:]))
        preds=int(name.split('.')[1])
        preys=int(name.split('.')[2])
        name=(motif,preds,preys)
      else:
        name='Other'
      mean=float(line.split()[1])
      motif_means[name]=mean

  return motif_means

def read_file(infile):
  lmdict={}
  f=open(infile,'r')
  for line in f:

    if line.split()[0]!='"Estimate"':
      name=line.split()[0][1:-1]
      if 'Intercept' in name:
        name='intercept'
      elif len(name.split('.'))==3:
        motif=str(int(name.split('.')[0][1:]))
        preds=int(name.split('.')[1])
        preys=int(name.split('.')[2])
        name=(motif,preds,preys)
      else:
        name='Other'
      effect=float(line.split()[1])
      pval=float(line.split()[-1])

      lmdict[name]=((effect,pval))

  return lmdict

def format_linegraph(graph,normtype,simple,motif):
  motdict={'6':'App. Comp.','36':'Dir. Comp.','38':'Omnivory','12':'Chain'}
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.xmin=0
  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.5,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  if normtype=='count':
    graph.world.xmax=1000
    graph.xaxis.tick.major=250
    graph.xaxis.label.configure(text='Count',char_size=1,just=2)
  elif normtype=='freq':
    graph.world.xmax=1
    graph.xaxis.tick.major=.2
    # graph.xaxis.label.configure(text='Frequency',char_size=1,just=2)
    graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  else:
    graph.world.xmin=-3.65
    graph.world.xmax=9.56
    graph.xaxis.tick.major=3
    graph.xaxis.label.configure(text='Z-score',char_size=1,just=2)
    graph.yaxis.label.configure(text=motdict[motif],char_size=1,just=2,place='opposite')

  if simple=='Deg':
    graph.world.ymin=0
    graph.world.ymax=125
    graph.yaxis.tick.major=25
  else:
    graph.world.ymin=0
    graph.world.ymax=3
    graph.yaxis.tick.major=1
  graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.5,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='iul',dy=0.02,dx=0.03)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))

  return graph

def populate_graph(graph,datadict,normtype,simple,motif,motif_means):
  if normtype=='count':
    xs=range(1,1101,100)
  elif normtype=='freq':
    bix=range(0,101,5)
    xs=[float(x)/100 for x in bix]
  else:
    bix=range(-9,21)
    xs=[float(x)/2 for x in bix]
  positions=[d for d in datadict.keys() if d[0]==motif]

  for p in sorted(positions):
    dats=[]
    for x in xs:
      motpart=datadict['intercept'][0]+x*datadict[p][0]
      if normtype!='freq':
        otherpart=sum([motif_means[d]*datadict[d][0] for d in datadict.keys() if d not in [p,'intercept']])
      else:
        otherpart=sum([(1-x)*float(motif_means[d]/sum(motif_means.values()))*datadict[d][0]for d in datadict.keys() if d not in [p,'intercept']])
      dats.append((x,motpart+otherpart))

    pointy=graph.add_dataset(dats)
    pointy.symbol.shape=0
    if p[1]==0: # Top
      pointy.line.linestyle=1
      j=12
      if normtype=='Z':
        pointy.legend='Top'
    elif p[2]==0: # Bottom
      pointy.line.linestyle=2
      j=21
      if normtype=='Z':
        pointy.legend='Bottom'
    else:
      pointy.line.linestyle=3
      j=14
      if normtype=='Z':
        pointy.legend='Middle'

    pointy.line.configure(linewidth=3,color=j)

  # graph.add_drawing_object(DrawText,text='STL',x=90,y=16,loctype='world',just=0,char_size=1)
  if normtype=='Z':
    if simple=='Deg':
      graph.legend.configure(box_linestyle=0,char_size=.5,loc=(0,120),loctype='world')
    else:
      graph.legend.configure(box_linestyle=0,char_size=.5,loc=(-3,1.5),loctype='world')

  # # # Add an arrow for S
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
simple='TL'
for motif in ['6','36','38','12']:
  for normtype in ['count','freq','Z']:
    datfile='../../data/summaries/positions_'+normtype+'_'+simple+'.tsv'
    datdict=read_file(datfile)
    motif_means=read_means('../../data/summaries/mean_positions_'+normtype+'.tsv')
    # Need to update populate_graph
    graph=grace.add_graph(Panel)
    graph=format_linegraph(graph,normtype,simple,motif)
    graph=populate_graph(graph,datdict,normtype,simple,motif,motif_means)


# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=4,cols=3,vgap=.03,hgap=.05)
if simple=='Deg':
  grace.set_col_yaxislabel(col=0,rowspan=(None,None),label="Degree",char_size=1.5,just=2,place='normal')
else:
  grace.set_col_yaxislabel(col=0,rowspan=(None,None),label="Trophic level",char_size=1.5,just=2,place='normal')
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/positions_vs_'+simple+'.eps')


grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['App. Comp.','Dir. Comp.','Omnivory','Chain','E','F'])
grace.set_label_scheme('dummy')
simple='Deg'
normtype='freq'
for motif in ['6','36','38','12']:
  datfile='../../data/summaries/positions_'+normtype+'_'+simple+'.tsv'
  datdict=read_file(datfile)
  motif_means=read_means('../../data/summaries/mean_positions_'+normtype+'.tsv')
  # Need to update populate_graph
  graph=grace.add_graph(Panel)
  graph=format_linegraph(graph,normtype,simple,motif)
  graph=populate_graph(graph,datdict,normtype,simple,motif,motif_means)

# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=2,cols=2,vgap=.03,hgap=.05)
grace.set_col_yaxislabel(col=0,rowspan=(None,None),label="Degree",char_size=1.5,just=2,place='normal')
grace.set_row_xaxislabel(row=1,colspan=(None,None),label="Frequency",char_size=1.5,just=2,place='normal')
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/positions_vs_Deg_freq.eps')


grace=MultiPanelGrace(colors=colors)
for motif in ['6','36','38','12']:
  for normtype in ['count','Z']:
    datfile='../../data/summaries/positions_'+normtype+'_'+simple+'.tsv'
    datdict=read_file(datfile)
    motif_means=read_means('../../data/summaries/mean_positions_'+normtype+'.tsv')
    # Need to update populate_graph
    graph=grace.add_graph(Panel)
    graph=format_linegraph(graph,normtype,simple,motif)
    graph=populate_graph(graph,datdict,normtype,simple,motif,motif_means)


# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=4,cols=2,vgap=.03,hgap=.05)
grace.set_col_yaxislabel(col=0,rowspan=(None,None),label="Degree",char_size=1.5,just=2,place='normal')
# grace.set_row_xaxislabel(row=0,colspan=(None,None),label="Count",char_size=1,just=2,place='normal',perpendicular_offset=0.05)
# grace.set_row_xaxislabel(row=1,colspan=(None,None),label="Z-score",char_size=1,just=2,place='normal',perpendicular_offset=0.05)
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/positions_vs_Deg_countZ.eps')




