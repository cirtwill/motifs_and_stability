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

def format_graph(graph,normtype,coefs):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  graph.world.ymax=19
  graph.yaxis.tick.major=0.5
  graph.xaxis.tick.major=1
  if normtype=='Network':
    graph.world.xmin=-2
    graph.world.xmax=2
  else:
    graph.world.xmin=-8
    graph.world.xmax=2.5

  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.yaxis.label.configure(text='Predictor',char_size=1,just=2)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.set_spec_ticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],[],
    tick_labels=['S:C','C','S','STL','Degree','D8','D7','D6','D5','D4','D3','D2','D1',
        'S5','S4','S3','S2','S1'])

  graph.xaxis.ticklabel.configure(format='decimal',prec=1,char_size=.75)
  graph.xaxis.label.configure(text='Sum of coefficients',char_size=1,just=2)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='ouc',dx=0,dy=.01)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))
  # graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)/
  return graph

def populate_graph(graph,normtype,coefs):
  if normtype=='Degree':
    predlist=['S:C','C','S','STL','Degree','cD8','cD7','cD6','cD5','cD4','cD3','cD2','cD1',
        'cS5','cS4','cS3','cS2','cS1']
  else:
    predlist=['S:C','C','S','STL','Degree','zD8','zD7','zD6','zD5','zD4','zD3','zD2','zD1',
        'zS5','zS4','zS3','zS2','zS1']

  # Add grey backgrounds for stable motifs
  upper=graph.add_dataset([(-15,18.5),(15,18.5)])
  upper.line.linestyle=0
  upper.fill.configure(color=70,type=2)
  S2=graph.add_dataset([(-15,16.5),(15,16.5)])  
  S2.fill.configure(color=0,type=2)
  S2.line.linestyle=0
  S4=graph.add_dataset([(-15,15.5),(15,15.5)])
  S4.fill.configure(color=70,type=2)
  S4.line.linestyle=0
  S5=graph.add_dataset([(-15,13.5),(15,13.5)])
  S5.fill.configure(color=0,type=2)
  S5.line.linestyle=0

  y=1
  for pred in predlist:
    bar=graph.add_dataset([(0,y),(coefs[pred],y)],)
    bar.symbol.shape=0
    bar.line.configure(linewidth=6,linestyle=1,color=3)
    y+=1


  # Add a line at 0
  graph.add_drawing_object(DrawLine,end=(0,0),start=(0,19),linewidth=1,linestyle=1,loctype='world')


  return graph

        # color = colorbar.z2color(pdf)
        # # you can change the opacity percentage of a single color, as well
        # # color.change_opacity(60)
        # graph.add_dataset([(x0,y0), (x1,y1)], SolidRectangle, color)

def read_coeffile(infile):
  subdict={}
  f=open(infile,'r')

  for line in f:
    if line.split('\t')[0]!='"metapreds$Persistence.1 comps"':
      pred=line.split()[0]
      if ')' in pred:
        pred=pred.split(')')[1]
      if '$' in pred:
        if ':' in pred:
          pred=pred.split('$')[1].split(':')[0]+':'+pred.split('$')[-1]
        else:
          pred=pred.split('$')[1]
      if '"' in pred:
        pred=''.join(pred.split('"'))
      if pred=='S:metapreds':
        print line.split()[0]
      cos=[]
      for val in line.split()[1:]:
        cos.append(float(val))        
      total=sum(cos)
      subdict[pred]=total

  f.close()

  return subdict

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

degfile='../../data/PLS_regression/degnorm_coefficients.tsv'
netfile='../../data/PLS_regression/netnorm_coefficients.tsv'

coefs={'Degree':{},'Network':{}}
coefs['Degree']=read_coeffile(degfile)
coefs['Network']=read_coeffile(netfile)

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A) Degree normalization','B) Network normalization','S2: omnivory','S4: direct competition','S5: apparent competition',''])
grace.set_label_scheme('dummy')
for normtype in ['Degree','Network']:
  graph=grace.add_graph(Panel)
  graph=format_graph(graph,normtype,coefs)
  graph=populate_graph(graph,normtype,coefs[normtype])

# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=2,cols=1)
grace.hide_redundant_labels()

grace.write_file('../../manuscript/figures/PLS/total_coefficients.eps')