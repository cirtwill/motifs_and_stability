import sys
import math
import random
from decimal import *
import numpy as np

#Pygrace libraries
from PyGrace.grace import Grace
from PyGrace.colors import RandomColorScheme, MarkovChainColorScheme, ColorBrewerScheme
from PyGrace.dataset import SYMBOLS
from PyGrace.Extensions.panel import NetworkPanel,Panel,MultiPanelGrace
from PyGrace.drawing_objects import DrawText, DrawLine, DrawBox
from PyGrace.axis import LINEAR_SCALE, LOGARITHMIC_SCALE
from PyGrace.Extensions.network import Network

from PyGrace.Extensions.distribution import CDFGraph, PDFGraph
from PyGrace.Extensions.latex_string import LatexString, CONVERT
from PyGrace.Extensions.colorbar import SolidRectangle, ColorBar
from PyGrace.Styles.el import ElGraph, ElLinColorBar, ElLogColorBar

colors=ColorBrewerScheme('Spectral')  # The blue is very beautiful but maybe harder to see.
colors.add_color(120,120,120,'grey')

def read_file(datafile):
  datadict={}
  f=open(datafile,'r')
  for line in f:
    if len(line.split())==4:
      position=line.split()[0][1:-1]
      ax1=float(line.split()[1])
      ax2=float(line.split()[2])
      ax3=float(line.split()[3])
      datadict[position]=((ax1,ax2,ax3))
  f.close()
  return datadict

# Positions and motifs
def format_netgraph(graph,form):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  if form=='paper':
    graph.world.ymax=.75
    graph.world.ymin=-12.75
    graph.world.xmin=-.5
    graph.world.xmax=1
  else:
    graph.world.xmin=-0.75
    graph.world.xmax=12.75
    graph.world.ymin=-.5
    graph.world.ymax=1
  graph.xaxis.onoff='on'
  graph.yaxis.onoff='on'
  graph.xaxis.ticklabel.configure(char_size=0,format='decimal',prec=1)
  graph.yaxis.ticklabel.configure(char_size=0,format='decimal',prec=1)
  graph.xaxis.tick.onoff='off'
  graph.yaxis.tick.onoff='off'

  graph.panel_label.configure(placement='iul',char_size=.75,dx=.01,dy=.01)
  return graph

def add_positions(graph,level):
  y=0.5
  blackset={'1a':(0,y),'2a':(1,0),'3a':(2,y),'4b':(3,0),'5a':(4,y),'6a':(5,0),
            '7a':(6,y),'8a':(7,0),'10a':(8.6,y),'11a':(10,y),'12a':(11,0)}
  whiteset={'1b':(-.4,0),'1c':(.4,0),'2b':(1.4,y),'3b':(1.6,0),'4a':(2.6,y),'4c':(3.4,y),'5b':(3.6,0),'6b':(4.6,y),'6c':(5.4,y),
            '7b':(5.6,0),'8b':(6.6,y),'8c':(7.4,y),'10b':(9,0),
            '11b':(9.6,0),'11c':(10.4,0),'12b':(10.6,y)}
  greyset={'2c':(.6,y),'3c':(2.4,0),'5c':(4.4,0),
           '7c':(6.4,0),'9a':(8,y),'9b':(7.6,0),'9c':(8.4,0),'10c':(9.4,y),
           '12c':(11.4,y),'13a':(12,y),'13b':(11.6,0),'13c':(12.4,0)
  }
  # Links go prey, pred
  links=[('1b','1a'),('1c','1a'),('2a','2b'),('2b','2c'),('3a','3c'),('3c','3a'),('3b','3a'),
    ('4b','4a'),('4b','4c'),('5b','5a'),('5c','5a'),('5c','5b'),('6b','6c'),('6c','6b'),('6a','6b'),('6a','6c'),
    ('7c','7a'),('7b','7c'),('7c','7b'),('8a','8b'),('8a','8c'),('8b','8a'),('8c','8a'),('9b','9a'),('9c','9b'),('9a','9c'),
    ('10b','10a'),('10b','10c'),('10c','10b'),('10a','10c'),
    ('11b','11a'),('11c','11a'),('11b','11c'),('11c','11b'),
    ('12a','12c'),('12c','12a'),('12b','12c'),('12c','12b'),('12a','12b'),
    ('13a','13b'),('13b','13a'),('13a','13c'),('13c','13a'),('13b','13c'),('13c','13b')
    ]

  graph=add_numbers(graph)

  graph.add_node_set(whiteset,color=0,size=.75)
  graph.add_node_set(blackset,color=1,size=.75) 
  graph.add_node_set(greyset,color='grey',size=.75)

  if level=='Ax1':
    rednodes={'1a':(0,y),'1b':(-.4,0),'1c':(.4,0),'5a':(4,y)}
  elif level=='Ax2':
    rednodes={'1b':(-.4,0),'1c':(.4,0),'2c':(.6,y)}
  elif level=='Ax3':
    rednodes={'1a':(0,y),'2b':(1.4,y),'4a':(2.6,y),'4c':(3.4,y),'5b':(3.6,0)}
  elif level=='allred':
    rednodes={'1a':(0,y),'1b':(-.4,0),'1c':(.4,0),'2b':(1.4,y),'2c':(.6,y),'4a':(2.6,y),'4c':(3.4,y),'5a':(4,y),'5b':(3.6,0)}
  else:
    rednodes={}
  graph.add_node_set(rednodes,color=3,size=.75)

  # Motifs in first rank are: 6, 12, 14, 36, 38, 46
  # Motifs in second rank are: 74, 78, 98, 102, 108, 110, 238
  # First rank: (6, 0, 2) (6, 1, 0) (12, 0, 1) (12, 1, 0) (12, 1, 1) (14, 1, 0) (14, 1, 1) (14, 1, 2) (36, 0, 1) (36, 2, 0) (38, 0, 2) (38, 1, 1) (38, 2, 0) (46, 1, 2) (46, 2, 0) 
  # Second rank: (74, 0, 1) (74, 1, 1) (74, 2, 1) (78, 1, 1) (78, 2, 2) (98, 1, 1) (102, 1, 1) (102, 1, 2) (102, 2, 1) (108, 0, 2) (108, 2, 1) (110, 1, 2) (110, 2, 1) (110, 2, 2) (238, 2, 2)

  lynx=graph.add_directed_link_set(links,curvature=0,avoid_crossing_nodes=False,color=1)

  return graph

def add_paper_positions(graph):
  x=0.5
  y=0.5
  blackset={'1a':(x,0),'2a':(0,-1),'3a':(x,-2),'4b':(0,-3),'5a':(x,-4),'6a':(0,-5),
            '7a':(x,-6),'8a':(0,-7),'10a':(x,-8.6),'11a':(x,-10),'12a':(0,-11)}
  whiteset={'1b':(0,-.4),'1c':(0,.4),'2b':(x,-1.4),'3b':(0,-1.6),'4a':(x,-2.6),'4c':(x,-3.4),'5b':(0,-3.6),'6b':(x,-4.6),'6c':(x,-5.4),
            '7b':(0,-5.6),'8b':(x,-6.6),'8c':(x,-7.4),'10b':(0,-9),
            '11b':(0,-9.6),'11c':(0,-10.4),'12b':(x,-10.6)}
  greyset={'2c':(x,-.6),'3c':(0,-2.4),'5c':(0,-4.4),
           '7c':(0,-6.4),'9a':(x,-8),'9b':(0,-7.6),'9c':(0,-8.4),'10c':(x,-9.4),
           '12c':(x,-11.4),'13a':(x,-12),'13b':(0,-11.6),'13c':(0,-12.4)
  }
  # Links go prey, pred
  links=[('1b','1a'),('1c','1a'),('2a','2b'),('2b','2c'),('3a','3c'),('3c','3a'),('3b','3a'),
    ('4b','4a'),('4b','4c'),('5b','5a'),('5c','5a'),('5c','5b'),('6b','6c'),('6c','6b'),('6a','6b'),('6a','6c'),
    ('7c','7a'),('7b','7c'),('7c','7b'),('8a','8b'),('8a','8c'),('8b','8a'),('8c','8a'),('9b','9a'),('9c','9b'),('9a','9c'),
    ('10b','10a'),('10b','10c'),('10c','10b'),('10a','10c'),
    ('11b','11a'),('11c','11a'),('11b','11c'),('11c','11b'),
    ('12a','12c'),('12c','12a'),('12b','12c'),('12c','12b'),('12a','12b'),
    ('13a','13b'),('13b','13a'),('13a','13c'),('13c','13a'),('13b','13c'),('13c','13b')
    ]

  graph=add_paper_numbers(graph)

  graph.add_node_set(blackset,color=1,size=.75) 
  graph.add_node_set(whiteset,color=0,size=.75)
  graph.add_node_set(greyset,color='grey',size=.75)

  lynx=graph.add_directed_link_set(links,curvature=0,avoid_crossing_nodes=False,color=1)

  return graph

def add_numbers(graph):
  # 1-6
  graph.add_drawing_object(DrawText,text='1',x=0,y=.65,char_size=.65,just=2,loctype='world')
  for x in [-.4,.4]:
    graph.add_drawing_object(DrawText,text='2',x=x,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='3',x=.6,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='4',x=1,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='5',x=1.4,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='6',x=1.6,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='7',x=2.4,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='8',x=2,y=.65,char_size=.65,just=2,loctype='world')
  for x in [2.6,3.4]:
    graph.add_drawing_object(DrawText,text='9',x=x,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='10',x=3,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='11',x=4,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='12',x=3.6,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='13',x=4.4,y=-.3,char_size=.65,just=2,loctype='world')
  for x in [4.6,5.4]:
    graph.add_drawing_object(DrawText,text='14',x=x,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='15',x=5,y=-.3,char_size=.65,just=2,loctype='world')
  # 7-13
  graph.add_drawing_object(DrawText,text='16',x=6,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='17',x=5.6,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='18',x=6.4,y=-.3,char_size=.65,just=2,loctype='world')
  for x in [6.6,7.4]:
    graph.add_drawing_object(DrawText,text='19',x=x,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='20',x=7,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='21',x=8,y=.65,char_size=.65,just=2,loctype='world')
  for x in [7.6,8.4]:
    graph.add_drawing_object(DrawText,text='21',x=x,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='22',x=8.6,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='23',x=9.4,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='24',x=9,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='25',x=10,y=.65,char_size=.65,just=2,loctype='world')
  for x in [9.6,10.4]:
    graph.add_drawing_object(DrawText,text='26',x=x,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='27',x=10.6,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='29',x=11.4,y=.65,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='28',x=11,y=-.3,char_size=.65,just=2,loctype='world')
  for x in [12.4,11.6]:
    graph.add_drawing_object(DrawText,text='30',x=x,y=-.3,char_size=.65,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='30',x=12,y=.65,char_size=.65,just=2,loctype='world')

  return graph

def add_paper_numbers(graph):
  # 1-6
  graph.add_drawing_object(DrawText,text='1',y=-0-.09,x=.65,char_size=.65,just=0,loctype='world')
  for x in [-.4,.4]:
    graph.add_drawing_object(DrawText,text='2',y=-x-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='3',y=-.6-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='4',y=-1-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='5',y=-1.4-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='6',y=-1.6-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='7',y=-2.4-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='8',y=-2-.09,x=.65,char_size=.65,just=0,loctype='world')
  for x in [2.6,3.4]:
    graph.add_drawing_object(DrawText,text='9',y=-x-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='10',y=-3-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='11',y=-4-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='12',y=-3.6-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='13',y=-4.4-.09,x=-.125,char_size=.65,just=1,loctype='world')
  for x in [4.6,5.4]:
    graph.add_drawing_object(DrawText,text='14',y=-x-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='15',y=-5-.09,x=-.125,char_size=.65,just=1,loctype='world')
  # 7-13
  graph.add_drawing_object(DrawText,text='16',y=-6-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='17',y=-5.6-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='18',y=-6.4-.09,x=-.125,char_size=.65,just=1,loctype='world')
  for x in [6.6,7.4]:
    graph.add_drawing_object(DrawText,text='19',y=-x-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='20',y=-7-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='21',y=-8-.09,x=.65,char_size=.65,just=0,loctype='world')
  for x in [7.6,8.4]:
    graph.add_drawing_object(DrawText,text='21',y=-x-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='22',y=-8.6-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='23',y=-9.4-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='24',y=-9-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='25',y=-10-.09,x=.65,char_size=.65,just=0,loctype='world')
  for x in [9.6,10.4]:
    graph.add_drawing_object(DrawText,text='26',y=-x-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='27',y=-10.6-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='29',y=-11.4-.09,x=.65,char_size=.65,just=0,loctype='world')
  graph.add_drawing_object(DrawText,text='28',y=-11-.09,x=-.125,char_size=.65,just=1,loctype='world')
  for x in [12.4,11.6]:
    graph.add_drawing_object(DrawText,text='30',y=-x-.09,x=-.125,char_size=.65,just=1,loctype='world')
  graph.add_drawing_object(DrawText,text='30',y=-12-.09,x=.65,char_size=.65,just=0,loctype='world')

  return graph

# PCA data
def format_graph(graph,graphtype,form):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.world.xmin=-.5
  graph.world.xmax=1
  graph.world.ymin=-.2
  graph.world.ymax=.2

  if form=='paper':
    if '1' in graphtype:
      graph.xaxis.label.configure(text='PCA axis 1 (31.5%)',char_size=.75,just=2)
      graph.xaxis.ticklabel.configure(char_size=0,format='decimal',prec=1)
    else:
      graph.xaxis.label.configure(text='PCA axis 2 (14.6%)',char_size=.75,just=2)
      graph.xaxis.ticklabel.configure(char_size=.5,format='decimal',prec=1)

    if graphtype[1]=='2':
      graph.yaxis.label.configure(text='PCA axis 2 (14.6%)',char_size=.75,just=2)
    else:
      graph.yaxis.label.configure(text='PCA axis 3 (12.1%)',char_size=.75,just=2)
    graph.yaxis.ticklabel.configure(char_size=.5,format='decimal',prec=1)
  else:
    if '1' in graphtype:
      graph.xaxis.label.configure(text='PCA axis 1 (31.5%)',char_size=.75,just=2)
    else:
      graph.xaxis.label.configure(text='PCA axis 2 (14.6%)',char_size=.75,just=2)

    if graphtype[1]=='2':
      graph.yaxis.label.configure(text='PCA axis 2 (14.6%)',char_size=.75,just=2)
    else:
      graph.yaxis.label.configure(text='PCA axis 3 (12.1%)',char_size=.75,just=2)
    graph.xaxis.ticklabel.configure(char_size=.5,format='decimal',prec=1)
    if graphtype=='12':
      graph.yaxis.ticklabel.configure(char_size=.5,format='decimal',prec=1)
    else:
      graph.yaxis.ticklabel.char_size=0
  graph.xaxis.tick.configure(major=.5,onoff='on',minor_ticks=1,major_size=.5,place='normal',minor_size=.3,major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.configure(major=.1,onoff='on',minor_ticks=1,major_size=.5,place='both',minor_size=.3,major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)

  return graph

def populate_graph(graph,graphtype):

  points=[]
  redpoints=[]
  for position in datadict:
    if position in ['Position_1','Position_2','Position_3','Position_5','Position_9','Position_11','Position_12']:
      if graphtype=='12':
        redpoints.append((datadict[position][0],datadict[position][1]))
      elif graphtype=='13':
        redpoints.append((datadict[position][0],datadict[position][2]))
      elif graphtype=='23':
        redpoints.append((datadict[position][1],datadict[position][2]))
    else:
      if graphtype=='12':
        points.append((datadict[position][0],datadict[position][1]))
      elif graphtype=='13':
        points.append((datadict[position][0],datadict[position][2]))
      elif graphtype=='23':
        points.append((datadict[position][1],datadict[position][2]))

  dots=graph.add_dataset(points,type='xy')
  dots.symbol.configure(fill_color=0,fill_pattern=0,color=1,size=.5)
  dots.line.linestyle=0

  reddots=graph.add_dataset(redpoints,type='xy')
  reddots.symbol.configure(fill_color=3,fill_pattern=1,color=1,size=.5)
  reddots.line.linestyle=0

  if graphtype=='12':
    for position in ['Position_1','Position_2','Position_3','Position_9','Position_11']:
      x=datadict[position][0]
      y=datadict[position][1]+0.0175
      graph.add_drawing_object(DrawText,text=position.split('_')[1],x=x,y=y,loctype='world',char_size=.5,just=2)

    for position in ['Position_5','Position_12']:
      x=datadict[position][0]
      y=datadict[position][1]-0.03
      graph.add_drawing_object(DrawText,text=position.split('_')[1],x=x,y=y,loctype='world',char_size=.5,just=2)

  elif graphtype=='13':
    for position in ['Position_1','Position_2','Position_3','Position_5','Position_9','Position_11']:
      x=datadict[position][0]
      y=datadict[position][2]+0.0175
      graph.add_drawing_object(DrawText,text=position.split('_')[1],x=x,y=y,loctype='world',char_size=.5,just=2)

    for position in ['Position_12']:
      x=datadict[position][0]
      y=datadict[position][2]-0.03
      graph.add_drawing_object(DrawText,text=position.split('_')[1],x=x,y=y,loctype='world',char_size=.5,just=2)
  else:
    for position in ['Position_1','Position_3','Position_12']:
      x=datadict[position][1]+0.04
      y=datadict[position][2]-0.01
      graph.add_drawing_object(DrawText,text=position.split('_')[1],x=x,y=y,loctype='world',char_size=.5,just=0)
    for position in ['Position_9']:
      x=datadict[position][1]
      y=datadict[position][2]+0.0175        
      graph.add_drawing_object(DrawText,text=position.split('_')[1],x=x,y=y,loctype='world',char_size=.5)
    for position in ['Position_2','Position_5','Position_11']:
      x=datadict[position][1]-0.05
      y=datadict[position][2]-0.01
      graph.add_drawing_object(DrawText,text=position.split('_')[1],x=x,y=y,loctype='world',char_size=.5,just=1)


  return graph

###########################################################################
## Main
###########################################################################



datafile='mean_loadings_PCAaxes.tsv'
datadict=read_file(datafile)


for form in ['talk','paper']:
  if form=='talk':
    levels=['axis','points','Ax1','Ax2','Ax3','allred']
  else:
    levels=['points']
  for level in levels:
    print form, level
    grace=MultiPanelGrace(colors=colors)
    grace.add_label_scheme('dummy',['A','B','C','D'])
    grace.set_label_scheme('dummy')
    graph12=grace.add_graph(Panel)
    graph13=grace.add_graph(Panel)
    graph23=grace.add_graph(Panel)
    positions=grace.add_graph(NetworkPanel)

    graph12=format_graph(graph12,'12',form)
    graph13=format_graph(graph13,'13',form)
    graph14=format_graph(graph23,'23',form)
    positions=format_netgraph(positions,form)

    if level!='axis':
      graph12=populate_graph(graph12,'12')
      graph13=populate_graph(graph13,'13')
      graph23=populate_graph(graph23,'23')

    if form=='talk':
      # grace.multi(rows=1,cols=3,vgap=.08,hgap=.04)
      graph12.set_view(0.15, 0.7, 0.4, 0.95)
      graph13.set_view(0.45, 0.7, 0.7, 0.95)
      graph23.set_view(0.75, 0.7, 1.0, 0.95)
      positions.set_view(0.15,0.475,1.0,0.625)
      positions=add_positions(positions,level)

    elif form=='paper':
      # grace.multi(rows=3,cols=1,vgap=.06,hgap=.05)
      for graph in grace.graphs:
        graph12.set_view(0.15, 0.7233333333333333, 0.4738095238095238, 0.95)
        graph13.set_view(0.15, 0.43666666666666654, 0.4738095238095238, 0.6633333333333332)
        graph23.set_view(0.15, 0.1499999999999999, 0.4738095238095238, 0.3766666666666666)
        positions=add_paper_positions(positions)
        positions.set_view(0.5,0.15,.65,0.95)

    grace.write_file('../../manuscript/figures/roles/roleplot_'+form+'_'+level+'.eps')
