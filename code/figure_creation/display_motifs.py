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

# Positions and motifs
def format_netgraph(graph,form):
  graph.yaxis.bar.configure(linewidth=0,color=0)
  graph.xaxis.bar.configure(linewidth=0,color=0)
  graph.frame.linewidth=0

  graph.world.ymin=-.25
  graph.world.ymax=1
  graph.world.xmin=0
  if form=='oneway':
    graph.world.xmin=-2.5
    graph.world.xmax=10.5
  else:
    graph.world.xmax=13

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

def add_oneway_positions(graph):
  x=0
  y=0
  blackset={'1a':(0.5,0),'2a':(2,0.5),'3a':(3.5,0),'3b':(4,0.5),'3c':(4.5,0),
    '5b':(6.5,0),'5c':(7.5,0)}
  whiteset={'1b':(1,0.5),'2b':(2.5,0),'4a':(5,0.5),'4c':(6,0.5)}
  greyset={'1c':(1.5,0),'2c':(3,0.5),'4b':(5.5,0),'5a':(7,0.5)}

  # Links go prey, pred
  links=[('1a','1b'),('1b','1c'),
    ('2a','2b'),('2b','2c'),('2a','2c'),
    ('3a','3b'),('3b','3c'),('3c','3a'),
    ('4b','4a'),('4b','4c'),
    ('5b','5a'),('5c','5a')]

  graph=add_oneway_numbers(graph)

  graph.add_node_set(blackset,color=1,size=.75) 
  graph.add_node_set(whiteset,color=0,size=.75)
  graph.add_node_set(greyset,color='grey',size=.75)

  lynx=graph.add_directed_link_set(links,curvature=0,avoid_crossing_nodes=False,color=1)

  return graph

def add_twoway_positions(graph):
  x=0
  y=0
  blackset={'1a':(0.5,0.5),'1b':(1.5,0.5),'3a':(3.5,0.5),'4a':(5,0),'5a':(6.5,0.5),
    '6a':(8,0),'6b':(8.5,0.5),'6c':(9,0),'7a':(9.5,0.5)}
  whiteset={'1c':(1,0),'2a':(2.5,0.5),'3b':(4,0),'4b':(5.5,0.5),'5b':(7,0),
    '7b':(10,0),'8a':(11.5,0.5)}
  greyset={'2b':(2,0),'2c':(3,0),'3c':(4.5,0.5),'4c':(6,0),'5c':(7.5,0.5),
    '7c':(10.5,0.5),'8b':(11,0),'8c':(12,0)}

  # blackset={'6a':(0,-5),'7a':(x,-6),'8a':(0,-7),'10a':(x,-8.6),'11a':(x,-10),'12a':(0,-11)}
  # whiteset={'6b':(x,-4.6),'6c':(x,-5.4),
  #           '7b':(0,-5.6),'8b':(x,-6.6),'8c':(x,-7.4),'10b':(0,-9),
  #           '11b':(0,-9.6),'11c':(0,-10.4),'12b':(x,-10.6)}
  # greyset={'7c':(0,-6.4),'9a':(x,-8),'9b':(0,-7.6),'9c':(0,-8.4),'10c':(x,-9.4),
  #              '12c':(x,-11.4),'13a':(x,-12),'13b':(0,-11.6),'13c':(0,-12.4)}

  # Links go prey, pred
  links=[('1a','1b'),('1b','1a'),('1c','1a'),('1c','1b'),
    ('2b','2a'),('2c','2a'),('2b','2c'),('2c','2b'),
    ('3a','3c'),('3c','3a'),('3b','3c'),
    ('4a','4b'),('4a','4c'),('4c','4a'),
    ('5c','5a'),('5a','5c'),('5a','5b'),('5b','5c'),
    ('6b','6c'),('6c','6b'),('6a','6b'),('6b','6a'),('6a','6c'),('6c','6a'),
    ('7b','7a'),('7a','7b'),('7b','7c'),('7c','7b'),('7a','7c'),
    ('8a','8b'),('8a','8c'),('8b','8a'),('8c','8a')
    ]

  graph=add_twoway_numbers(graph)

  graph.add_node_set(blackset,color=1,size=.75) 
  graph.add_node_set(whiteset,color=0,size=.75)
  graph.add_node_set(greyset,color='grey',size=.75)

  lynx=graph.add_directed_link_set(links,curvature=0,avoid_crossing_nodes=False,color=1)

  return graph

def add_oneway_numbers(graph):
  # 1-6
  graph.add_drawing_object(DrawText,text='S1',y=.75,x=1,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='S2',y=.75,x=2.5,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='S3',y=.75,x=4,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='S4',y=.75,x=5.5,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='S5',y=.75,x=7,char_size=1,just=2,loctype='world')

  return graph

def add_twoway_numbers(graph):
  # 1-6
  graph.add_drawing_object(DrawText,text='D1',y=.75,x=1,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='D2',y=.75,x=2.5,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='D3',y=.75,x=4,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='D4',y=.75,x=5.5,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='D5',y=.75,x=7,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='D6',y=.75,x=8.5,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='D7',y=.75,x=10,char_size=1,just=2,loctype='world')
  graph.add_drawing_object(DrawText,text='D8',y=.75,x=11.5,char_size=1,just=2,loctype='world')

  return graph

# PCA data
def format_graph(graph):
  graph.yaxis.bar.configure(linewidth=1,color=0)
  graph.xaxis.bar.configure(linewidth=1,color=0)
  graph.frame.configure(linewidth=0,color=0)
  graph.xaxis.tick.onoff='off'
  graph.xaxis.ticklabel.onoff='off'
  graph.yaxis.tick.onoff='off'
  graph.yaxis.ticklabel.onoff='off'
  graph.world.xmin=0
  graph.world.xmax=10
  graph.world.ymin=-14
  graph.world.ymax=0

  return graph

def populate_graph(graph):

  graph.add_drawing_object(DrawText,text='S1: Three-species chain',loctype='world',x=0.5,y=-1,just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='S2: Omnivory',x=0.5,y=-2,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='S3: One-way three-species loop',x=0.5,y=-3,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='S4: Direct competition',x=0.5,y=-4,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='S5: Apparent competition',x=0.5,y=-5,loctype='world',just=0,char_size=1)

  graph.add_drawing_object(DrawText,text='D1: S4 with mutual predation among predators',loctype='world',x=0.5,y=-6,just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='D2: S5 with mutual predation among prey',x=0.5,y=-7,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='D3: S1 with mutual predation among top and intermediate spp.',x=0.5,y=-8,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='D4: S1 with mutual predation among bottom and intermediate spp.',x=0.5,y=-9,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='D5: S2 with mutual predation among top and bottom spp.',x=0.5,y=-10,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='D6: Two-way three-species loop',x=0.5,y=-11,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='D7: Three-species loop with two two-way links',x=0.5,y=-12,loctype='world',just=0,char_size=1)
  graph.add_drawing_object(DrawText,text='D8: S1 with mutual predation along both links.',x=0.5,y=-13,loctype='world',just=0,char_size=1)

  return graph

###########################################################################
## Main
###########################################################################

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['','','','D'])
grace.set_label_scheme('dummy')
positions1=grace.add_graph(NetworkPanel)
positions2=grace.add_graph(NetworkPanel)
graph12=grace.add_graph(Panel)

positions1=format_netgraph(positions1,'oneway')
positions1=add_oneway_positions(positions1)

positions2=format_netgraph(positions2,'twoway')
positions2=add_twoway_positions(positions2)

graph12=format_graph(graph12)
graph12=populate_graph(graph12)

grace.multi(rows=3,cols=1,vgap=.06,hgap=.05)
for graph in grace.graphs:
  positions1.set_view(0.05, 0.8, 0.95, 0.95)
  positions2.set_view(0.05,0.65,0.95,0.8)
  graph12.set_view(0.05,0.05,0.95,0.65)

grace.write_file('../../manuscript/figures/motifs.eps')
