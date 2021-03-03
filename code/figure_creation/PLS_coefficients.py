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
colors.add_color(204,102,119,'Tol7')
colors.add_color(136,34,85,'Tol8')
colors.add_color(170,68,153,'Tol9')
colors.add_color(221,221,221,'Tol10')


def format_graph(graph,normtype,coefs):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  graph.world.ymax=19
  graph.xaxis.tick.major=2
  graph.world.xmin=-8
  graph.world.xmax=3
  # elif normtype=="Degree":
  #   graph.world.xmin=-8
  #   graph.world.xmax=2.5
  # else:
  #   graph.world.xmin=-2
  #   graph.world.xmax=2

  if normtype=='Raw':
    graph.add_drawing_object(DrawText,text='Scaled',char_size=.75,just=2,x=-2.5,y=20,loctype='world')

  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
  graph.yaxis.label.configure(text='Predictor',char_size=1,just=2)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.5,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.set_spec_ticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],[],
    tick_labels=['S:C','C','S','STL','Degree','D8','D7','D6','D5','D4','D3','D2','D1',
        'S5','S4','S3','S2','S1'])

  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.5)
  # graph.xaxis.label.configure(text='Sum of coefficients',char_size=1,just=2)
  graph.xaxis.tick.configure(onoff='on',major_size=.5,minor_size=.3,place='both',major_linewidth=.75,minor_linewidth=.75)

  graph.panel_label.configure(char_size=0.75,placement='iul',dx=.03,dy=.01)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))
  # graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)/
  return graph

def format_Sgraph(graph,normtype,coefs):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  graph.world.ymax=19
  if normtype=='Network':
    graph.world.xmin=-6
    graph.world.xmax=9
    graph.xaxis.tick.major=3
  elif normtype=='Raw':
    graph.world.xmin=-175
    graph.world.xmax=125
    graph.xaxis.tick.major=50
  else:
    graph.world.xmin=-275
    graph.world.xmax=100
    graph.xaxis.tick.major=100
  # elif normtype=="Degree":
  #   graph.world.xmin=-8
  #   graph.world.xmax=2.5
  # else:
  #   graph.world.xmin=-2
  #   graph.world.xmax=2

  if normtype=='Raw':
    graph.yaxis.label.configure(text='Raw roles',char_size=.75,just=2,place='opposite')
    graph.add_drawing_object(DrawText,text='De-scaled',char_size=.75,just=2,x=-25,y=20,loctype='world')
  elif normtype=='Network':
    graph.yaxis.label.configure(text='Network normalization',char_size=.75,just=2,place='opposite')
  else:
    graph.yaxis.label.configure(text='Degree normalization',char_size=.75,just=2,place='opposite')

  graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=0)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.5,minor_size=.3,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.set_spec_ticks([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],[],
    tick_labels=['S:C','C','S','STL','Degree','D8','D7','D6','D5','D4','D3','D2','D1',
        'S5','S4','S3','S2','S1'])

  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.5)
  # graph.xaxis.label.configure(text='Sum of coefficients',char_size=1,just=2)
  graph.xaxis.tick.configure(onoff='on',minor_ticks=1,major_size=.5,minor_size=.3,place='both',major_linewidth=.75,minor_linewidth=.75)

  graph.panel_label.configure(char_size=0.75,placement='iul',dx=.03,dy=.01)
  # graph.legend.configure(char_size=0.75,box_linestyle=0,loctype='world',loc=(750,300))
  # graph.panel_label.configure(placement='iul',char_size=.75,dx=.02,dy=.01)/
  return graph

def populate_graph(graph,normtype,coefs):
  if normtype in ['Degree','Raw']:
    predlist=['S:C','C','S','STL','Degree','cD8','cD7','cD6','cD5','cD4','cD3','cD2','cD1',
        'cS5','cS4','cS3','cS2','cS1']
  elif normtype=='Network':
    predlist=['S:C','C','S','STL','Degree','zD8','zD7','zD6','zD5','zD4','zD3','zD2','zD1',
        'zS5','zS4','zS3','zS2','zS1']

  # Add grey backgrounds for stable motifs
  upper=graph.add_dataset([(-50,18.5),(50,18.5)])
  upper.line.linestyle=0
  upper.fill.configure(color='Tol1',type=2)
  S2=graph.add_dataset([(-50,16.5),(50,16.5)])  
  S2.fill.configure(color=0,type=2)
  S2.line.linestyle=0
  S4=graph.add_dataset([(-50,15.5),(50,15.5)])
  S4.fill.configure(color='Tol1',type=2)
  S4.line.linestyle=0
  S5=graph.add_dataset([(-50,13.5),(50,13.5)])
  S5.fill.configure(color=0,type=2)
  S5.line.linestyle=0

  y=1
  for pred in predlist:
    bar=graph.add_dataset([(0,y),(coefs[pred],y)],)
    bar.symbol.shape=0
    bar.line.configure(linewidth=6,linestyle=1,color='Tol4')
    y+=1


  # Add a line at 0
  graph.add_drawing_object(DrawLine,end=(0,0),start=(0,19),linewidth=1,linestyle=1,loctype='world')


  return graph

        # color = colorbar.z2color(pdf)
        # # you can change the opacity percentage of a single color, as well
        # # color.change_opacity(60)
        # graph.add_dataset([(x0,y0), (x1,y1)], SolidRectangle, color)

def populate_Sgraph(graph,normtype,coefs,scales):
  if normtype in ['Degree','Raw']:
    predlist=['S:C','C','S','STL','Degree','cD8','cD7','cD6','cD5','cD4','cD3','cD2','cD1',
        'cS5','cS4','cS3','cS2','cS1']
  elif normtype=='Network':
    predlist=['S:C','C','S','STL','Degree','zD8','zD7','zD6','zD5','zD4','zD3','zD2','zD1',
        'zS5','zS4','zS3','zS2','zS1']

  # Add grey backgrounds for stable motifs
  upper=graph.add_dataset([(-500,18.5),(500,18.5)])
  upper.line.linestyle=0
  upper.fill.configure(color='Tol2',type=2)
  S2=graph.add_dataset([(-500,16.5),(500,16.5)])  
  S2.fill.configure(color=0,type=2)
  S2.line.linestyle=0
  S4=graph.add_dataset([(-500,15.5),(500,15.5)])
  S4.fill.configure(color='Tol2',type=2)
  S4.line.linestyle=0
  S5=graph.add_dataset([(-500,13.5),(500,13.5)])
  S5.fill.configure(color=0,type=2)
  S5.line.linestyle=0

  y=1
  for pred in predlist:
    if pred!='S:C':
      bar=graph.add_dataset([(0,y),(coefs[pred]*scales[pred],y)],)
    else:
      bar=graph.add_dataset([(0,y),(coefs[pred]*scales['S']*scales['C'],y)],)
    bar.symbol.shape=0
    bar.line.configure(linewidth=6,linestyle=1,color='Tol3')
    y+=1


  # Add a line at 0
  graph.add_drawing_object(DrawLine,end=(0,0),start=(0,19),linewidth=1,linestyle=1,loctype='world')


  return graph

        # color = colorbar.z2color(pdf)
        # # you can change the opacity percentage of a single color, as well
        # # color.change_opacity(60)
        # graph.add_dataset([(x0,y0), (x1,y1)], SolidRectangle, color)

# Loadings are not for each component separately but all components summed. So, take only trailing value.
def read_coeffile(infile):
  subdict={}
  f=open(infile,'r')

  for line in f: 
    if line.split('\t')[0]!='"metadata$Persistence.1 comps"':
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
      coef=float(line.split()[-1]) # Overall coefficient in model including X components
      # cos=[]
      # for val in line.split()[1:]:
      #   cos.append(float(val))        
      # total=sum(cos)
      subdict[pred]=coef

  f.close()

  return subdict

def read_scalefile(scalefile):
  scales={}
  f=open(scalefile,'r')
  for line in f:
    if len(line.split())==2:
      pred=line.split()[0][1:-1]
      scale=float(line.split()[1])
      scales[pred]=scale
  f.close()
  return scales

###############################################################################################
###############################################################################################
#
#           Assemble the plots
#
###############################################################################################
###############################################################################################

rawfile='../../data/PLS_regression/raw_coefficients.tsv'
degfile='../../data/PLS_regression/degnorm_coefficients.tsv'
netfile='../../data/PLS_regression/netnorm_coefficients.tsv'
scalefile='../../data/PLS_regression/scales.tsv'

scales=read_scalefile(scalefile)
coefs={'Raw':{},'Degree':{},'Network':{}}
coefs['Raw']=read_coeffile(rawfile)
coefs['Degree']=read_coeffile(degfile)
coefs['Network']=read_coeffile(netfile)

grace=MultiPanelGrace(colors=colors)
grace.add_label_scheme('dummy',['A','B','C','D','E','F'])
grace.set_label_scheme('dummy')
for normtype in ['Raw','Degree','Network']:
  for scaltype in ['normal','descaled']:
    graph=grace.add_graph(Panel)
    if scaltype=='normal':
      graph=format_graph(graph,normtype,coefs)
      graph=populate_graph(graph,normtype,coefs[normtype])
    else:
      graph=format_Sgraph(graph,normtype,coefs)
      graph=populate_Sgraph(graph,normtype,coefs[normtype],scales)

# graph.set_view(0.1,0.45,0.9,0.95)
grace.multi(rows=3,cols=2,vgap=.03,hgap=.02)
grace.hide_redundant_yaxislabels()
grace.set_col_yaxislabel(rowspan=(None,None),col=0,label='Predictor',char_size=1,angle=90)
grace.set_row_xaxislabel(colspan=(None,None),row=2,label='Sum of coefficients',char_size=1,angle=90,perpendicular_offset=.04)

grace.write_file('../../manuscript/figures/PLS/total_coefficients.eps')

# Sgrace=MultiPanelGrace(colors=colors)
# Sgrace.add_label_scheme('dummy',['A) Raw roles','B) Degree normalization','C) Network normalization','S2: omnivory','S4: direct competition','S5: apparent competition',''])
# Sgrace.set_label_scheme('dummy')
# for normtype in ['Raw','Degree','Network']:
#   graph=Sgrace.add_graph(Panel)

# # graph.set_view(0.1,0.45,0.9,0.95)
# Sgrace.multi(rows=2,cols=2,vgap=.08,hgap=.04)
# Sgrace.hide_redundant_labels()
# # grace.set_col_yaxislabel(rowspan=(None,None),col=0,label='Predictor',char_size=1,angle=90)

# Sgrace.write_file('../../manuscript/figures/PLS/total_coefficients_descaled.eps')

