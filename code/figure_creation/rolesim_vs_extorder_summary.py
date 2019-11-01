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

def read_lmfile(infile):
  f=open(infile,'r')
  for line in f:
    if line.split()[0]=='"Dissim"':
      dissim=float(line.split()[1])      
      dissim_p=float(line.split()[-1])
    elif line.split()[0]=='"Path"':
      path=float(line.split()[1])
      path_p=float(line.split()[-1])
    elif line.split()[0]=='"(Intercept)"':
      inter=float(line.split()[1])
      inter_p=float(line.split()[-1])
    elif line.split()[0]=='"Dissim:Path"':
      interact=float(line.split()[1])
      interact_p=float(line.split()[-1])
  f.close()

  lmdict={'dissim':(dissim,dissim_p),'path':(path,path_p),
    'interact':(interact,interact_p),'intercept':(inter,inter_p)}

  return lmdict

def format_graph(graph,form,beta):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1

  graph.world.ymin=0
  if beta=='dissim':
    graph.world.ymax=10
    graph.world.ymin=-2
    graph.yaxis.tick.major=2
  else:
    graph.world.ymin=-3
    graph.world.ymax=2
    graph.yaxis.tick.major=1    
  graph.world.xmin=45
  graph.world.xmax=105
  graph.yaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)
  graph.xaxis.tick.major=10
  graph.xaxis.ticklabel.configure(format='decimal',prec=0,char_size=.75)

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  graph.panel_label.configure(char_size=0.75,placement='ouc',dx=0.0,dy=0.01)
  if beta=='dissim':
    graph.yaxis.label.configure(text='Coefficient',char_size=1,just=2)
  if beta=='path':
    graph.xaxis.label.configure(text='Species richness',char_size=1,just=2)
  return graph


def populate_graph(graph,level,datadict,beta):
  if level=='point':
    c=datadict[50].keys()[0]
    pointy=(50,datadict[50][c]['dissim'][0])
    dat=graph.add_dataset([pointy])
    if datadict[50][c]['dissim'][1]>0.05:
      shap=3
    else:
      shap=1
    dat.symbol.configure(fill_color=colorbar.z2color(c),shape=shap,size=0.75,color=1)

  else:
    for S in datadict:
      for c in datadict[S]:
        slope=datadict[S][c][beta][0]
        pval=datadict[S][c][beta][1]
        point=(S,slope)
        dat=graph.add_dataset([point])
        dat.line.linestyle=0
        if pval>0.05:
          shap=3
        else:
          shap=1
        dat.symbol.configure(fill_color=colorbar.z2color(c),shape=shap,size=0.75,color=1)
        if beta=='dissim':
          if S==50 and c==0.08 and shap==1:
            dat.legend='Significant'
          if shap==3:
            dat.legend='Non-significant'
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
for s in os.listdir(datadir):
  datadict[int(s)]={}
  options=os.listdir(datadir+s)
  lmfiles=[x for x in options if '_lm_' in x]
  for lmfile in lmfiles:
    c=float(lmfile.split('_')[-1].split('.tsv')[0])
    datadict[int(s)][c]=read_lmfile(datadir+s+'/'+lmfile)

for form in ['talk']:
  for level in ['axis','point','full']:
    grace=MultiPanelGrace(colors=colors)
    grace.add_label_scheme('dummy',['','Role dissimilarity','Path length','Interaction'])
    grace.set_label_scheme('dummy')
    colorbar = grace.add_graph(ElLogColorBar,domain=(0,0.2),
                               scale=LINEAR_SCALE,autoscale=False)
    colorbar.yaxis.label.configure(text="Connectance",char_size=1,just=2)
    colorbar.yaxis.tick.configure(major=0.02,major_size=.5,minor_ticks=0)
    colorbar.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.75)
    dissimgraph=grace.add_graph(Panel)
    dissimgraph=format_graph(dissimgraph,form,'dissim')

    pathgraph=grace.add_graph(Panel)
    pathgraph=format_graph(pathgraph,form,'path')

    intergraph=grace.add_graph(Panel)
    intergraph=format_graph(intergraph,form,'interaction')

    if level=='point':
      dissimgraph=populate_graph(dissimgraph,level,datadict,'dissim')
    if level=='full':
      dissimgraph=populate_graph(dissimgraph,level,datadict,'dissim')
      pathgraph=populate_graph(pathgraph,level,datadict,'path')
      intergraph=populate_graph(intergraph,level,datadict,'interact')

      dissimgraph.legend.configure(char_size=.5,box_linestyle=0,fill=0,fill_pattern=0,
        loc=(47,1.5),loctype='world')
   # graph.legend.configure(box_linestyle=0,fill=0,fill_pattern=0,char_size=.75,
   #  loc=(100,.75),loctype='world')
    grace.multi(rows=1,cols=4,vgap=.08,hgap=.04)
    dissimgraph.set_view(0.05,0.2,0.3,0.5)
    pathgraph.set_view(0.35,0.2,0.6,0.5)
    intergraph.set_view(0.65,0.2,0.9,0.5)
    colorbar.set_view(0.925,0.2,0.975,0.5)

    # grace.set_row_xaxislabel(colspan=(0,1),row=0,label="Species richness",just=2,char_size=1,perpendicular_offset=0.05)


    grace.write_file('../../manuscript/figures/extinction_order/dissimilarity_lm_summary_'+form+'_'+level+'.eps')
