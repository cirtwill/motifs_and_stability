import sys
import os
import random

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

colors=ColorBrewerScheme('RdYlBu',n=15,reverse=True)  # The blue is very beautiful but maybe harder to see.

def prep_table_anova(dispdir,outfiles,roletype):
  anovadir={}
  for fil in outfiles:  # F, p by S and C
    S=int(fil.split('_')[-2])
    C=float(fil.split('_')[-1].split('.tsv')[0])
    f=open(dispdir+roletype+'/betadisper/'+fil,'r')
    for line in f:
      if line.split()[0]=='"Groups"':
        df=line.split()[1]
        SS=float(line.split()[2])
        MS=float(line.split()[3])
        F=float(line.split()[4])
        p=float(line.split()[5])
    f.close()
    anovadir[(S,C)]=((F,p))
  return anovadir

def prep_table_linmod(dispdir,outfiles,roletype):
  lmdir={}
  for fil in outfiles:  # beta(group), p by S and C
    S=int(fil.split('_')[-2])
    C=float(fil.split('_')[-1].split('.tsv')[0])
    f=open(dispdir+roletype+'/linmod/'+fil,'r')
    for line in f:
      if line.split()[0]=='"as.numeric(as.character(group))"':
        beta=float(line.split()[1])
        pval=float(line.split()[-1])
    f.close()
    lmdir[(S,C)]=((beta,pval))
  return lmdir

# Think I want to see what the differences are in significant vs. non-significant pairs...
def prep_table_Tukey(dispdir,outfiles,roletype):
  Tdir={}
  for fil in outfiles:  # beta(group), p by S and C
    S=int(fil.split('_')[-2])
    C=float(fil.split('_')[-1].split('.tsv')[0])
    Tdir[(S,C)]={'sig':[],'insig':[]}
    f=open(dispdir+roletype+'/Tukey/'+fil,'r')
    for line in f:
      if line.split()[0]!='"diff"':
        groups=line.split()[0][1:-1].split('-')
        diff=float(line.split()[1])
        lower=float(line.split()[2])
        upper=float(line.split()[3])
        p=float(line.split()[4])
        if p<0.05:
          Tdir[(S,C)]['sig'].append((groups,diff,p))
        else:
          Tdir[(S,C)]['insig'].append((groups,diff,p))
    f.close()

  return Tdir

def write_supptable(anovadir,lmdir,roletype):
  f=open('../../manuscript/tables/betadisper_'+roletype+'.tsv','w')
  f.write("&&\\multicolumn{2}{c|}{ANOVA}&\\multicolumn{2}{c||}{Regression}&&&\\multicolumn{2}{c|}{ANOVA}&\\multicolumn{2}{c|}{Regression}\\\\\n")
  f.write("S&C&F&$p$-value&$"+"\\"+"beta$&$p$-value&S&C&F&$p$-value&$"+"\\"+"beta$&$p$-value"+ "\\"+"\\"+"\n")
  f.write('\\hline')
  for (S,C) in sorted(anovadir):
    if S<80:
      f.write(str(S)+'&'+str(C)+'&'+str(round(anovadir[(S,C)][0],2))+'&'+str(round(anovadir[(S,C)][1],3))+'&')
      f.write(str(lmdir[(S,C)][0])+'&'+str(round(lmdir[(S,C)][1],3))+'&')
      S2=S+30
      f.write(str(S2)+'&'+str(C)+'&'+str(round(anovadir[(S2,C)][0],2))+'&'+str(round(anovadir[(S2,C)][1],3))+'&')
      f.write(str(lmdir[(S2,C)][0])+'&'+str(round(lmdir[(S2,C)][1],3))+'\\'+'\\'+'\n')
  f.close()

def format_graph(graph,key):
  graph.yaxis.bar.linewidth=1
  graph.xaxis.bar.linewidth=1
  graph.frame.linewidth=1
  graph.panel_label.configure(char_size=0.75,placement='iul',dx=0.03,dy=0.03)
  graph.world.ymin=0
  graph.world.xmin=0
  graph.world.xmax=11
  graph.xaxis.tick.major=2
  if key=='network':
    graph.world.ymax=1.75
    graph.yaxis.tick.major=0.4
    graph.yaxis.ticklabel.configure(format='decimal',prec=1,char_size=.5)
  else:
    graph.world.ymax=.12
    graph.yaxis.ticklabel.configure(format='decimal',prec=2,char_size=.5)
    graph.yaxis.tick.major=0.03
  graph.yaxis.label.configure(text='Difference between groups',char_size=.75,just=2)
  graph.xaxis.ticklabel.configure(char_size=.5,format='decimal',prec=0)
  graph.xaxis.label.configure(text='Difference between deciles',char_size=.75,just=2)

  graph.xaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)
  graph.yaxis.tick.configure(onoff='on',minor_ticks=0,major_size=.7,minor_size=.5,place='both',major_linewidth=1,minor_linewidth=1)

  return graph

def populate_graph(Tdict,graph,key):
  sigs=[]
  insigs=[]
  for (S,C) in Tdict:
    for ([g1,g2],diff,p) in Tdict[(S,C)]['sig']:
      gr1=float(g1)
      gr2=float(g2)
      x=max([gr1,gr2])-min([gr1,gr2])
      sigs.append((10*x-.2+random.uniform(0,0.15),abs(diff)))
    for ([g1,g2],diff,p) in Tdict[(S,C)]['insig']:
      gr1=float(g1)
      gr2=float(g2)
      x=max([gr1,gr2])-min([gr1,gr2])
      insigs.append((10*x+.2-random.uniform(0,0.15),abs(diff)))

  sdat=graph.add_dataset(sigs)
  sdat.line.linestyle=0
  sdat.symbol.configure(size=.1,color=13,fill_color=13,shape=1)
  if key=='species':
    sdat.legend='Significant differences'

  idat=graph.add_dataset(insigs)
  idat.line.linestyle=0
  idat.symbol.configure(size=.1,color=2,fill_color=2,shape=3)
  if key=='species':
    idat.legend='Non-significant differences'

  graph.legend.configure(box_linestyle=0,char_size=.5,loc=(1,0.11),loctype='world')

  return graph

def plot_differences(all_Tukeys):
  grace=MultiPanelGrace(colors=colors)
  grace.add_label_scheme('dummy',['Count','Network normalization','Species normalization'])
  grace.set_label_scheme('dummy')

  for key in sorted(all_Tukeys):
    graph=grace.add_graph(Panel)
    graph=format_graph(graph,key)
    graph=populate_graph(all_Tukeys[key],graph,key)
    graph.panel_label.configure(placement='ouc',char_size=.75,dy=.01)

  grace.multi(rows=3,cols=1,vgap=.05)
  grace.hide_redundant_labels()
  grace.set_col_yaxislabel(label='Difference between groups',char_size=1,perpendicular_offset=.07,col=0,rowSpan=(None,None))
  grace.write_file('../../manuscript/figures/Tukey_differences.eps')




def main():

  all_Tukeys={}
  dispdir='../../data/summaries/extorder_disps/'
  for roletype in os.listdir(dispdir): # count, network norm, species norm
    for testtype in os.listdir(dispdir+roletype): # betadisper (anova), Tukey HSD, linear model (linmod)
      outfiles=os.listdir(dispdir+roletype+'/'+testtype)
      if testtype=='betadisper':
        anovadir=prep_table_anova(dispdir,outfiles,roletype)
      elif testtype=='linmod':
        lmdir=prep_table_linmod(dispdir,outfiles,roletype)
      else:
        Tukeydir=prep_table_Tukey(dispdir,outfiles,roletype)
        all_Tukeys[roletype]=Tukeydir
    write_supptable(anovadir,lmdir,roletype)

  plot_differences(all_Tukeys)

if __name__ == '__main__':
  main()
