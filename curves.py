#!/usr/bin/env python

import os
import numpy
from matplotlib.pyplot import *
import sys
import getopt

varfilepos={
    'genomesize':('stat_fitness',3),
    'fitness':('stat_fitness',2),
    'metabolism':('stat_fitness',6),
    'secretion':('stat_fitness',9),
    'cRNA':('stat_genes',1),
    'ncRNA':('stat_genes',2),
    'fgenes':('stat_genes',5),
    'nfgenes':('stat_genes',6),
    'sizecRNA':('stat_genes',3),
    'sizencRNA':('stat_genes',4),
    'sizefgenes':('stat_genes',7),
    'sizenfgenes':('stat_genes',8)
}

guname={
    0:'',
    1:'_chromosome',
    2:'_plasmids'
}

class experiment:
    def __init__(self):
        self.data=dict()
    def load_exp(self,path,listvars=['genomesize','fitness','metabolism','secretion','cRNA','ncRNA','fgenes','nfgenes','sizecRNA','sizencRNA','sizefgenes','sizenfgenes'],start=None,end=None,gu=0):
        allrep=[x for x in os.listdir(path) if os.path.isdir(os.path.join(path,x))]
        self.nbrep=len(allrep)
        if start==None:
            start=0
        if end==None:
            self.nbgen=int(open(path + '/' + allrep[0] + '/last_gener.txt','r').read())+1
            end=start+self.nbgen-1
        else:
            self.nbgen=end-start+1
            nbgeninfile=int(open(path + '/' + allrep[0] + '/last_gener.txt','r').read())+1
            assert(nbgeninfile>=self.nbgen)
        for var in listvars:
            self.data[var]=numpy.zeros((self.nbrep,self.nbgen))
        fname=dict()
        for rep,replicate in enumerate(allrep):
            fcontent=dict()
            for fprefix in set([varfilepos[var][0] for var in listvars]):
                fname[fprefix]=path + '/' + replicate + '/stats/' + fprefix  + guname[gu] + '_glob.out'
                fcontent[fprefix]=numpy.loadtxt(fname[fprefix])
            for var in listvars:
                self.data[var][rep,0:self.nbgen]=fcontent[varfilepos[var][0]][start:end+1,varfilepos[var][1]]

    def output_mean_std(self,path,listvars=['genomesize','fitness','metabolism','secretion','cRNA','ncRNA','fgenes','nfgenes','sizecRNA','sizencRNA','sizefgenes','sizenfgenes'],start=None,end=None,gu=0,computeonly=False):
        self.load_exp(path,listvars,start,end,gu)
        mean=dict()
        std=dict()
        for var in listvars:
            mean[var]=numpy.mean(self.data[var],0)
            std[var]=numpy.std(self.data[var],0)/numpy.sqrt(self.nbrep)
        x=None
        if start==None:
            x=numpy.arange(0,self.nbgen,1)
        elif end==None:
            x=numpy.arange(start,start+self.nbgen,1)
        else:
            x=numpy.arange(start,end+1,1)

        if (computeonly):
            return (x,{var:(mean[var],std[var]) for var in listvars})

        # Produce figures
        for var in listvars:
            figure()
            plot(x,mean[var])
            fill_between(x,mean[var]-std[var],mean[var]+std[var],alpha=0.5,facecolor='blue',edgecolor='none')
            savefig(path+'/'+var+guname[gu]+'.pdf')
            close()

    def output_all(self,path,listvars=['genomesize','fitness','metabolism','secretion','cRNA','ncRNA','fgenes','nfgenes','sizecRNA','sizencRNA','sizefgenes','sizenfgenes'],start=None,end=None,gu=0):
        self.load_exp(path,listvars,start,end,gu)
        x=None
        if start==None:
            x=numpy.arange(0,self.nbgen,1)
        elif end==None:
            x=numpy.arange(start,start+self.nbgen,1)
        else:
            x=numpy.arange(start,end+1,1)

        # Produce figures
        for var in listvars:
            figure()
            plot(x,self.data[var].T)
            savefig(path+'/'+var+guname[gu]+'_all.pdf')
            close()

    def output_prctile(self,path,listvars=['genomesize','fitness','metabolism','secretion','cRNA','ncRNA','fgenes','nfgenes','sizecRNA','sizencRNA','sizefgenes','sizenfgenes'],start=None,end=None):
        load_exp(path,listvars,start,end)
        x=None
        if start==None:
            x=numpy.arange(0,len(g),1)
        elif end==None:
            x=numpy.arange(start,start+len(g),1)
        else:
            x=numpy.arange(start,end+1,1)


def print_help():
    print "Help available on http://github.com/frenoy/aevol_python"

opts, args = getopt.getopt(sys.argv[1:], "s:e:g:cao:h", ["start", "end","genunit","compare","all","offset","help"])
start=None
end=None
all=False
gus=set([])
compare=False
offsets=[]

for o, a in opts:
    if o=='-s' or o=='--start':
        start=int(a)
    elif o=='-e' or o=='--end':
        end=int(a)
    elif o=='-c' or o=='--compare':
        compare=True
    elif o=='-a' or o=='--all':
        all=True
    elif o=='-g' or o=='--genunit':
        gus.add(int(a))
    elif o=='-o' or o=='--offset':
        offsets.append(int(a))
    elif o=='-h' or o=='--help':
        print_help()
        exit()

if len(gus)==0:
    gus.add(0)

if compare:
    epath=os.getcwd()
    #epath=os.path.dirname(args[0])
    for gu in gus:
        r=dict()
        for (i,arg) in enumerate(args):
            e=experiment()
            (x,r[args])=e.output_mean_std(arg,start=start+offsets[i],end=end+offsets[i],gu=gu,computeonly=True)
            #(x,r[arg])=output_mean_std(arg,start+offsets[i],end+offsets[i],gu,True)
        for feature in r[args[0]].keys():
            figure()
            for arg in args:
                data=r[arg][feature]
                h=plot(x,data[0])
                color=h[0].get_color()
                fill_between(x,data[0]-data[1],data[0]+data[1],alpha=0.5,facecolor=color,edgecolor='none')
            foldernames=[os.path.basename(arg) for arg in args]
            #foldernames=['M-S 5k','M-S 420k','N-S']
            legend(foldernames,loc='upper left')
            xlabel('Generations')
            ylabel('Average secretion')
            savefig(epath+'/'+feature+'.pdf')
            close()
else:
    for arg in args:
        if all==True:
            for gu in gus:
                e=experiment()
                e.output_all(arg,start=start,end=end,gu=gu)
        else:
            for gu in gus:
                e=experiment()
                e.output_mean_std(arg,start=start,end=end,gu=gu)

