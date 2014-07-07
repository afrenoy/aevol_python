#!/usr/bin/env python

import os
import numpy
from matplotlib.pyplot import *
import sys
import getopt

genomesize=None
fitness=None
metabolism=None
secretion=None
cRNA=None
ncRNA=None
fgenes=None
nfgenes=None
sizecRNA=None
sizencRNA=None
sizefgenes=None
sizenfgenes=None
nbrep=None
nbgen=None

def load_exp(path,start,end,gu=0):
    global genomesize
    global fitness
    global metabolism
    global secretion
    global cRNA
    global ncRNA
    global fgenes
    global nfgenes
    global sizecRNA
    global sizencRNA
    global sizefgenes
    global sizenfgenes
    global nbrep
    global nbgen
    
    allrep=[x for x in os.listdir(path) if os.path.isdir(os.path.join(path,x))]
    nbrep=len(allrep)
    
    if start==None:
        start=0
    if end==None:
        nbgen=int(open(path + '/' + allrep[0] + '/last_gener.txt','r').read())+1
        end=start+nbgen-1
    else:
        nbgen=end-start+1
        nbgeninfile=int(open(path + '/' + allrep[0] + '/last_gener.txt','r').read())+1
        assert(nbgeninfile>=nbgen)
        
    genomesize=numpy.zeros((nbrep,nbgen))
    fitness=numpy.zeros((nbrep,nbgen))
    metabolism=numpy.zeros((nbrep,nbgen))
    secretion=numpy.zeros((nbrep,nbgen))
    cRNA=numpy.zeros((nbrep,nbgen))
    ncRNA=numpy.zeros((nbrep,nbgen))
    fgenes=numpy.zeros((nbrep,nbgen))
    nfgenes=numpy.zeros((nbrep,nbgen))
    sizecRNA=numpy.zeros((nbrep,nbgen))
    sizencRNA=numpy.zeros((nbrep,nbgen))
    sizefgenes=numpy.zeros((nbrep,nbgen))
    sizenfgenes=numpy.zeros((nbrep,nbgen))
    
    for rep,replicate in enumerate(allrep):
        if gu==0:
            stat_fitness_glob=path + '/' + replicate + '/stats/stat_fitness_glob.out'
        elif gu==1:
            stat_fitness_glob=path + '/' + replicate + '/stats/stat_fitness_chromosome_glob.out'
        elif gu==2:
            stat_fitness_glob=path + '/' + replicate + '/stats/stat_fitness_plasmids_glob.out'
        a=numpy.loadtxt(stat_fitness_glob)
        genomesize[rep,0:nbgen]=a[start:end+1,3]
        fitness[rep,0:nbgen]=a[start:end+1,2]
        metabolism[rep,0:nbgen]=a[start:end+1,6]
        secretion[rep,0:nbgen]=a[start:end+1,9]
        if gu==0:
            stat_genes_glob=path + '/' + replicate + '/stats/stat_genes_glob.out'
        elif gu==1:
            stat_genes_glob=path + '/' + replicate + '/stats/stat_genes_chromosome_glob.out'
        elif gu==2:
            stat_genes_glob=path + '/' + replicate + '/stats/stat_genes_plasmids_glob.out'
        b=numpy.loadtxt(stat_genes_glob)
        cRNA[rep,0:nbgen]=b[start:end+1,1]
        ncRNA[rep,0:nbgen]=b[start:end+1,2]
        fgenes[rep,0:nbgen]=b[start:end+1,5]
        nfgenes[rep,0:nbgen]=b[start:end+1,6]
        sizecRNA[rep,0:nbgen]=b[start:end+1,3]
        sizencRNA[rep,0:nbgen]=b[start:end+1,4]
        sizefgenes[rep,0:nbgen]=b[start:end+1,7]
        sizenfgenes[rep,0:nbgen]=b[start:end+1,8]


def output_mean_std(path,start,end,gu=0):
    global genomesize
    global fitness
    global metabolism
    global secretion
    global cRNA
    global ncRNA
    global fgenes
    global nfgenes
    global sizecRNA
    global sizencRNA
    global sizefgenes
    global sizenfgenes
    global nbrep
    
    load_exp(path,start,end,gu)
    
    g=numpy.mean(genomesize,0)
    f=numpy.mean(fitness,0)
    m=numpy.mean(metabolism,0)
    s=numpy.mean(secretion,0)
    c=numpy.mean(cRNA,0)
    nc=numpy.mean(ncRNA,0)
    fonc=numpy.mean(fgenes,0)
    nfonc=numpy.mean(nfgenes,0)
    sizec=numpy.mean(sizecRNA,0)
    sizenc=numpy.mean(sizencRNA,0)
    sizefonc=numpy.mean(sizefgenes,0)
    sizenfonc=numpy.mean(sizenfgenes,0)
    g_s=numpy.std(genomesize,0)/numpy.sqrt(nbrep)
    f_s=numpy.std(fitness,0)/numpy.sqrt(nbrep)
    m_s=numpy.std(metabolism,0)/numpy.sqrt(nbrep)
    s_s=numpy.std(secretion,0)/numpy.sqrt(nbrep)
    c_s=numpy.std(cRNA,0)/numpy.sqrt(nbrep)
    nc_s=numpy.std(ncRNA,0)/numpy.sqrt(nbrep)
    fonc_s=numpy.std(fgenes,0)/numpy.sqrt(nbrep)
    nfonc_s=numpy.std(nfgenes,0)/numpy.sqrt(nbrep)
    sizec_s=numpy.std(sizecRNA,0)/numpy.sqrt(nbrep)
    sizenc_s=numpy.std(sizencRNA,0)/numpy.sqrt(nbrep)
    sizefonc_s=numpy.std(sizefgenes,0)/numpy.sqrt(nbrep)
    sizenfonc_s=numpy.std(sizenfgenes,0)/numpy.sqrt(nbrep)
    
    if gu==0:
        epath=path+'/'
    elif gu==1:
        epath=path+'/chr_'
    elif gu==2:
        epath=path+'/pla_'

    x=None
    if start==None:
        x=numpy.arange(0,len(g),1)
    elif end==None:
        x=numpy.arange(start,start+len(g),1)
    else:
        x=numpy.arange(start,end+1,1)
    # genome size
    figure()
    plot(x,g)
    fill_between(x,g-g_s,g+g_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'genomesize.pdf')
    close()
    # fitness
    figure()
    plot(x,f)
    fill_between(x,f-f_s,f+f_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'fitness.pdf')
    close()
    # metabolism
    figure()
    plot(x,m)
    fill_between(x,m-m_s,m+m_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'metabolism.pdf')
    close()
    # secretion
    figure()
    plot(x,s)
    fill_between(x,s-s_s,s+s_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'secretion.pdf')
    close()
    # coding RNA
    figure()
    plot(x,c)
    fill_between(x,c-c_s,c+c_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'codingRNA.pdf')
    close()
    # non coding RNA
    figure()
    plot(x,nc)
    fill_between(x,nc-nc_s,nc+nc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'noncodingRNA.pdf')
    close()
    # fonctional genes
    figure()
    plot(x,fonc)
    fill_between(x,fonc-fonc_s,fonc+fonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'fonctionalGenes.pdf')
    close()
    # non fonctional genes
    figure()
    plot(x,nfonc)
    fill_between(x,nfonc-nfonc_s,nfonc+nfonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'nonfonctionalGenes.pdf')
    close()
    # size coding RNA
    figure()
    plot(x,sizec)
    fill_between(x,sizec-sizec_s,sizec+sizec_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'sizecodingRNA.pdf')
    close()
    # size non coding RNA
    figure()
    plot(x,sizenc)
    fill_between(x,sizenc-sizenc_s,sizenc+sizenc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'sizenoncodingRNA.pdf')
    close()
    # size fonctional genes
    figure()
    plot(x,sizefonc)
    fill_between(x,sizefonc-sizefonc_s,sizefonc+sizefonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'sizefonctionalGenes.pdf')
    close()
    # size non fonctional genes
    figure()
    plot(x,sizenfonc)
    fill_between(x,sizenfonc-sizenfonc_s,sizenfonc+sizenfonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(epath+'sizenonfonctionalGenes.pdf')
    close()

def output_all(path,start,end,gu=0):
    global genomesize
    global fitness
    global metabolism
    global secretion
    global cRNA
    global ncRNA
    global fgenes
    global nfgenes
    global sizecRNA
    global sizencRNA
    global sizefgenes
    global sizenfgenes
    
    load_exp(path,start,end,gu)
    x=None
    if start==None:
        x=numpy.arange(0,numpy.shape(genomesize)[1],1)
    elif end==None:
        x=numpy.arange(start,start+numpy.shape(genomesize)[1],1)
    else:
        x=numpy.arange(start,end+1,1)
        
    if gu==0:
        epath=path+'/'
    elif gu==1:
        epath=path+'/chr_'
    elif gu==2:
        epath=path+'/pla_'
    
    # genome size
    figure()
    plot(x,genomesize.T)
    savefig(epath+'genomesize_all.pdf')
    close()
    # fitness
    figure()
    plot(x,fitness.T)
    savefig(epath+'fitness_all.pdf')
    close()
    # metabolism
    figure()
    plot(x,metabolism.T)
    savefig(epath+'metabolism_all.pdf')
    close()
    # secretion
    figure()
    plot(x,secretion.T)
    savefig(epath+'secretion_all.pdf')
    close()
    # coding RNA
    figure()
    plot(x,cRNA.T)
    savefig(epath+'codingRNA_all.pdf')
    close()
    # non coding RNA
    figure()
    plot(x,ncRNA.T)
    savefig(epath+'noncodingRNA_all.pdf')
    close()
    # fonctional genes
    figure()
    plot(x,fgenes.T)
    savefig(epath+'fonctionalGenes_all.pdf')
    close()
    # non fonctional genes
    figure()
    plot(x,nfgenes.T)
    savefig(epath+'nonfonctionalGenes_all.pdf')
    close()
    # size coding RNA
    figure()
    plot(x,sizecRNA.T)
    savefig(epath+'sizecodingRNA_all.pdf')
    close()
    # size non coding RNA
    figure()
    plot(x,sizencRNA.T)
    savefig(epath+'sizenoncodingRNA_all.pdf')
    close()
    # size fonctional genes
    figure()
    plot(x,sizefgenes.T)
    savefig(epath+'sizefonctionalGenes_all.pdf')
    close()
    # size non fonctional genes
    figure()
    plot(x,sizenfgenes.T)
    savefig(epath+'sizenonfonctionalGenes_all.pdf')
    close()

def output_prctile(path,start,end):
    load_exp(path,start,end)
    x=None
    if start==None:
        x=numpy.arange(0,len(g),1)
    elif end==None:
        x=numpy.arange(start,start+len(g),1)
    else:
        x=numpy.arange(start,end+1,1)
        
        
opts, args = getopt.getopt(sys.argv[1:], "s:e:g:a", ["start", "end"])
start=None
end=None
all=False
gus=set([])

for o, a in opts:
    if o=='-s':
        start=int(a)
    elif o=='-e':
        end=int(a)
    elif o=='-a':
        all=True
    elif o=='-g':
        gus.add(int(a))

if len(gus)==0:
    gus.add(0)

for arg in args:
    if all==True:
        for gu in gus:
            output_all(arg,start,end,gu)
    else:
        for gu in gus:
            output_mean_std(arg,start,end,gu)

