#!/usr/bin/env python

import os
import numpy
from matplotlib.pyplot import *
import sys
import getopt

def get_stats_exp(path,start,end):
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
        stat_fitness_glob=path + '/' + replicate + '/stats/stat_fitness_glob.out'
        a=numpy.loadtxt(stat_fitness_glob)
        genomesize[rep,0:nbgen]=a[start:end+1,3]
        fitness[rep,0:nbgen]=a[start:end+1,2]
        metabolism[rep,0:nbgen]=a[start:end+1,6]
        secretion[rep,0:nbgen]=a[start:end+1,9]
        stat_genes_glob=path + '/' + replicate + '/stats/stat_genes_glob.out'
        b=numpy.loadtxt(stat_genes_glob)
        cRNA[rep,0:nbgen]=b[start:end+1,1]
        ncRNA[rep,0:nbgen]=b[start:end+1,2]
        fgenes[rep,0:nbgen]=b[start:end+1,5]
        nfgenes[rep,0:nbgen]=b[start:end+1,6]
        sizecRNA[rep,0:nbgen]=b[start:end+1,3]
        sizencRNA[rep,0:nbgen]=b[start:end+1,4]
        sizefgenes[rep,0:nbgen]=b[start:end+1,7]
        sizenfgenes[rep,0:nbgen]=b[start:end+1,8]
    return (
            numpy.mean(genomesize,0),numpy.mean(fitness,0),numpy.mean(metabolism,0),numpy.mean(secretion,0),
            numpy.mean(cRNA,0),numpy.mean(ncRNA,0),numpy.mean(fgenes,0),numpy.mean(nfgenes,0),numpy.mean(sizecRNA,0),numpy.mean(sizencRNA,0),numpy.mean(sizefgenes,0),numpy.mean(sizenfgenes,0),
            numpy.std(genomesize,0)/numpy.sqrt(nbrep),numpy.std(fitness,0)/numpy.sqrt(nbrep),numpy.std(metabolism,0)/numpy.sqrt(nbrep),numpy.std(secretion,0)/numpy.sqrt(nbrep),
            numpy.std(cRNA,0)/numpy.sqrt(nbrep),numpy.std(ncRNA,0)/numpy.sqrt(nbrep),numpy.std(fgenes,0)/numpy.sqrt(nbrep),numpy.std(nfgenes,0)/numpy.sqrt(nbrep),numpy.std(sizecRNA,0)/numpy.sqrt(nbrep),numpy.std(sizencRNA,0)/numpy.sqrt(nbrep),numpy.std(sizefgenes,0)/numpy.sqrt(nbrep),numpy.std(sizenfgenes,0)/numpy.sqrt(nbrep)
           )

def output_stats(path,start,end):
    [g,f,m,s,c,nc,fonc,nfonc,sizec,sizenc,sizefonc,sizenfonc,g_s,f_s,m_s,s_s,c_s,nc_s,fonc_s,nfonc_s,sizec_s,sizenc_s,sizefonc_s,sizenfonc_s]=get_stats_exp(path,start,end)
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
    savefig(path+'/genomesize.pdf')
    close()
    # fitness
    figure()
    plot(x,f)
    fill_between(x,f-f_s,f+f_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/fitness.pdf')
    close()
    # metabolism
    figure()
    plot(x,m)
    fill_between(x,m-m_s,m+m_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/metabolism.pdf')
    close()
    # secretion
    figure()
    plot(x,s)
    fill_between(x,s-s_s,s+s_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/secretion.pdf')
    close()
    # coding RNA
    figure()
    plot(x,c)
    fill_between(x,c-c_s,c+c_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/codingRNA.pdf')
    close()
    # non coding RNA
    figure()
    plot(x,nc)
    fill_between(x,nc-nc_s,nc+nc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/noncodingRNA.pdf')
    close()
    # fonctional genes
    figure()
    plot(x,fonc)
    fill_between(x,fonc-fonc_s,fonc+fonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/fonctionalGenes.pdf')
    close()
    # non fonctional genes
    figure()
    plot(x,nfonc)
    fill_between(x,nfonc-nfonc_s,nfonc+nfonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/nonfonctionalGenes.pdf')
    close()
    # size coding RNA
    figure()
    plot(x,sizec)
    fill_between(x,sizec-sizec_s,sizec+sizec_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/sizecodingRNA.pdf')
    close()
    # size non coding RNA
    figure()
    plot(x,sizenc)
    fill_between(x,sizenc-sizenc_s,sizenc+sizenc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/sizenoncodingRNA.pdf')
    close()
    # size fonctional genes
    figure()
    plot(x,sizefonc)
    fill_between(x,sizefonc-sizefonc_s,sizefonc+sizefonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/sizefonctionalGenes.pdf')
    close()
    # size non fonctional genes
    figure()
    plot(x,sizenfonc)
    fill_between(x,sizenfonc-sizenfonc_s,sizenfonc+sizenfonc_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/sizenonfonctionalGenes.pdf')
    close()

opts, args = getopt.getopt(sys.argv[1:], "s:e:", ["start", "end"])
start=None
end=None

for o, a in opts:
    if o=='-s':
        start=int(a)
    elif o=='-e':
        end=int(a)
    
for arg in args:
    output_stats(arg,start,end)
