#!/opt/local/bin/python

import os
import numpy
from matplotlib.pyplot import *
import sys

def get_stats_exp(path):
    allrep=[x for x in os.listdir(path) if os.path.isdir(os.path.join(path,x))]
    nbrep=len(allrep)
    nbgen=int(open(path + '/' + allrep[0] + '/last_gener.txt','r').read())+1
    genomesize=numpy.zeros((nbrep,nbgen))
    fitness=numpy.zeros((nbrep,nbgen))
    metabolism=numpy.zeros((nbrep,nbgen))
    secretion=numpy.zeros((nbrep,nbgen))
    for rep,replicate in enumerate(allrep):
        stat_fitness_glob=path + '/' + replicate + '/stats/stat_fitness_glob.out'
        a=numpy.loadtxt(stat_fitness_glob)
        genomesize[rep,:]=a[:,3]
        fitness[rep,:]=a[:,2]
        metabolism[rep,:]=a[:,6]
        secretion[rep,:]=a[:,9]
    return numpy.mean(genomesize,0),numpy.mean(fitness,0),numpy.mean(metabolism,0),numpy.mean(secretion,0),numpy.std(genomesize,0)/numpy.sqrt(nbrep),numpy.std(fitness,0)/numpy.sqrt(nbrep),numpy.std(metabolism,0)/numpy.sqrt(nbrep),numpy.std(secretion,0)/numpy.sqrt(nbrep)

def output_stats(path):
    [g,f,m,s,g_s,f_s,m_s,s_s]=get_stats_exp(path)
    x=numpy.arange(0,len(g),1)
    # genome size
    figure()
    plot(x,g)
    fill_between(x,g-g_s,g+g_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/genomesize.pdf')
    # fitness
    figure()
    plot(x,f)
    fill_between(x,f-f_s,f+f_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/fitness.pdf')
    # metabolism
    figure()
    plot(x,m)
    fill_between(x,m-m_s,m+m_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/metabolism.pdf')
    # secretion
    figure()
    plot(x,s)
    fill_between(x,s-s_s,s+s_s,alpha=0.5,facecolor='blue',edgecolor='none')
    savefig(path+'/secretion.pdf')

for arg in sys.argv[1:]:
    output_stats(arg)