#!/usr/bin/env python

import os
import numpy
import sys
import getopt
import png

def treat_generation(gfile,name,num,xsize,ysize,minv=0.,maxv=1.):
    table=numpy.zeros((xsize,ysize))
    f = open(gfile,'r')
    for line in f:
        if len(line.split())==3:
            [x,y,z]=line.split()
            value=float(z)
            if value<minv:
                value=minv
            if value>maxv:
                value=maxv
            table[int(x),int(y)]=value
    f.close()
    pix=[tuple(reduce(lambda x,y:x+y, [(int((j-minv)/(maxv-minv)*255),0,0) for j in i])) for i in table]

    namepng='_%s%04d.png'%(name,num)
    pngfile=open(namepng,'wb')
    wr=png.Writer(xsize,ysize)
    wr.write(pngfile,pix)
    pngfile.close()
        
opts, args = getopt.getopt(sys.argv[1:], "i:s:e:p:", ["inputname","start","end","patchsize"])

start=None
end=None
inputname=None
patchsize=20

for o, a in opts:
    if o=='-i' or o=='--inputname': # Name of the output to treat (basename of the matching dump files, for example 'fitness_metabolic'). 
        inputname=a
    elif o=='-s' or o =='--start': # (optional) First generation to treat
        start=int(a)
    elif o=='-e' or o=='--end': # (optional) Last generation to treat
        end=int(a)
    elif o=='-p' or o=='--patchsize': # (optional) Size of each location in pixels (default 20)
        patchsize=int(a)

if not inputname:
    print 'Error: you must specifiy option -i (--inputname)'
    exit(-1)

for arg in args:
    first=True
    for filename in os.listdir(arg+'/stats/dump/'):
        if filename.startswith(inputname):
            ngen=int(filename.split('_').pop().split('.')[0])
            if first: # Use this first file to infer size of the grid
                mx=-1
                my=-1
                first=False
                temp=open(arg+'/stats/dump/'+inputname+'_%04d.out'%ngen)
                for line in temp:
                    if len(line)<=3 or line.startswith('#'):
                        continue
                    x=int(line.split()[0])
                    y=int(line.split()[0])
                    if x>mx:
                        mx=x
                    if y>my:
                        my=y
                temp.close()
                print "Detected x size: %d, detected y size: %d"%(mx+1,my+1)
            if (start and ngen<start) or (end and ngen>end): # do not treat generations that are not in given interval if any
                continue
            treat_generation(arg+'/stats/dump/'+inputname+'_%04d.out'%ngen,inputname,ngen,mx+1,my+1)
    os.system("ffmpeg -pattern_type glob -i '*.png' -r 6  -vcodec png -s "+str((mx+1)*patchsize)+"x"+str((my+1)*patchsize)+" -sws_flags neighbor -sws_dither none " + arg + '/' + inputname + ".mov")
    os.system("rm *.png")

