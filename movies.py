#!/usr/bin/env python

import os
import numpy
import sys
import getopt
import png

xsize=0
ysize=0

def treat_generation(gfile,name,num):
    global xsize
    global ysize

    table=numpy.zeros((xsize,ysize))
    f = open(gfile,'r')
    for line in f:
        if len(line.split())==3:
            [x,y,z]=line.split()
            table[int(x),int(y)]=int(float(z))
    f.close()
    
    pix=[tuple(reduce(lambda x,y:x+y, [(int(j)*255,0,0) for j in i])) for i in table]

    namepng='_%s%04d.png'%(name,num)
    pngfile=open(namepng,'wb')
    wr=png.Writer(xsize,ysize)
    wr.write(pngfile,pix)
    pngfile.close()

        
opts, args = getopt.getopt(sys.argv[1:], "", [])

start=None
end=None
all=False
gus=set([])

#for o, a in opts:
#    if o=='-s':
#        start=int(a)
#    elif o=='-e':
#        end=int(a)
#    elif o=='-a':
#        all=True
#    elif o=='-g':
#        gus.add(int(a))


xsize=32
ysize=32
for arg in args:
    for i in range(0,1001):
        treat_generation(arg+'/stats/dump/secreted_amount_%04d.out'%i,'secretion',i)
    os.system("ffmpeg -pattern_type glob -i '*.png' -r 6  -vcodec png -s 320x320 -sws_flags neighbor -sws_dither none " + arg + "_secretion.mov")
    os.system("rm *.png")

