#!/usr/bin/env python

import dataObj as do
import os
import numpy as np
import logging as log
import time
import pylab
import argparse
    
def readDNAfna(fname):
    if not os.path.exists(fname):
        print "ERROR: input file not found", fname
        exit()
        
    f = open(fname, "r")
    seq=""
    info=""
    # read every line of input file
    for l in f.readlines():
        if l.strip() != "":
            if l.strip()[0]==">":
                # extract title
                info=l.strip(">")
            else:
                seq += l.strip()
    data = do.dataObj(name=info, fname=fname, dataRaw=list(seq))
    return data

def specgram(data, n):
    steps = len(data)/n
    specs=np.zeros((steps,n/2))
    for s in range(steps):
        # FFT of the real part, A/T
        d1 = np.fft.rfft(data[(s*n):(s*n)+n].real).real
        
        # FFT of the imag part, C/G
        d2 = np.fft.rfft(data[(s*n):(s*n)+n].imag).real
        
        # add the two pieces together
        d=d1+d2
        
        # cut off half of the frequency range, since it repeats
        specs[s][:]=d[0:n/2]
    return specs

def chooseAnchors(spec, method):
    if method == 1:
        # original
        spec *= spec
        cutoff = args.anchorThresh * np.std(spec)
        
        t3=np.where(spec<cutoff)
        t4=np.where(spec>=cutoff)
        spec[t3]=0
        print "cutoff val:",cutoff
        print "below cutoff:",len(t3[0])
        print "over cutoff:",len(t4[0])
        
        if args.showPlots:
            spec2 = spec.copy()
            spec2[np.where(spec2>0)] = np.log(spec2[np.where(spec2>0)])
            img = pylab.imshow(np.transpose(spec2))
            pylab.colorbar(img)
            pylab.show()
        
        print spec.flatten().shape
        tmp=np.delete(spec.flatten(),np.where(spec.flatten()<cutoff))
        print tmp.flatten().shape
        print np.min(tmp), np.max(tmp)
        
        if args.showPlots:
            h = np.histogram(tmp,bins=100)
            pylab.bar(h[1][1:],h[0])
            pylab.autoscale()
            pylab.show()
        
        # return anchor points
        tmp[np.where(tmp>0)]=1    
        return tmp
    elif method == 2:
        # non squared
        cutoff = args.anchorThresh * np.std(spec)
        
        # swing everything positive
        spec = abs(spec)
        
        belowCut=np.where(spec<cutoff)
        aboveCut=np.where(spec>=cutoff)
        spec[belowCut]=0
        print "cutoff val:",cutoff
        #print "below cutoff:",len(belowCut[0])
        #print "over cutoff(Anchors):",len(aboveCut[0])
        
        if args.showPlots:
            spec2 = spec.copy()
            spec2[np.where(spec2>0)] = np.log(spec2[np.where(spec2>0)])
            img = pylab.imshow(np.transpose(spec2))
            pylab.colorbar(img)
            pylab.show()
        
        print spec.shape
        #tmp=np.delete(spec,np.where(spec<cutoff))
        spec[np.where(spec<cutoff)] = 0
        spec[np.where(spec>=cutoff)] = 1
        #print tmp.flatten().shape
        #print np.min(tmp), np.max(tmp)
        
        if args.showPlots:
            print "Final distribution of Anchors:"
            h = np.histogram(tmp,bins=100)
            pylab.bar(h[1][1:],h[0])
            pylab.autoscale()
            pylab.show()
        
        # return anchor points
        #tmp[np.where(tmp>0)]=1
        return spec
    else:
        print "you fail"
        exit()

def getPointsInBox(x_time,y_freq, anchorMap, searchBox):
    res=[]
    x=x_time
    print "x_time:%i  y_freq:%i" % (x_time,y_freq)
    while x <= (x_time+searchBox[0]):
        y=y_freq
        while y <= (y_freq+searchBox[1]):
            if anchorMap[x,y] == 1:
                # Distance to points
                # sqrt( (x1-x2)**2 + (y1-y2)**2 )
                # hash f1:f2:Dist: t1
                dist = np.sqrt( (x-x_time)**2 + (y-y_freq)**2 )
                #TODO: should the distance be rounded? to allow for small distance differences?
                f1 = y_freq
                f2 = y
                t1 = x_time
                #             [ Hash this | time1 ]
                res.append( ( (f1,f2,dist), t1 ) )
            y += 1
        x += 1
    return res

def getConstellations(anchorMap, searchBox=[10,10]):
    points = np.where(anchorMap==1)
    #print points, points, anchorMap.shape
    pointList = zip(points[0],points[1])
    
    
    for n,p in enumerate(pointList):
        #print len(anchorMap[0])
        #print anchorMap[p[0]:(p[0]+searchBox[0])][p[1]:(p[1]+searchBox[1])]
        #print np.where(anchorMap[p[0]:(p[0]+searchBox[0])][p[1]:(p[1]+searchBox[1])]>0)
        constellation = getPointsInBox(p[0],p[1],anchorMap, searchBox)
        print constellation
        print len(constellation)
        # p[0] is the time step
        # p[1] is the frequncy
        print n, p, anchorMap[p[0],p[1]]
        if n > 2:
            exit()

def main(args):
    dnaSeq = None
    stime=time.time()
    if args.rawInput is not None:
        dnaSeq = readDNAfna(args.rawInput)
        dnaSeq.dataTrans=dnaSeq.dataTrans[0:512*1000]
        dnaSeq.save("preprocessed/test2.b")
    elif args.preprocInput is not None:
        dnaSeq = do.dataObj(loadFile="preprocessed/test2.b")
    else:
        log.error("Must either have rawInput or preprocInput set")
        exit()
        
    log.info("loadtime:" + str(time.time()-stime) )    
    
    n=args.windowSize
    spec=specgram(dnaSeq.dataTrans[:n*args.specSize],n)
    
    if args.showPlots:
        print "Spectrogram:"
        img = pylab.imshow(np.transpose(spec))
        pylab.colorbar(img)
        pylab.show()
        
        print "Distribution of Spectrogram"
        h = np.histogram(spec,bins=100)
        pylab.bar(h[1][1:],h[0])
        pylab.autoscale()
        pylab.show()
    
    if True is False:
        spec2 = spec.copy()
        #cutoff = 55.75#3 * np.std(spec2)
        cutoff = args.anchorThresh * np.std(spec2)
        print "std:", np.std(spec2)
        print "cutoff:", args.anchorThresh * np.std(spec2)
    
        # create a set to subtract the center of the distribution
        t1 = np.where(spec2.flatten() > (-cutoff))
        t2 = np.where(spec2.flatten() < cutoff)
        
        # get the intersection and delete from spec
        c = np.intersect1d(t1[0],t2[0])
        tmp=np.delete(spec2.flatten(), c)
        
        # tmp now contains the left over pieces
        print "tmplen:",len(tmp)
        print np.min(tmp), np.max(tmp)
        h = np.histogram(tmp,bins=100)
        pylab.bar(h[1][1:],h[0])
        pylab.ylim([0,4000])
        pylab.show()
    
    # the anchorMap is a map of Ones and Zeros, every one representing an anchorpoint
    anchorMap = chooseAnchors(spec, 2)
    print "Selected %i Anchor points" % len(np.where(anchorMap>0)[0])
    
    getConstellations(anchorMap, searchBox=[10,10])
    
    
if __name__ == "__main__":
    #main("/home/madmaze/trash/DNA/Bacillus_anthracis/NC_003997.fna")
    parser = argparse.ArgumentParser(description="DNA fingerprinting")
    parser.add_argument("-r", "--rawInput", dest="rawInput", default=None, help="Raw Input file (Default: ./data/hs_alt_HuRef_chr22.fa)")
    parser.add_argument("-p", "--preprocInput", dest="preprocInput", default=None, help="Preprocessed Input file (Default: ./preprocessed/test3.b)")
    parser.add_argument("-w", "--windowSize", dest="windowSize", default=1024, type=int, help="Window size (Default: 1024)")
    parser.add_argument("-o", "--overlap", dest="overlap", default=512, type=int, help="Overlap size (Default: 512)")
    parser.add_argument("-s", "--specSize", dest="specSize", default=1000, type=int, help="SpecSize (Default: 1000)")
    parser.add_argument("-a", "--anchorThresh", dest="anchorThresh", default=10, type=int, help="Anchor threshold in sigmas (Default: 10)")
    parser.add_argument("--showPlots", dest="showPlots", action="store_true", help="Show plots (Default: False)")
    args = parser.parse_args()
    main(args)
    
