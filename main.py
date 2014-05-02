#!/usr/bin/env python

import dataObj as do
import os
import numpy as np
import logging as log
import time
import pylab
import argparse
import hashlib
import db
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import (generate_binary_structure,
                                      iterate_structure, binary_erosion)
    
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
    print "Length of genetic seq:", len(list(seq))
    data = do.dataObj(name=info, fname=fname, dataRaw=list(seq))
    return data

def specgram(data, n, overlap=512):
    #TODO: implement overlap of windows
    if False is True:
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
    else:
        # Method with overlap
        x = 0
        winSize = n
        specs = []
        while x+winSize <= len(data):
            # FFT of the real part, A/T
            d1 = np.fft.rfft(data[x:x+winSize].real).real
            
            # FFT of the imag part, C/G
            d2 = np.fft.rfft(data[x:x+winSize].imag).real
            
            # add the two pieces together
            d=d1+d2
            
            # cut off half of the frequency range, since it repeats
            specs.append(d[0:n/2].tolist())
            
            x += winSize-overlap
        specs=np.asarray(specs)
    return specs

def chooseAnchors(spec, method):
    stime=time.time()
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
        log.info("ChooseAnchor(1) time: %fs" % (time.time()-stime))
        tmp[np.where(tmp>0)]=1
        return tmp
    elif method == 2:
        # non squared (Better)
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
        log.info("ChooseAnchor(2) time: %fs" % (time.time()-stime))
        return spec
    elif method==3:
        # method as described in github worldveil/dejavu/fingerprint.py
        # generate binary mask
        binMask = generate_binary_structure(2,1)
        grownBinMask = iterate_structure(binMask, 15)
        
        filter = maximum_filter(spec, footprint=grownBinMask)
        #print filter
        local_max = filter == spec
        #print local_max
        background = (spec == 0)
        eroded_background = binary_erosion(background, structure=grownBinMask,
                                       border_value=1)

        # Boolean mask of arr2D with True at peaks
        detected_peaks = local_max - eroded_background
        return detected_peaks.astype(int)
    else:
        print "you fail"
        exit()

def getPointsInBox(x_time,y_freq, anchorMap, searchBox):
    res=[]
    x=x_time
    #print "x_time:%i  y_freq:%i" % (x_time,y_freq)
    #print anchorMap.shape
    while x <= anchorMap.shape[0]-1 and x <= (x_time+searchBox[0]):
        y=y_freq
        while y <= anchorMap.shape[1]-1 and y <= (y_freq+searchBox[1]):
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
                #res.append( ( (f1,f2,dist), t1 ) )
                #print str(f1),str(f2),str(dist)
                hash = hashlib.sha1("%s;%s;%s" % (str(f1),str(f2),str(dist)))
                res.append( ( hash.hexdigest(), t1 ) )
            y += 1
        x += 1
    return res

def getConstellations(anchorMap, searchBox=[10,10]):
    stime=time.time()
    points = np.where(anchorMap==1)
    #print points, points, anchorMap.shape
    pointList = zip(points[0],points[1])
    constellations=[]
    
    for n,p in enumerate(pointList):
        # p[0] is the time step
        # p[1] is the frequncy
        #print n, p, anchorMap[p[0],p[1]]
        
        # from current point, get all points/hashes within 
        constellations.extend(getPointsInBox(p[0],p[1],anchorMap, searchBox))
        #print constellation
        #print len(constellation)
        
        #if n > 10:
        #    return constellations
        
    log.info("getConstellations time: %fs" % (time.time()-stime))
    return constellations

def main(args):
    if args.reinitDB:
        """
        create database dnaindex; grant all privileges on dnaindex.* to 'dnafinger'@'localhost' identified by 'testpw';
        """
        dbcon = db.dbconn()
        dbcon.clearTable()
        dbcon.createTable()
        print "reinitialized DB.."
        exit()
        
    dnaSeq = None
    stime=time.time()
    
    if args.overlap == None:
        args.overlap = args.windowSize/2
    
    if args.rawInput is not None:
        dnaSeq = readDNAfna(args.rawInput)
        if args.specSize != 0:
            dnaSeq.dataTrans=dnaSeq.dataTrans[0:args.windowSize*args.specSize]

        #    dnaSeq.save("preprocessed/test2.b")
        #elif args.preprocInput is not None:
        #    dnaSeq = do.dataObj(loadFile="preprocessed/test2.b")
        #else:
        #    log.error("Must either have rawInput or preprocInput set")
        #    exit()
            
        log.info("loadtime: %fs" % (time.time()-stime) )
        stime=time.time()
        
        n=args.windowSize
        if args.specSize != 0:
            spec=specgram(dnaSeq.dataTrans[:n*args.specSize],n, args.overlap)
        else:
            spec=specgram(dnaSeq.dataTrans[:],n, args.overlap)
        
        log.info("Spectrogram time: %fs" % (time.time()-stime))
        stime=time.time()
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
        anchorMap = chooseAnchors(spec, args.anchorSelect)
        log.info("Selected %i Anchor points" % len(np.where(anchorMap>0)[0]))
        
        fingerprints = getConstellations(anchorMap, searchBox=[args.searchBox,args.searchBox])
        log.info("Generated %i fingerprints from anchorMap" % len(fingerprints))
        dbcon = db.dbconn()
        # insert into DB with filename as ID for now
        dbcon.bulkInset(fingerprints, args.rawInput)
    elif args.searchSeq != None:
        stime=time.time()
        # Lets search some stuff
        dnaSeq = readDNAfna(args.searchSeq)
        
        spec=specgram(dnaSeq.dataTrans[:],args.windowSize, args.overlap)
        
        anchorMap = chooseAnchors(spec, args.anchorSelect)
        
        fingerprints = getConstellations(anchorMap, searchBox=[args.searchBox,args.searchBox])
        
        dbcon = db.dbconn()
        
        res = dbcon.searchIndex(fingerprints)
        
        log.info("Total time: %fsec" % (time.time()-stime))
    elif args.DBstats is True:
        dbcon = db.dbconn()
        res = dbcon.getDBstats()
        log.info("total:%i unique:%i" % res)
    else:
        "Derp nothing to do.."
        exit()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DNA fingerprinting")
    parser.add_argument("-i", "--index", dest="rawInput", default=None, help="Raw Input file (Default: ./data/hs_alt_HuRef_chr22.fa)")
    parser.add_argument("--searchSeq", dest="searchSeq", default=None, help="Search a sequence (Default: ./test/test.seq)")
    parser.add_argument("-w", "--windowSize", dest="windowSize", default=1024, type=int, help="Window size (Default: 1024)")
    parser.add_argument("-o", "--overlap", dest="overlap", default=None, type=int, help="Overlap size (Default: windowSize/2)")
    parser.add_argument("-s", "--specSize", dest="specSize", default=0, type=int, help="SpecSize, number of window sizes. if zero then everything(Default: 0)")
    parser.add_argument("-a", "--anchorThresh", dest="anchorThresh", default=3, type=int, help="Anchor threshold in sigmas (Default: 3)")
    parser.add_argument("--anchorSelect", dest="anchorSelect", default=3, type=int, help="Anchor Selection method (Default: 3)")
    parser.add_argument("--searchBox", dest="searchBox", default=10, type=float, help="Search box size (Default: 10)")
    parser.add_argument("--showPlots", dest="showPlots", action="store_true", help="Show plots (Default: False)")
    parser.add_argument("--reinitDB", dest="reinitDB", action="store_true", help="Reinit DB (Default: False)")
    parser.add_argument("--DBstats", dest="DBstats", action="store_true", help="DB Stats (Default: False)")
    args = parser.parse_args()
    
    # setup logging
    log.basicConfig(level=log.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    main(args)
    
