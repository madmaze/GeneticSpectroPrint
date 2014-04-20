#!/usr/bin/env python
import numpy as np
import logging as log
import math
import cPickle as pickle

class dataObj:
    dataTrans=np.array([])
    
    def __init__(self, name = None, fname = None, dataRaw = [], loadFile=None):
       
        
        if loadFile==None:
            self.name = name.strip()
            self.fname = fname.strip()
            self.dataRaw = dataRaw
            
            log.debug("starting transcription of " + self.name + " " + str(len(dataRaw)))
            
            self.dataTrans=self.transcribe(self.dataRaw).astype(np.complex64)
            # pad to power of two
            #self.dataTrans = np.zeros(2**int(math.ceil(math.log(len(dataRaw),2))),dtype=np.complex32)
            #self.dataTrans
            log.debug("done transcription of " + self.name)
        else:
            self.loadPreprocessed(loadFile)
            
    
    def save(self, fname):
        pickle.dump((self.name, self.fname, self.dataTrans), open(fname,"wb"))
    
    def loadPreprocessed(self,fname):
        (self.name, self.fname, self.dataTrans) = pickle.load(open(fname,"rb"))
    
    def __len__(self):
        return len(self.dataRaw)
    
    def _transcribe(inst,i):
        i=i.lower()
        if i == "a":
            return np.complex(1)
        elif i == "u" or i == "t":
            return np.complex(-1)
        elif i == "c":
            return np.complex(1j)
        elif i == "g":
            return np.complex(-1j)
        elif i == "n":
            return np.complex(0)
        else:
            print "ERROR: transcoding."
            print type(i),i
            exit()
        return 0
    
    def _verifyTranscription(inst, i):
        if abs(i) == 1:
            return True
        elif abs(i) == 1j:
            return True
        elif i == 0:
            return True
        else:
            return False
    
    def verifyTranscription(self):
        log.debug("verifying transcription of " + self.name)
        i=0
        v = np.vectorize(self._verifyTranscription)
        tmp = v(self.dataTrans)
        log.debug("done verifying transcription")
        if np.count_nonzero(tmp) == len(self.dataTrans):
            return True
        else:
            return False
        
    def getTransPadded(self, padLen):
        padded=np.zeros(padLen,dtype=np.complex64)
        padded[0:len(self.dataTrans)]=self.dataTrans
        return padded[:]
    
    def getRawPadded(self, padLen):
        padded=["-"]*padLen
        padded[0:len(self.dataRaw)]=self.dataRaw
        return padded[:]
        
    def transcribe(self, raw):
        v = np.vectorize(self._transcribe)
        tmp = v(raw)
        return tmp 
        
    def __repr__(self):
        s=self.fname+"\n"
        return s