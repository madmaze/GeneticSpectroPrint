import MySQLdb
import time
import logging as log
#from MySQLdb.cursors import DictCursor

class dbconn():
    db=None
    
    def __init__(self):
        self.db = MySQLdb.connect(host="localhost",
                                  user="dnafinger",
                                  passwd="testpw",
                                  db="dnaindex")
        
    def clearTable(self):
        cur = self.db.cursor()
        cur.execute("DROP TABLE IF EXISTS fingerprints")
        self.db.commit()
        cur.close()
        
    def createTable(self):
        cur = self.db.cursor()
        cur.execute("""
        CREATE TABLE fingerprints(
                    hash varchar(40) not null,
                    offset int not null,
                    dnaname varchar(128) not null,
                    PRIMARY KEY (hash, offset, dnaname),
                    INDEX (hash)
                    ) ENGINE=INNODB""")
        self.db.commit()
        cur.close()
        
    def insertFingerprint(self, hash, offset, identifier):
        cur = self.db.cursor()
        cur.execute("""INSERT IGNORE INTO fingerprints (hash, offset, dnaname) values ("%s",%i,"%s")""" % (hash, offset, identifier))
        self.db.commit()
        cur.close()
        
    def bulkInset(self, fingerprints, identifier):
        stime=time.time()
        log.info("Inserting fingerprints into database..")
        cur = self.db.cursor()
        for h,o in fingerprints:
            cur.execute("""INSERT IGNORE INTO fingerprints (hash, offset, dnaname) values ("%s",%i,"%s")""" % (h, o, identifier))
        self.db.commit()
        cur.close()
        log.info("fingerprintInsert time:"+str(time.time()-stime))
    
    def getDBstats(self):
        cur=self.db.cursor()
        query="select count(hash) as hashes, count(distinct(hash)) as uhashes from fingerprints;"
        cur.execute(query)
        
        for res in cur:
            return res
        
    def searchIndex(self, fingerprints):
        log.info("Number of fingerprints to search: %i" % len(fingerprints))
        stime=time.time()
        cur = self.db.cursor()
        
        fpHashtable = {}
        for fp,offset in fingerprints:
            fpHashtable[fp]=offset
        
        hashList,_ = zip(*fingerprints)
        hashString = "','".join(hashList)
        query = "SELECT * FROM fingerprints WHERE hash IN ('%s')" %hashString 
        
        cur.execute(query)
        #r=[]
        #for (hash,off,id) in cur:
        #    r.append((hash,off,id))

        #log.info("searchIndex[db query] time: %fs" % (time.time()-stime))
        #stime=time.time()
        
        cnt=0
        results={}
        for (hash,off,id) in cur:
            cnt+=1
            if not results.has_key(id):
                results[id]={}
                
            t_dist=off-fpHashtable[hash]
            if not results[id].has_key(t_dist):
                results[id][t_dist]=[]
            results[id][t_dist].append((hash,off))
            
                
        log.info("Num of Results: %i" % cnt)
        bestRes={}
        # if retBestAll == True, then multiple matches from one ID are allowed (slower)
        # if retBestAll == False, then only the best match for each ID is returned
        retBestAll=False
        for id in results.keys():
            m=(0,0)
            if retBestAll is True:
                for dist in results[id].keys():
                    bestRes[(id,dist)]=(len(results[id][dist]),dist)
            else:
                for dist in results[id].keys():
                    if m[0] < len(results[id][dist]):
                        m = (len(results[id][dist]),dist)
                bestRes[id]=m
        
        sortedRes = sorted(bestRes.items(), key=lambda t: t[1][0], reverse=True)
        print "TOP 5:"
        print "Score\tOffset\tID"
        for n,r in enumerate(sortedRes):
            if n<5:
                print "%i\t%i\t%s" % (r[1][0],r[1][1],r[0])
                
        log.info("searchIndex[anaylze result] time: %fs" % (time.time()-stime))
        