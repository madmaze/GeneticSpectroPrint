import MySQLdb
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
        cur.execute("""CREATE TABLE fingerprints(
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
        print "Inserting fingerprints into database.."
        cur = self.db.cursor()
        for h,o in fingerprints:
            cur.execute("""INSERT IGNORE INTO fingerprints (hash, offset, dnaname) values ("%s",%i,"%s")""" % (h, o, identifier))
        self.db.commit()
        cur.close()
        
    def searchIndex(self, fingerprints):
        print "Number of fingerprints:",len(fingerprints)
        cur = self.db.cursor()
        
        fpHashtable = {}
        for fp,offset in fingerprints:
            fpHashtable[fp]=offset
        
        hashList,_ = zip(*fingerprints)
        hashString = "','".join(hashList)
        query = "SELECT * FROM fingerprints WHERE hash IN ('%s')" %hashString 
        
        cur.execute(query)
        r=[]
        for (hash,off,id) in cur:
            r.append((hash,off,id))
        print "done with db.."
        
        cnt=0
        results={}
        for (hash,off,id) in r:
            cnt+=1
            if not results.has_key(id):
                results[id]={}
                
            t_dist=off-fpHashtable[hash]
            if not results[id].has_key(t_dist):
                results[id][t_dist]=[]
            results[id][t_dist].append((hash,off))
            
            #if cnt%10000 == 0 :
            #    print cnt
                
        print "results:",cnt
        print len(results.keys())
        bestRes={}
        for id in results.keys():
            m=(0,0)
            for dist in results[id].keys():
                if m[0] < len(results[id][dist]):
                    m = (len(results[id][dist]),dist)
            bestRes[id]=m
        
        print bestRes
        