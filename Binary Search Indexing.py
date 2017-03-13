import bisect

class Index(object):
    ## preprocesses reference genome s into index table
    def __init__(self,s,k):
        self.k = k
        self.index = []
        for i in range(len(s)-k+1):    #adds every k-mer index to list
            self.index.append((s[i:i+k], i))    #appends a tuple with two items per entry: the kmer string, s[i:i+k], and its offset, i 
        self.index.sort()  
    
    def query(self,t):  #searches for a kmer from t in index of s
        kmer= t[:self.k]    
        i = bisect.bisect_left(self.index,(kmer, -1))  #this ensures that kmer occurs at leftmost offset
        hits = []                           # all the locations of the kmer in t matches s
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1]) 
            i+=1
        return hits


def queryIndex(t,s,index):
    k = index.k
    offsets = []
    for i in index.query(t):
        if t[k:]==s[i+k:i+len(t)]:   #verification if all of t exactly matches a portion in s
            offsets.append(i)
    return offsets
    
index = Index(s,k)
print(queryIndex(t,s,index))
