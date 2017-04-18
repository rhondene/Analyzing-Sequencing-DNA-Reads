#finds the minimum overlap between prefix of b in a, returns lenght of longest suffix
def overlap(a,b, min_length=3):  
    start = 0
    
    while True:
        start = a.find(b[:min_length], start ) #look in a for min prefix in b
        if start == -1: #no occurence of b in a, reach last part of a -1
            return 0
        if b.startswith(a[start:]): #verifies is prefix of b is qual to suffix of a
            return len(a)- start  #returns length of overlap
        start +=1
        
##Implement naive overlap map
from itertools import permutations as perm

def naiveOverlapMap(reads, k):  #k is min overlap length
    olaps = {}   
    for a,b in perm(reads,2):  #gets all possible pairs of reads in reads set
        olen = overlap(a,b, min_length = k)
        if olen > 0:
            olaps[(a,b)] = olen  #store each overlap as a tuble and assign a key
        return olaps
    
