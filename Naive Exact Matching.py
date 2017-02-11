# Sequencing-Reads
Python code for reading and basic sequencing of DNA reads
###A code I wrote that the leftmost offset in a reference genome s that matches a pattern t and Reverse complement of t

##Matches only t with s
def naiveMatch(t,s):                    ## s is the reference genome and t is pattern 
	occur = []
	for i in range(len(s)- len(t) + 1):    ##iterates over alignment
		match = True                        ## intialise match as true
		for j in range(len(t)):            ##iterates over each char in t and s
			if s[i+j] != t[j]:         #compares each char of s and t
				match = False
				break
		if match:
		    occur.append(i)
	return occur
  
  ###matches t and its reverse complement in a reference genome
  #### Strandaware Naive Exact Match: Searches for occurences of pattern t
#### as well as reverse complement of t in reference genome s. Reports
####  occurences where reverse complement of t = t, whihc is treated as one occurrences.

def NaiveMatchRC(tc, t, s):
    occur_t = []
    occur_tc = []
    for i in range(len(s)-len(t)+1):        #loops over number of alignments
        match_t = True
        match_tc = True
        for j in range(len(t)):
            if s[i+j] != t[j]:
                match_t = False
                break
        for k in range(len(tc)):
            if s[i+k] != tc[k]:
                match_tc = False
                break
        if match_t:
            occur_t.append(i)
        if match_tc == True and tc != t:         #does not report occurence if t and its reverse complement are the same
            occur_tc.append(i)
    return occur_t, occur_tc
