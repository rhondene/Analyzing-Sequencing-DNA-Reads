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


def NaiveMatchRC(t,tc,s):
	occur_t = []
	occur_tc = []
	for i in range(len(s)-len(t)+1):
		match_t = True
		for j in range(len(t)):
			if s[i+j] != t[j]:
				match_t = False
				break
		if match_t:
			occur_t.append(i)
	for a in range(len(s)-len(tc)+1):
		match_tc = True
		for k in range(len(tc)):
			if s[a+k] != tc[k]:
				match_tc = False
				break
		if match_tc == True and tc != t:    ##will not double count if t is same as its reveerse complement
			occur_tc.append(a)
	return occur_t, occur_tc, len(occur_t), len(occur_tc)  ##pass these to 4 variables

x,y,a,b = NaiveMatchRC(t, tc, refgenome)

