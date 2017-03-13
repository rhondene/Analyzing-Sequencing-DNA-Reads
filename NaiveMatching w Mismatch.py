# Sequencing-Reads
Python code for reading and basic sequencing of DNA reads
##My modified Naive Matching algorithm that allows for two mismatch per alignment

def revComp(t):       ##for reverse complement of a pattern t
	complement={ 'A': 'T', 'G': 'C', 'T':'A', 'C':'G' , 'N':'N' }
	tc = '' ##stores reverse complement of sequence t
	for base in t:
		
	        tc = complement[base] + tc
	return tc



def naiveTwomm(t,s):
	occur_t = []
	occur_tc = []
	tc = revComp(t)
	mm_t = 0
	mm_tc = 0

	for i in range(len(s)-len(t)+1):
		match_t = True
		for j in range(len(t)):
			if s[i+j] != t[j]: mm_t+=1
		if mm_t<3:
			occur_t.append(i)
		mm_t = 0                                  	##reset counter for each successive alignment	
	if tc != t:
		for a in range(len(s)-len(tc)+1):
			match_tc = True
			for k in range(len(tc)):
				if s[a+k] != tc[k]: mm_tc +=1
			if mm_tc<3 :
				occur_tc.append(a)
			mm_tc = 0	
	return occur_t, occur_tc, len(occur_t), len(occur_tc)   ##pass into 4 variables, especially for large genomes
