# Sequencing-Reads
Python code for reading and simple analyzing DNA reads from Sequencer, in fulfillment of the course Algorithms for DNA Sequencing by John's Hopkin via Coursera
import bm_preproc    ## script with preprocessing tables for pattern t
exec ('bm_preproc')

def boyer_moore(t, t_bm, s):  ##t is pattern and s is reference genome
    i = 0  ##keeps track number of shifts in text t
    occur = []
    num = 0
    x = 0
    total = len(s)-len(t)+1
    while i < len(s)-len(t)+1:  ##so that alignments of where t could start do not exceed len of s
        shift = 1  ## indicate the amount alignment shifts after comparision
        matched = True
        for j in range(len(t)-1, -1, -1): #loop through pattern  t , backwards
            if t[j] != s[i+j]:    ##upon mismatch
                skip_bc = t_bm.bad_character_rule(j, s[i+j])       ##calculate # shift for bc, pass in index and mismacth char
                skip_gs = t_bm.good_suffix_rule(j)                  ##calculates #of shift for gs rule, pass in index of mismatch
                ##we will shift by the largest shift from either rule to save most time
                shift = max(shift, skip_bc, skip_gs) 
                matched = False
                break                               #if mismatched then no further comparison needed, exit current alignement
        if matched:  
            occur.append(i)   
            skip_gs = t_bm.match_skip()            ##still do more skips with gs rule
            shift = max(shift, skip_gs) 
        i+=shift  ## update i by shift
        x+=1
    num = total - (i-x)
    return occur, num
    
    
 ###Approximate Matching using modified Boyer-Moore

def approx_match(t, s, n ):##n is number of mismtaches
    
    segment_length = int(round(len(t) / (n+1)))   ##computes length of each n+1 segment
    all_matches = set()            ##does not allow duplicates
    for i in range(n+1):      ## the outerloop handles matches in each ith segment
        start = i*segment_length
        end = min((i+1)*segment_length, len(t))  ##so that segment does not run past len of t
        t_bm = BoyerMoore(t[start:end], alphabet = 'AGCT')  ##does preprocessing tables for gs rule and bc rule
        matches = boyer_moore(t[start:end], t_bm, s)   ##returns where segment of t matches refeence text s
        
        for m in matches:
            if m < start or m-start+len(t) > len(s):
                continue
                
            mismatches = 0
            for j in range(0,start):     ##for the jth alignment
                if t[j] != s[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            for j in range(end,len(t)):   
                if t[j] != s[m-start+j]:
                    mismatches += 1
                    if mismatches > n:
                        break
            if mismatches <= n:
                all_matches.add(m - start)
    return list(all_matches)  ##turns the set all_match into an indexed list
    
