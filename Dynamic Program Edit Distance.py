#Edit distance includes all possible edits subs and indels
#implementing edit distance recursively slooow
#longest prefix means all bases upto but not include last base

n=0
def editDistRecur(x,y):
    global n
    if len(x) ==0 :
        return  len(y)  #edit dist between b and a empty string is the length of b
    elif len(y)==0:
        return len(x)
    else:
        distHor = editDistRecur(x[:-1], y) +1  #edit longest prefix of x into y then insert the final base in x to match y
        distVer = editDistRecur(x, y[:-1]) + 1
        if x[-1] == y[-1]:
            distDiag = editDistRecur(x[:-1], y[:-1]) #edit longest prefices of both x and y
        else:
            distDiag = editDistRecur(x[:-1], y[:-1])+1   #edit x nto y then do one more sub to turn last base
    return min(distHor, distVer, distDiag)


#dynamic programming proves a faster method to find edit distance
def editDistDP(x,y): 
    D = [] 
    for i in range(len(x)+1) :  #build the x+1 by y+1 matrix
        D.append([0]*(len(y)+1))
    
    for i in range(len(x)+1):  #initialise first row and col as as ascending 0,1,2,3,4... in case an empty string is given
        D[i][0] = i             #if len(y) = zero then len of of x is return 
    for i in range(len(y)+1):
        D[0][i] = i
    
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1]+1
            distVer = D[i-1][j]+1
            if x[i-1]==y[j-1]:
                distDiag = D[i-1][j-1]+1
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    return D[-1][-1]  #last row and column is edit dist
