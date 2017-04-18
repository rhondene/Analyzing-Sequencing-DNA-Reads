#Implementing global alignnment which applies differential penalty for edits, by modifying editDistDP
alphabet = ['A','C', 'G','T']
# 4 is for transversion, 2 for transition, 8 is for skips or deletions
score = [[0,4,2,4,8],\
         [4,0,4,2,8],\
         [2,4,0,4,8],\
         [4,2,4,0,8],\
         [8,8,8,8,8]]
def globalAlign(x,y):
    G= []
    for i in range(len(x)+1): #fill up G with zero
        G.append([0]*(len(y)+1))
    
    ##reinitialise first row and first column with penalty 
    for i in range(1,len(x)+1):
        G[i][0] = G[i-1][0] + score[alphabet.index(x[i-1])][-1] #penalises all skips for empty y, alpha.index tells us what base to look for in  penalty matrix
    for i in range(1, len(y)+1):
        G[0][i] = G[0][i-1] + score[-1][alphabet.index(y[i-1])] # penalise all skips for empty x
    
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = G[i][j-1] + score[-1][alphabet.index(y[j-1])] #penalty for skipping current base in y
            distVer = G[i-1][j] + score[alphabet.index(x[i-1])][-1]  # penalty for skipping char in x
            if x[i-1] == y[j-1]:
                distDiag = G[i-1][j-1]
            else:
                distDiag = G[i-1][j-1] + score[alphabet.index(x[i-1])][alphabet.index(y[j-1])] #penalty for mismatch for current char in x and y
            G[i][j] = min(distHor, distVer, distDiag)
    return G[-1][-1]   
