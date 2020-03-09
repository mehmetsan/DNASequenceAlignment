import sys
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 22:23:24 2019

@author: MehmetSanisoglu
"""
match = 2
mismatch = -3
def writeToFile(fileName, seq1, seq2):
    out = open(fileName,"w")
    maxChars = max(len(seq1),len(seq2))
    sections = (maxChars / 60) + 1
    firstText = "my_first_sequence\t\t"
    secondText = "another_sequence\t\t"
    print(maxChars, sections, len(seq1),len(seq2))
    for i in range(sections):
        firstContent = ""
        secondContent = ""
        for j in range(60):
            if(maxChars == (i*60)+j):
                break
            firstContent += seq1[(i*60)+j]
        out.write(firstText + firstContent + "\n" )
        for j in range(60):
            if(maxChars == (i*60)+j):
                break
            secondContent += seq2[(i*60)+j]
        out.write(secondText + secondContent + "\n\n" ) 
    
    out.close()
    return
    
def retrieveSequences(filename):
    seq1, seq2 = "", ""
    with open(filename, 'r') as file:
        first = False
        for each in file:
            if (">" not in each ):
                    if(first):
                        if("\n" in each):
                            seq1 += each[0:-1]
                        else:
                            seq1 += each
                    else:
                        if("\n" in each):
                            seq2 += each[0:-1]
                        else:
                            seq2 += each
            else:
                first = not first
    seq1 = seq1.replace("\r","")
    seq2 = seq2.replace("\r","")
    file.close()
    return seq1, seq2


def localBacktrack(matrix, S1, S2 ):
    M = matrix
               
    constructedSeq1 = []
    constructedSeq2 = []  
    len1 = len(M[0])
    len2 = len(M)
         
    
    indexX = 1
    indexY = 1
    currentMax = float("-inf")

    for y in range(1, len2):
        for x in range(1, len1):
            if M[y][x][0] > currentMax:
                currentMax = M[y][x][0]
                indexY = y
                indexX = x
                
    y = indexY
    x = indexX

    while(x != 0 or y !=0):
        if(M[y][x][1] == "d"):
            constructedSeq1.insert(0,S1[x-1]) 
            constructedSeq2.insert(0,S2[y-1]) 
            x -= 1
            y -= 1
        elif(M[y][x][1] == "l"):
            constructedSeq1.insert(0,S1[x-1]) 
            constructedSeq2.insert(0,"_")
            x -= 1
        else:
            constructedSeq1.insert(0,"_")
            constructedSeq2.insert(0,S2[y-1]) 
            y -= 1
    return [constructedSeq1,constructedSeq2], currentMax
          
def globalBacktrack(matrix, S1, S2 ):
    M = matrix
               
    constructedSeq1 = []
    constructedSeq2 = []
        
    y = len(matrix)-1
    x = len(matrix[0])-1
     
    score = M[y][x][0]
    while(x != 0 or y !=0):
        if(M[y][x][1] == "d"):
            constructedSeq1.insert(0,S1[x-1]) 
            constructedSeq2.insert(0,S2[y-1]) 
            x -= 1
            y -= 1
        elif(M[y][x][1] == "l"):
            constructedSeq1.insert(0,S1[x-1]) 
            constructedSeq2.insert(0,"_")
            x -= 1
        else:
            constructedSeq1.insert(0,"_")
            constructedSeq2.insert(0,S2[y-1]) 
            y -= 1
    return [constructedSeq1,constructedSeq2], score   
 
    
def globalMethod(gapopen, s1, s2):
    gapPenalty = gapopen

    seq1 = s1
    seq2 = s2
    
    len1 = len(seq1) 
    len2 = len(seq2)        
    Matrix = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    
    #CREATIN THE BORDERS OF THE TABLE
    for i in range(len1):
        Matrix[0][i+1] = (gapPenalty * (i+1),"l")
    for j in range(len2):   
        Matrix[j+1][0] = (gapPenalty * (j+1),"u")
    Matrix[0][0] = (0,"")
    
    #FILLING THE INSIDE OF THE TABLE    
    for y in range(1,len2+1):
        for x in range(1,len1+1):
            diag = 0
            if(seq1[x-1] == seq2[y-1]):
                diag = Matrix[y-1][x-1][0] + 2
            else:
                diag = Matrix[y-1][x-1][0] - 3
            
            up = Matrix[y-1][x][0] + gapPenalty
            left = Matrix[y][x-1][0] + gapPenalty
            
            result = max(diag,up,left)
            #USING TUPLES TO STORE THE DIRECTION TO FOLLOW DURING BACKTRACKING
            if(result == diag):
                temp = (result,"d")
                Matrix[y].insert(x,temp)
            elif(result == up):
                temp = (result,"u")
                Matrix[y].insert(x,temp)
            else:
                temp = (result,"l")
                Matrix[y].insert(x,temp)
            #KEEPING THE SIZE SAME
            Matrix[y].pop(len(Matrix[y])-1)
    seqs, score = globalBacktrack(Matrix, seq1, seq2)
    writeToFile("global-naiveGap.aln",seqs[0],seqs[1])
    return score


#############################################
def aglobalMethod(gapopen, gapextend,s1, s2):
    Wg, We = gapopen, gapextend
    
    seq1 = s1
    seq2 = s2
    
    #INITIALIZATION
    len1 = len(seq1) 
    len2 = len(seq2) 
    E = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    F = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    G = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    V = [[0 for x in range(len1+1)] for y in range(len2+1)] 
      
    for j in range(len1+1):  
        E[0][j] = 0
        F[0][j] = Wg + j * We
        G[0][j] = 0
        V[0][j] = (Wg + j * We,"l")
    for i in range(len2+1):
        E[i][0] = Wg + i * We
        F[i][0] = 0
        G[i][0] = 0
        V[i][0] = (Wg + i * We,"u")
    E[0][0] = float('-inf')
    F[0][0] = float('-inf')
    G[0][0] = float('-inf')
    V[0][0] = (0,"l")
    
    for j in range(1,len1+1):
        for i in range(1,len2+1):
            E[i][j] = max( E[i][j-1]+We, V[i][j-1][0]+Wg+We)
            F[i][j] = max( F[i-1][j]+We, V[i-1][j][0]+Wg+We )          
            if(seq1[j-1] == seq2[i-1]):
                G[i][j] = V[i-1][j-1][0] + match
            else:
                G[i][j] = V[i-1][j-1][0] + mismatch   
            result = max(G[i][j], E[i][j], F[i][j])
            if(result == E[i][j]):
                temp = (result,"l")
                V[i].insert(j,temp)
            elif(result == F[i][j]):
                temp = (result,"u")
                V[i].insert(j,temp)
            else:
                temp = (result,"d")
                V[i].insert(j,temp)
            V[i].pop(len(V[i])-1)    
    
    seqs, score = globalBacktrack(V, seq1, seq2 )
    writeToFile("global-affineGap.aln",seqs[0],seqs[1])
    return score
#############################################
    
def localMethod(gapopen, s1, s2):
    gapPenalty = gapopen

    seq1 = s1
    seq2 = s2
    
    len1 = len(seq1) 
    len2 = len(seq2)        
    Matrix = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    
    #CREATIN THE BORDERS OF THE TABLE
    for i in range(len1):
        Matrix[0][i+1] = (gapPenalty * (i+1),"l")
    for j in range(len2):   
        Matrix[j+1][0] = (gapPenalty * (j+1),"u")
    Matrix[0][0] = (0,"")
    
    #FILLING THE INSIDE OF THE TABLE    
    for y in range(1,len2+1):
        for x in range(1,len1+1):
            diag = 0
            if(seq1[x-1] == seq2[y-1]):
                diag = Matrix[y-1][x-1][0] + 2
            else:
                diag = Matrix[y-1][x-1][0] - 3
            
            up = Matrix[y-1][x][0] + gapPenalty
            left = Matrix[y][x-1][0] + gapPenalty
            
            result = max(0,diag,up,left)
            #USING TUPLES TO STORE THE DIRECTION TO FOLLOW DURING BACKTRACKING
            if(result == diag):
                temp = (result,"d")
                Matrix[y].insert(x,temp)
            elif(result == up):
                temp = (result,"u")
                Matrix[y].insert(x,temp)
            else:
                temp = (result,"l")
                Matrix[y].insert(x,temp)
            #KEEPING THE SIZE SAME
            Matrix[y].pop(len(Matrix[y])-1)
    seqs, score = localBacktrack(Matrix, seq1, seq2)
    writeToFile("local-naiveGap.aln",seqs[0],seqs[1])
    return score

#############################################

def alocalMethod(gapopen, gapextend, s1, s2):
    
    Wg, We = gapopen, gapextend
    seq1 = s1
    seq2 = s2
    
    #INITIALIZATION
    len1 = len(seq1) 
    len2 = len(seq2) 
    E = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    F = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    G = [[0 for x in range(len1+1)] for y in range(len2+1)] 
    V = [[0 for x in range(len1+1)] for y in range(len2+1)] 
      
    for j in range(len1+1):  
        E[0][j] = 0
        F[0][j] = Wg + j * We
        G[0][j] = 0
        V[0][j] = (Wg + j * We,"l")
    for i in range(len2+1):
        E[i][0] = Wg + i * We
        F[i][0] = 0
        G[i][0] = 0
        V[i][0] = (Wg + i * We,"u")
    E[0][0] = float('-inf')
    F[0][0] = float('-inf')
    G[0][0] = float('-inf')
    V[0][0] = (0,"l")
    
    for j in range(1,len1+1):
        for i in range(1,len2+1):
            E[i][j] = max( E[i][j-1]+We, V[i][j-1][0]+Wg+We)
            F[i][j] = max( F[i-1][j]+We, V[i-1][j][0]+Wg+We )          
            if(seq1[j-1] == seq2[i-1]):
                G[i][j] = V[i-1][j-1][0] + match
            else:
                G[i][j] = V[i-1][j-1][0] + mismatch   
            result = max(0,G[i][j], E[i][j], F[i][j])
                       
            
            if(result == G[i][j]):
                temp = (result,"d")
                V[i].insert(j,temp)
            elif(result == F[i][j]):
                temp = (result,"u")
                V[i].insert(j,temp)
            else:
                temp = (result,"l")
                V[i].insert(j,temp)
            V[i].pop(len(V[i])-1)    

    seqs, score = globalBacktrack(V, seq1, seq2 )
    writeToFile("local-affineGap.aln",seqs[0],seqs[1])
    return score

def main():  
    
    mode = sys.argv[2]
    inputFile = sys.argv[4]
    gapOpen = int(sys.argv[6])
    gapExt = -1     #DEFAULT VALUE
    if(len(sys.argv)>8):
        gapExt = int(sys.argv[8])  
    seq1 , seq2 = retrieveSequences(inputFile)  
    if( mode == "global"):
        print("Global Score =",globalMethod(gapOpen,seq1,seq2))
    elif( mode == "aglobal"):
        print("Aglobal Score =",aglobalMethod(gapOpen,gapExt,seq1,seq2))
    elif( mode == "local"):
        print("Local Score =",localMethod(gapOpen,seq1,seq2))
    elif( mode == "alocal"):
        print("Alocal Score =",alocalMethod(gapOpen,gapExt,seq1,seq2))
    else:
        print("Incorrect mode typed")
  
if __name__== "__main__":
    main()






       

        
        
        
        
    