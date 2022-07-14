import sympy as sym
from sympy.abc import G
import numpy as np

matrix_height=6  #矩阵总高度（宽度固定为6）
target_stat=7         #强化目标
eff=2            #初始有效词条
G_length=target_stat-eff-5+1       #初始G占有的高度（也就是多少以上词条可以接受）
H=np.array([3780,520,-5165,-13650,-26880,-50315]) #经验损失矩阵

#generate the Snell envelope matrix
S=np.zeros((6,6),dtype=object)
S[...,5]=np.concatenate((np.repeat(G,G_length),np.repeat(H[5],matrix_height-G_length)))
for i in np.arange(4,-1,-1):
    for j in np.arange(matrix_height-1-i,6):
        S[j,i]=S[j-1,i+1]*eff/4+S[j,i+1]*(4-eff)/4
        if isinstance(S[j,i],float):
            S[j,i]=H[i]
print(S)

#determine the critical value
def critValue(S):
    for i in np.arange(4,-1,-1):
        for j in np.arange(matrix_height-1-i,6):
            if type(S[j,i])==sym.core.add.Add:
                S[j,i]=sym.solve(S[j,i]-H[i])[0]
            else:
                S[j,i]=0
    print(S)
    
#determine the numeric value of Snel envelope matrix(with a pre-determined G)
def numValue(G, target_num=target_stat, eff_num=eff):
    G_length_num=target_num-eff_num-5+1

    S=np.zeros((6,6),dtype=object)
    S[...,5]=np.concatenate((np.repeat(G,G_length_num),np.repeat(H[5],matrix_height-G_length)))
    for i in np.arange(4,-1,-1):
        for j in np.arange(matrix_height-1-i,6):
            S[j,i]=S[j-1,i+1]*eff/4+S[j,i+1]*(4-eff)/4
            if S[j,i]<H[i]:
                S[j,i]=H[i]
    return S


