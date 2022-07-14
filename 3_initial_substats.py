import sympy as sym
from sympy.abc import G
import numpy as np

#now only applicable for problem with following setting
matrix_height=6  #矩阵总高度（宽度固定为6）
target_stat=5    #强化目标
eff=1            #初始有效词条
G_length=target_stat-eff-5+1       #初始G占有的高度（也就是多少以上词条可以接受）
H=np.array([3780,520,-5165,-13650,-26880,-50315]) #经验损失矩阵


'''
generate the Snell envelope matrix S(3-dimensional)
'''
S=np.zeros((6,6,6),dtype=object)

#S[t,i,j] stand for the Snell envelope at time t 
#with i effetive substat 1 and j effetive substat 2
#determine the terminal values(t=5)
for i in np.arange(1,6):      
    for j in np.arange(0,6):
        if i+j<target_stat:     #unaccepted result
            S[5,i,j]=H[5]
        elif i+j>eff+5:         #impossible occassions
            S[5,i,j]=0
        else:
            S[5,i,j]=G

#backward iterartion(until t=1)
for t in np.arange(4,0,-1):
    for i in np.arange(1,t+1):
        for j in np.arange(0,t+1):
            if i+j>eff+t:       #impossible occassions
                S[t,i,j]=0
            else:
                if j==0:        #only one effective substat
                    S[t,i,j]=S[t+1,i+1,0]/4+S[t+1,i,0]*3/4
                else:           #two effective substat
                    S[t,i,j]=S[t+1,i+1,j]/4+S[t+1,i,j+1]/4+S[t+1,i,j]/2
                if not isinstance(S[t,i,j],sym.core.add.Add):       #entries without G
                    S[t,i,j]=H[t]

#backward iteration(t=0)
#direct calculation
S[0,1,0]=S[1,1,1]*75/(100*4+150*2+75)+S[1,1,0]*(100*4+150*2)/(100*4+150*2+75)

'''
determine critical value
'''
def critValue(S):
    c=np.array([])
    for t in np.arange(0,6):
        for i in np.arange(0,6):
            for j in np.arange(0,6):
                if type(S[t,i,j])==sym.core.add.Add:
                    c=np.append(c,(sym.solve(S[t,i,j])[0],(t,i,j)))
    for i in c:
        print(i)
        print("\n")
