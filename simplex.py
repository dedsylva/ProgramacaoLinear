import numpy as np
from duasfases import DuasFases 
from bigM import BigM 

'''
Problema Infactível
A = np.array(
    [[1, 1,-1,0],
    [-1, 1,0,-1]])

b = np.array([2,1])
c = np.array([-1,2,0,0])
'''

A = np.array(
    [[0, 1, 1, 0],
     [1, 1, 0, 1]]
)

b = np.array([3, 6])
c = np.array([0, -1, 0, 0])


df = DuasFases(A,b,c)
df.resolve()

df = BigM(A,b,c)
df.resolve()