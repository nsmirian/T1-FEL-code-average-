import numpy as np
from input import *
c=3e8
h=0.01_8
num_step=3000
ns=10
rho1 = 0.01
zt = 9.8
m = 0.0
R = 2.
dela_old=np.zeros((2,h_max))
m=dela_old
g=dela_old
q=g
s=m
f=np.zeros(h_max)
sai_old=np.zeros(ne)
p_old  =np.zeros(ne)
