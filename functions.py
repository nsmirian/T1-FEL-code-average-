#***********************************************************
def average(sai):
    from input import ne, h_max
    import numpy as np

    aver=np.zeros((2,h_max))

    for hn in range(h_max):
        aver[0,hn]=np.mean(np.cos(sai*(hn+1)))
        aver[1,hn]=np.mean(np.sin(sai*(hn+1)))
    return aver
#-***************************************************
def faw(zprime, m, zt, a):
	#if (zprime < zt) then
    fw = a
	#else
	 # faw = a - m * (zprime - zt)
    return fw
#*********************************************
def pr(zprime, m, zt, a, rho):
	#if (zprime < zt) then
    pr = 0.0
	#else
	  #pr = -(1 - fb(zprime, m, zt, a))/rho
	#end if
    return pr


#-***************************************************
def fb(zprime, m, zt, a):
	#if (zprime < zt) then
    fbb =1.0
	#else
	  #fb = 1 - (m/a) * (zprime - zt)
	#endif
    return fbb

#===============================================================================================================

def bsself(v, m, x):
    vf=1.
    for i in range(v):
        vf=vf*(i+1)
    A=1./vf
    besslf=0
    for j in range(m+2):
        i=j+1
        besslf=A*((x/2)**(v+2*(i-1)))+besslf
        A=A*(-1)/((v+i)*i)
    return besslf
