#-------------------------------------------------------------------------------#
# 1392/06/16: NOTapering wiggler one-beam FEL					#
#-------------------------------------------------------------------------------#
# N S MIRIAN 		8TH SEP 2013		Simple code			#
#-------------------------------------------------------------------------------#
# dela1r : real part of field amplitude						#
# dela1i : imaginary part of field amplitude					#
# dela3r : real part of field amplitude						#
# dela3i : imaginary part of field amplitude					#
# ave1 : average used in Eq.(3)							#
# ave2 : average used in Eq.(3)							#
# sai : ponderomotive phase							#
# p  : conjugate momentun							#
# ne : number of particles							#
# delta : detuning parameter							#
# bunching_1 : bunching parameter for lower  energy beam			#
# bunching_3 : bunching parameter for lower  energy beam			#
# bunching   : bunching parameter for higher energy beam			#
# 		m=0								#
#-------------------------------------------------------------------------------#
from input import *
import numpy as np
from vectors import *
from functions import *
from mpmath import *
import scipy.special as sp
pi=np.pi
#initial valuecd
for i in range(2):   # real and imaginery
    for hn in range(h_max):   # harmonices
        dela_old[i,hn]=0.0001/(0.01**(hn)) #0.01_8

k=lambdau1/lambdar
omega=k
for j in range(ne):
   sai_old[j]=2.0*pi*(j-1.0)/(ne)
   p_old[j]=0.0
def faw(zprime, m, zt, a):
	#if (zprime < zt) then
    fw = a
	#else
	 # faw = a - m * (zprime - zt)
    return fw
#Runge-Kutta
#Yn(i+1)=Yn(i)+h/6*(Kn1+2Kn2+2Kn3+Kn4)
z=0.0

aw=e*B1*lambdau1/(2.0*np.sqrt(2.)*pi*m_e*c**2.0)
wp=(1./r)* np.sqrt(4.*gamma0*e*I*amp_to_statamp/(c*m_e*np.sqrt(gamma0**2-1.)))
#(1/(r*k_w*c))*dsqrt(4*gamma0*e*I*amp_to_statamp/(m_e*c*dsqrt(gamma0**2-1)))
k_w=2.*pi/(lambdau1)
ro=(1./gamma0)*((aw/4.)*(wp/(c*k_w)))**(2./3.)

print("wp=", wp)
print( "ro=", ro)
for t in range(num_step):
    dela=dela_old
    sai=sai_old
    p= p_old
    z = z + h
    kissi = (faw(z, m, zt, aw) ** 2.) /(2. * (1. + faw(z, m, zt, aw) ** 2))
    for hn in range(h_max):
        f[hn]=bsself(0,ns,(hn+1)*kissi)-\
                bsself(1,ns,(hn+1)*kissi)

#-----------------------------------------------------------------------------
# Kn1=m
    ave= average(sai)
    for i in range(2):   # real and imaginery
       for hn in range(h_max):   # harmonices
            m[i, hn]= h * ( ave[i,hn]*f[hn] )*(-1)**i



    msai=h*( p - pr(z, m, zt, aw, rho1)) / fb(z, m, zt, aw)
    mp=h*(-2.0*dela[0,0]*np.cos(sai)+2.0*dela[1,0]*np.sin(sai))*f[0]# +\
	         #h*(-2.*dela3r*np.cos(3.*sai(j))+2.*dela3i*np.sin(3.*sai(j)))*f3

    sai  = sai_old   + 0.5 * msai
    p   = p_old + 0.5 * mp

    dela   = dela_old   + 0.5 * m

#-----------------------------------------------------------------------------
#Kn2=g
    ave= average(sai)
    for i in range(2):   # real and imaginery
       for hn in range(h_max):   # harmonices
            g[i, hn]= h * ( ave[i,hn]*f[hn] )*(-1)**i



    gsai=h*( p - pr(z, m, zt, aw, rho1)) / fb(z, m, zt, aw)
    gp=h*(-2.0*dela[0,0]*np.cos(sai)+2.0*dela[1,0]*np.sin(sai))*f[0]# +\
	         #h*(-2.*dela3r*np.cos(3.*sai(j))+2.*dela3i*np.sin(3.*sai(j)))*f3

    sai  = sai_old   + 0.5 * gsai
    p   = p_old + 0.5 * gp

    dela   = dela_old   + 0.5 * g

##-----------------------------------------------------------------------------
#Kn3=q
    ave= average(sai)
    for i in range(2):   # real and imaginery
       for hn in range(h_max):   # harmonices
            q[i, hn]= h * ( ave[i,hn]*f[hn] )*(-1)**i



    qsai=h*( p - pr(z, m, zt, aw, rho1)) / fb(z, m, zt, aw)
    qp=h*(-2.0*dela[0,0]*np.cos(sai)+2.0*dela[1,0]*np.sin(sai))*f[0]# +\
	         #h*(-2.*dela3r*np.cos(3.*sai(j))+2.*dela3i*np.sin(3.*sai(j)))*f3

    sai  = sai_old   +  qsai
    p   = p_old +  qp

    dela   = dela_old   +  q
##-----------------------------------------------------------------------------
#Kn4=s
    ave= average(sai)
    for i in range(2):   # real and imaginery
       for hn in range(h_max):   # harmonices
            s[i, hn]= h * ( ave[i,hn]*f[hn] )*(-1)**i



    ssai=h*( p - pr(z, m, zt, aw, rho1)) / fb(z, m, zt, aw)
    sp=h*(-2.0*dela[0,0]*np.cos(sai)+2.0*dela[1,0]*np.sin(sai))*f[0]# +\
	         #h*(-2.*dela3r*np.cos(3.*sai(j))+2.*dela3i*np.sin(3.*sai(j)))*f3
#
#-----------------------------------------------------------------------------
#Kn1=mn
#Kn2=pn
#Kn3=qn
#Kn4=sn
    sai  = sai_old   + (msai+2*gsai+2*qsai+ssai)/6.
    p   = p_old +  (mp+2*gp+2*qp+sp)/6.

    dela   = dela_old   + (m+2*g+2*q+s)/6.

    dela_old= dela
    sai_old=sai
    p_old=p

    print(t)
