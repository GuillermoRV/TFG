# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 17:18:46 2021

@author: GuillermoRV
"""
import numpy as np
from matplotlib import pyplot as plt
M=1.98*10**30;l=9.1*10**38;w0=1;m=3.285*10**23;G=6.67*10**-11;e=0.205;mu=-G*M;p=l**2/(G*(m+M)*m**2)
rmax=p/(1-e)
rmin=p/(1+e)
t=np.linspace(0,100,1000)
x=np.zeros(len(t))
y=np.zeros(len(t))
u=np.zeros(len(t))
v=np.zeros(len(t))
r=np.zeros(len(t))
R=np.zeros(len(t))
phi=np.linspace(min(t),w0*max(t),len(t))
h=phi[1]-phi[0]
def f(u):
    return -u+G*M*m**2/l**2
u[0]=1/rmax
#for i in range(len(t)):
#    R[i]=l**2/(G*M*m**2)*(1/(1+e*np.cos(phi[i])))
for i in range(len(t)):
    if i!=0:
        k1=h*v[i-1]
        l1=h*f(u[i-1])
    
        k2=h*(v[i-1]+l1/2)
        l2=h*f(u[i-1]+k1/2)
    
        k3=h*(v[i-1]+l2/2)
        l3=h*f(u[i-1]+k2/2)
    
        k4=h*(v[i-1]+l3)
        l4=h*f(u[i-1]+k3)
    
        v[i]=v[i-1]+(l1+2*l2+2*l3+l4)/6
        u[i]=u[i-1]+(k1+2*k2+2*k3+k4)/6
for i in range(len(t)):
    r[i]=1/u[i]
for i in range(len(t)):
    x[i]=r[i]*np.cos(phi[i]);y[i]=r[i]*np.sin(phi[i])
    #x[i]=R[i]*np.cos(phi[i]);y[i]=R[i]*np.sin(phi[i])
plt.plot(x,y,label='r\N{SUBSCRIPT ZERO}='+str(1/u[0]),color='blue')
plt.title('M='+str(M)+'kg ,l='+str(l)+'kg m2/s ,m='+str(m)+'kg ,eccentricity='+str(e)+' ;(SI units)')
plt.xlabel('x(m)');plt.ylabel('y(m)')
plt.xlim(-6*10**10,8*10**10);plt.ylim(-7*10**10,7*10**10)
plt.scatter(0,0)
plt.grid()