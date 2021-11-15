# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 17:18:46 2021

@author: GuillermoRV
"""
import numpy as np
from matplotlib import pyplot as plt
#Unidades Geometricas G=1; c=1
M=1;R_S=2*M;R_photon=3*M
x0=1000;y0=0
phi_max=1000
u=np.zeros(phi_max)
v=np.zeros(phi_max)
x=np.zeros(phi_max)
y=np.zeros(phi_max)
r=np.zeros(phi_max)
phi=np.linspace(np.arctan(y0/x0),7.5,phi_max)
h=phi[1]-phi[0]
def f(u):
    return 3*M*u**2-u
u[0]=1/(np.sqrt(x0**2+y0**2));v[0]=0.1911
#fig, ax = plt.subplots();BH=plt.Circle((0,0),R_S,color='black');ax.add_patch(BH);plt.plot(R_photon*np.cos(phi),R_photon*np.sin(phi),'red')
#plt.grid()
plt.xlim(-28,27);plt.ylim(-7,17)
i=1
while i<phi_max:
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
    i+=1
for i in range(phi_max):
    r[i]=1/u[i]

for i in range(phi_max):
    x[i]=r[i]*np.cos(phi[i])
    y[i]=r[i]*np.sin(phi[i])

plt.plot(x,y,'b')
plt.xlabel("x")
plt.ylabel("y")
