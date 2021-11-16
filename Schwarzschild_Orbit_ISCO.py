# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 17:18:46 2021

@author: GuillermoRV
"""
import numpy as np
from matplotlib import pyplot as plt
#Unidades Geometricas G=1; c=1
M=1;l=5;R=2*M
Rmas=l**2*(1+np.sqrt(1-12*M**2/l**2))/(2*M)
Rmenos=l**2*(1-np.sqrt(1-12*M**2/l**2))/(2*M)
t=np.linspace(0,50,1000)
u=np.zeros(len(t))
v=np.zeros(len(t))
x=np.zeros(len(t))
y=np.zeros(len(t))
r=np.zeros(len(t))
phi=np.linspace(min(t),max(t),len(t))
h=phi[1]-phi[0]
def f(u):
    return 3*M*u**2-u+M/l**2
u[0]=1/(30)
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
    x[i]=r[i]*np.cos(phi[i])
    y[i]=r[i]*np.sin(phi[i])
fig, ax = plt.subplots()
BH=plt.Circle((0,0),R,color='black');ax.add_patch(BH)
plt.plot(x,y,label='R='+str(1/u[0],color='orange'))
#ax.plot(x9,y9,label='R=9',color='purple')
#ax.plot(xmenos,ymenos,label='R-≈3.48',color='b')
#ax.plot(x3,y3, label='R='+str(1/u[0]),color='green')
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=True, ncol=5)
plt.xlabel('x(Geometrical Units)');plt.ylabel('y(Geometrical Units)')
plt.xlim(-35,35);plt.ylim(-35,35)
plt.title('M='+str(M)+', l='+str(l))
plt.grid()
print(l**2>12*M**2)