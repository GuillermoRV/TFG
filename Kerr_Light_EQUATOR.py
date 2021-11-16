# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 17:18:46 2021

@author: GuillermoRV
"""
import numpy as np
from matplotlib import pyplot as plt
#Unidades Geometricas G=1; c=1
M=1;a=M*0.9;l=6;e=0.999
b=abs(l/e)
Rmas=M+np.sqrt(M**2-a**2)
Rmenos=M-np.sqrt(M**2-a**2)
x0=1000;y0=0
phi_max=1000
u=np.zeros(phi_max)
v=np.zeros(phi_max)
x=np.zeros(phi_max)
y=np.zeros(phi_max)
r=np.zeros(phi_max)
phi=np.linspace(np.arctan(y0/x0),3.7,phi_max)
h=phi[1]-phi[0]

def f(u):
    return 3*M*(1-np.sign(l)*a/b)**2*u**2+u*((a/b)**2-1)
u[0]=1/(np.sqrt(x0**2+y0**2));v[0]=0.18
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

tau=np.linspace(0,10,phi_max)
Ry=np.zeros(phi_max);Rz=np.zeros(phi_max)
for i in range(len(tau)):
    Ry[i]=Rmas*np.sin(tau[i])
    Rz[i]=Rmas*np.cos(tau[i])
#plt.plot(Ry,Rz,'black')
fig, ax = plt.subplots()
BH=plt.Circle((0,0),Rmas,color='black');ax.add_patch(BH)
plt.plot(x,y,'b');plt.xlabel("x(Geometrical Units)");plt.ylabel("y(Geometrical Units)")
#ax.plot(x,y,label='a='+str(a),color='blue');ax.legend();plt.xlabel('x');plt.ylabel('y')
plt.grid()
plt.xlim(-28,27);plt.ylim(-7,17)
plt.title('M='+str(M)+'; l='+str(l)+'; a='+str(a)+'; e='+str(e))