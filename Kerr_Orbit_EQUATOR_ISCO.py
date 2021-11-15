# -*- coding: utf-8 -*-
"""
Created on Sun Sep 19 17:18:46 2021

@author: GuillermoRV
"""
import numpy as np
from matplotlib import pyplot as plt
#Unidades Geometricas G=1; c=1
M=1;l=6;a=M*0.9;e=0.999;E=(e**2-1)/2
Rmas=M+np.sqrt(M**2-a**2);Rmenos=M-np.sqrt(M**2-a**2)
Risco=(l**2+a**2*(1-e**2)+np.sqrt((l**2+a**2*(1-e**2))**2-12*M**2*(l-a*e)**2))/(2*M)
#Risco2=(l**2+a**2*(1-e**2)-np.sqrt((l**2+a**2*(1-e**2))**2-12*M**2*(l-a*e)**2))/(2*M)

t=np.linspace(0,10,1000)
u=np.zeros(len(t))
v=np.zeros(len(t))
x=np.zeros(len(t))
y=np.zeros(len(t))
r=np.zeros(len(t))
phi=np.linspace(min(t),max(t),len(t))
h=phi[1]-phi[0]
def B(u):
    return M*(7*l**3-6*a*e*l**2+a**2*l*(11-3*e**2)+2*a**3*e*(e**2-1))+4*M**3*(l-a*e)+u*(-3*a**2*l**3+3*a**4*l*(e**2-1)+2*M**2*(l-a*e)*(-8*l**2+a*(7*l*e+a*(e**2-5)))-M*u*(l-a*e)*(a**2*(-11*l**2+a*(7*l*e+4*a*(e**2-1)))-12*M**2*(l-a*e)**2+10*a**2*M*u*(l-a*e)**2))

def f(u):
    return (1-2*M*u+a**2*u**2)*(M*(l-2*a*e*(e**2-1))+u*(l*(-l**2+3*a**2*(e**2-1))-2*M**2*(2*l+2*a*e)+u*B(u)))/(l-2*M*u*(l-a*e))**3
u[0]=1/(Risco);v[0]=0
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

r=1/u
for i in range(len(t)):
    x[i]=r[i]*np.cos(phi[i])
    y[i]=r[i]*np.sin(phi[i])

#fig, ax = plt.subplots()
Ry=np.zeros(len(t));Rz=np.zeros(len(t))
for i in range(len(t)):
    Ry[i]=Rmas*np.sin(t[i])
    Rz[i]=Rmas*np.cos(t[i])
plt.plot(Ry,Rz,'black')
ax.plot(x50,y50,label='R=50',color='purple')
ax.plot(xisco,yisco,label='R≈33.68',color='b')
ax.plot(x30,y30, label='R='+str(20),color='green')
ax.legend()
plt.plot(x,y);plt.xlabel('x');plt.ylabel('y')
plt.xlim(-45,55);plt.ylim(-40,55)
plt.title('M='+str(M)+'; l='+str(l)+'; a='+str(a)+'; (Geometrical Units)')
#plt.scatter(0,0)
plt.grid()
print(l**2>12*M**2)
