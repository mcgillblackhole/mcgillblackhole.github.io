#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  7 06:39:55 2021

@author: vitor
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
import matplotlib.colors as colors

from matplotlib import rc
rc('font',**{'family':'serif','serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#%%
# =============================================================================
#### EIH equations of motion

np.random.seed(5)

def dist(x1, x2):
    
    d = np.linalg.norm(x1 - x2)
    
    return d

def force(m2, x1, x2):
    
    f0 = m2*(x2 - x1)/dist(x1, x2)**3
    
    return f0

def force2(m1, m2, m3, x1, x2, x3):
    
    dis = m2/dist(x1, x2)
    f2 = force(m1, x2,x1) + force(m3, x2,x3)
    
    return (7/2)*dis*f2


def force3(m2, x1, x2, v1, v2):
    
    f = m2/(dist(x1,x2)**2)*(np.dot((x1 - x2)/dist(x1,x2),(4*v1 - 3*v2)) * (v1-v2))
    
    return f

def force4(m1, m2, m3, x1, x2, x3, v1, v2):
    
    f = force(m2, x1, x2)*(np.dot(v1, v1) + 2*np.dot(v2,v2) - 4*np.dot(v1,v2) - 
                       3/2 * np.dot((x1-x2)/dist(x1,x2), v2)**2 + 
                       1/2 * np.dot(x2 - x1, force(m1, x2, x1) + force(m3, x2, x3)))
    return f

def force5(m2, x1, x2, x3):
    
    r12 = 1/dist(x1, x2)
    r23 = 1/dist(x2, x3)
    r13 = 1/dist(x1, x3)
    f1 = force(m2, x1,x2)*(r12 + r13) 
    f2 = force(m2, x1,x2)*(r12 + r23)
    
    return (-4*f1 -f2)


def dydx(t, u):
    
    """
    x1, y1, z1, vx1, vy1, vz1 = u[:6]
    x2, y2, z2, vx2, vy2, vz2 = u[6:12]
    x3, y3, z3, vx3, vy3, vz3 = u[12:]
    """
    m1, m2, m3 = 1, 1, 1
    
    x1 = u[:3]
    x2 = u[6:9]
    x3 = u[12:15]
    
    v1 = u[3:6]
    v2 = u[9:12]
    v3 = u[15:]
    
    f01 = force(m2, x1, x2) + force(m3, x1, x3)
    f11 = force2(m1, m2, m3, x1, x2, x3) + force2(m1, m3, m2, x1 ,x3, x2)
    f21 = force3(m2, x1, x2, v1, v2) + force3(m3, x1, x3, v1, v3)
    f31 = force4(m1, m2, m3, x1, x2, x3, v1, v2) + force4(m1, m3, m2, x1, x3, x2, v1, v3)
    f41 = force5(m2, x1,x2,x3) + force5(m3, x1,x3,x2)

    
    f02 = force(m1, x2,x1) + force(m3, x2, x3)
    f12 = force2(m2, m1, m3, x2, x1, x3) + force2(m2, m3, m1, x2, x3, x1)
    f22 = force3(m1, x2, x1, v2, v1) + force3(m3, x2, x3, v2, v3)
    f32 = force4(m2, m1, m3, x2, x1, x3, v2, v1) + force4(m2, m3, m1, x2, x3, x1, v2, v3)
    f42 = force5(m1, x2, x1, x3) + force5(m3, x2, x3, x1)

    
    f03 = force(m1, x3, x1) + force(m2, x3, x2)
    f13 = force2(m3, m1, m2, x3, x1, x2) + force2(m3, m2, m1, x3, x2, x1)
    f23 = force3(m1, x3, x1, v3, v1) + force3(m2, x3, x2, v3, v2)
    f33 = force4(m3, m2, m1, x3, x2, x1, v3, v2) + force4(m3, m1, m2, x3, x1, x2, v3, v1)
    f43 = force5(m1, x3,x1,x2) + force5(m2, x3,x2,x1)
    
    tol2 = 1e-3*0
    tol3 = 1e-2*0
    tol4 = 1e-2*0
    tol4 = 1e-3*0
    
    return np.array([*v1,*(f01 + tol2*f11 + tol3*f21 + tol4*f31 + tol4*f41 ),  ## 1,2 + 1,3
                     *v2,*(f02 + tol2*f12 + tol3*f22 + tol4*f32 + tol4*f42 ),  ## 2,1 + 2,3
                     *v3,*(f03 + tol2*f13 + tol3*f23 + tol4*f33 + tol4*f43 )]) ## 3,1 + 3,2

                    
L = 5
G = M = 1


# =============================================================================


#### Eight figure
Pos_0 = np.array([ [-0.97000436, 0.24308753,0],
                  [0,0,0],
                  [0.97000436, -0.24308753,0]   ])


Vel_0 = np.array([ [ 0.4662036850, 0.4323657300,0],[-0.93240737, -0.86473146,0],
                                      [ 0.4662036850, 0.4323657300,0]])


##### Circle
# Pos_0 = np.array([ [-L/2,-((np.sqrt(3)*L/2) - L/np.sqrt(3)),0.0],
#                    [0.0,(np.sqrt(3)*L/2)-((np.sqrt(3)*L/2) - L/np.sqrt(3)),1.0],
#                    [L/2,-((np.sqrt(3)*L/2) -L/np.sqrt(3)),2.0]   ])

# Vel_0 = (np.sqrt( G*M/L ) )*np.array([ [ 1/2 , -np.sqrt(3)/2 , 0.0],[-1,0.0,0.0],[ 1/2 , np.sqrt(3)/2 ,0.0]])

# N_part = 3

# Pos_0 = np.random.randn(N_part,3)
# Vel_0 = 0.5*np.random.randn(N_part,3)


u0 = np.zeros(18)

u0[:3] = Pos_0[0,:]
u0[3:6] = Vel_0[0,:]
u0[6:9] = Pos_0[1,:]
u0[9:12] = Vel_0[1,:]
u0[12:15] = Pos_0[2,:]
u0[15:] = Vel_0[2, :]

N_steps = 6000

dt = 100/N_steps

sol = solve_ivp(dydx, [0, 100], u0, t_eval = np.linspace(0, 100, N_steps ), method = 'Radau', max_step = 1e-2)

# =============================================================================

# for a given time instant idt 
x1 = sol.y[0, :]
y1 = sol.y[1, :]
z1 = sol.y[2, :]

Pos_1 = x1,y1,z1


x2 = sol.y[6, :]
y2 = sol.y[7, :]
z2 = sol.y[8, :]

Pos_2 = x2,y2,z2


x3 = sol.y[12, :]
y3 = sol.y[13, :]
z3 = sol.y[14, :]

Pos_3 = x3,y3,z3


r = 1
m1 = 1
m2 = m1
m3 = m1

# %%
# for i in range(100):
    
#     ax = plt.axes(projection = '3d')
#     ax.scatter3D(sol.y[0, i], sol.y[1, i], sol.y[2, i], color = 'C0')
#     ax.scatter3D(sol.y[6, i], sol.y[7, i], sol.y[8, i], color = 'C1')
#     ax.scatter3D(sol.y[12, i], sol.y[13, i], sol.y[14, i], color = 'C2')
#     # ax.plot(sol.y[0, max(0, i-10):i], sol.y[1, max(0, i-10):i], sol.y[2, max(0, i-10):i], color = 'C0')
#     # ax.plot(sol.y[6, max(0, i-10):i], sol.y[7, max(0, i-10):i], sol.y[8, max(0, i-10):i], color = 'C1')
#     # ax.plot(sol.y[12, max(0, i-10):i], sol.y[13, max(0, i-10):i], sol.y[14, max(0, i-10):i], color = 'C2')
#     ax.plot(sol.y[0, :i], sol.y[1, :i], sol.y[2, :i], color = 'C0')
#     ax.plot(sol.y[6, :i], sol.y[7, :i], sol.y[8, :i], color = 'C1')
#     ax.plot(sol.y[12, :i], sol.y[13, :i], sol.y[14, :i], color = 'C2')
    
    
#     ax.set_xlim3d(-5, 5)
#     ax.set_ylim3d(-5, 5)
#     ax.set_zlim3d(-5, 5)
#     plt.show()

#%%

# def h(i,j,r):
#     return -(2/r)*(m1*(Accel_1[i]*Pos_1[j] + 2*Vel_1[i]*Vel_1[j] + Pos_1[i]*Accel_1[j]) + m2*(Accel_2[i]*Pos_2[j] + 2*Vel_2[i]*Vel_2[j] + Pos_2[i]*Accel_2[j]) + m3*(Accel_3[i]*Pos_3[j] + 2*Vel_3[i]*Vel_3[j] + Pos_3[i]*Accel_3[j]) )


# =============================================================================
# def L(i,j):
#     return (m1*(Accel_1[i]*Pos_1[j] + 2*Vel_1[i]*Vel_1[j] + Pos_1[i]*Accel_1[j]) + m2*(Accel_2[i]*Pos_2[j] + 2*Vel_2[i]*Vel_2[j] + Pos_2[i]*Accel_2[j]) + m3*(Accel_3[i]*Pos_3[j] + 2*Vel_3[i]*Vel_3[j] + Pos_3[i]*Accel_3[j]) )

def L(i,j):
    return m1*np.diff(Pos_1[i]*Pos_1[j],2) + m2*np.diff(Pos_2[i]*Pos_2[j],2  ) + m3*np.diff(Pos_2[i]*Pos_2[j],2  ) 



def tr(t,r):
    return t - r

# def h_time(i,j,r,t):
#     n = [int(w) for w in t/dt]
#    # k = [int(p) for p in r/dt]
#     m = np.array(n) - int(r/dt)
#     Wave = np.zeros(len(t)+1)
#     Wave[n] = -(2/r)*L(i,j)[m]*(1 + np.sign(tr(t,r)))/2
#     return Wave


t_fixed = np.array([50])
dr = 0.1
r = np.arange(1,100-2,dr)
h_space = []



def h_space(i,j,r,t):
    #n = [int(w) for w in t/dt]
    n = int(t/dt)
    k = [int(p) for p in r/dt]
    a = [int(b) for b in r/dr]
    a = np.array(a) - a[0]
    #a = np.array(a) - 50
    
    #m = np.array(n) - np.array(k)
    m = n - np.array(k)
    Wave = np.zeros(len(r))
    Wave[a] = -(2/r)*L(i,j)[m]*(1 + np.sign(tr(t,r)))/2
    return Wave

#%%
for t_fixed in np.linspace(0, 100, int(30) ):
    
    plt.plot(r, h_space(1,0,r,t_fixed), label = '+ Polarization', color = 'black')
    # print(h_space(1,0,r,t_fixed))
    #plt.plot(r, h_space(0,0,r,t_fixed), label = 'X Polarization', color = 'red')
    plt.xlabel(r'Radial distance to the 3-body system $r$')
    plt.ylabel(r'Gravitational Wave amplitude $h_{+}(t,r)$')
    plt.ylim(-0.001,0.001)
    plt.legend()
    plt.show()   


#%%

# t_fixed = 4

# plt.plot(r, h_space(1,0,r,t_fixed), label = '+ Polarization', color = 'black')
# print(h_space(1,0,r,t_fixed))
# #plt.plot(r, h_space(0,0,r,t_fixed), label = 'X Polarization', color = 'red')
# plt.xlabel(r'Radial distance to the 3-body system $r$')
# plt.ylabel(r'Gravitational Wave amplitude $h_{+}(t,r)$')
# plt.ylim(-6,6)
# plt.legend()
# plt.show()   


#%%

### Contour Plot:

    
#### scale:
xmin = -50
xmax = 50
ymin = -50
ymax = 50

for t_fixed in np.linspace(0, 100, int(30) ):
    fig = plt.figure(figsize = (7,6))# ww w  .d e m  o2s.co m
    
    ax = fig.add_subplot()

    ax.axis([xmin, xmax, ymin ,ymax])
    

    plt.suptitle(r'Normalized amplitude of Gravitational Waves', fontsize='16')
    ax.set_xlabel(r'$ x $',fontsize=17, fontweight='bold')
    ax.set_ylabel(r'$ y $',fontsize=17, fontweight='bold')

    # r = np.arange(1,100,dr)
    p = np.linspace(0, 2*np.pi, int(900*dr/dr))
    R, P = np.meshgrid(r, p)
    X, Y = R*np.cos(P), R*np.sin(P)
    Z0 = np.zeros((900,990))
    
    
    cut = np.array(h_space(1,0,r,t_fixed))
    
    Z = 1E+4*np.tile(cut,(900,1))
    
    # Plot the surface.
    lvl = 1e+4*np.linspace(-4e-4,4e-4,50)
    # ax.contourf(X, Y, Z, zdir='z', offset=0.0, cmap='inferno', levels = lvl) #Contorno no eixo xy    # ax.plot
    # ax.contourf(X, Y, Z, zdir='z', offset=0.0, cmap='hsv', levels = lvl) #Contorno no eixo xy    # ax.plot
    
    
    plotmax = xmax
    # Z = np.log(Z)
    cs = ax.contourf(X, Y, Z, zdir='z', offset=0.0, cmap='magma', levels = lvl) #Contorno no eixo xy    # ax.plot
    fig.colorbar(cs, ax=ax, shrink=0.9)    
   
    ax.tick_params(axis='both', which='both', labelsize=14, direction="in")
    ax.minorticks_on()
    
    fig1 = plt.gcf()
    plt.show()
    plt.draw()
    fig1.savefig('%d.pdf'%(t_fixed))

    










