#!/usr/bin/env python3

# Generates multiple figures showing the region of absolute
# stability of foward Euler, or its real part [-2,0], overlain with
# z = k lam, for eigenvalues lam of the FD MOL approximation
#   U_j'(t) = (U_{j-1}(t) - 2 U_j(t) + U_{j+1}(t)) / h^2
# of the heat equation  u_t = u_xx,  on the spatial interval [0,1]
# and with time step k.  Note this is the MOL system
#   U'(t) = A^h U(t)
# where A^h is the matrix in equation (2.10) of LeVeque (2007).
# The eigenvalues lam are from equation (2.23).

import numpy as np
import matplotlib.pyplot as plt

def writeout(outname):
    print('writing file ' + outname)
    plt.savefig(outname,bbox_inches='tight')

def zplane(box):
    x, y = np.linspace(box[0],box[1],250), np.linspace(box[2],box[3],250)
    xx, yy = np.meshgrid(x,y)
    zz = xx + yy * 1j
    fig = plt.figure()
    plt.contourf(x,y,abs(zz + 1.0),[0.0, 1.0],colors='C1',alpha=0.25)
    plt.xlabel('Re(z)')
    plt.ylabel('Im(z)')
    plt.grid(True)
    plt.plot(box[0:2],[0.0, 0.0],'C0')
    plt.plot([0.0, 0.0], box[2:],'C0')
    plt.axis('scaled')
    plt.axis(box)

# z-plane figures
for m in [3, 4, 9, 19, 39]:
    h = 1.0 / (m + 1)
    p = np.arange(m) + 1   # = 1, ..., m
    lam = (2.0/h**2) * (np.cos(p*np.pi*h) - 1.0)
    k = 2.0 / (-lam[-1])
    z = k * lam
    zplane([-2.5, 1.0, -1.5, 1.5])
    for j in range(m):
        plt.plot(z[j],0.0,'o',color='C3',markersize=5)
        if m < 10 or j == m-1:
            plt.text(z[j],0.12,r'$z_{%d}$' % (j+1)) 
    plt.title(r'$m=%d$ and $h=%.3f$: stability requires $k$ <= $%.4f$' \
              % (m, h, k))
    plt.show()

# Re(z) values depend on h figure
fig = plt.figure()
plt.xlabel(r'Re(z)')
plt.ylabel(r'$h$')
plt.grid(True)
box = [-2.5, 1.0, -0.02, 0.3]
plt.plot(box[0:2],[0.0, 0.0],'C0')
plt.plot([0.0, 0.0], box[2:],'C0')
plt.plot([-2.0, 0.0],[0.0, 0.0],linewidth=5.0,color='C1',alpha=0.45)
m = np.array([3, 4, 9, 19, 39])
h = 1.0 / (m + 1)
lamleft = (2.0/h**2) * (np.cos(m*np.pi*h) - 1.0)
k = 2.0 / (-lamleft)
for p in range(1,39+1):
    lamp = (2.0/h**2) * (np.cos(p*np.pi*h) - 1.0)
    zp = k * lamp
    plt.plot(zp[p<=m],h[p<=m],'-o',color='C3',markersize=4, linewidth=0.5)
    if p <= 3:
        plt.text(zp[0],h[0]+0.01,r'$z_{%d}$' % p)
for j in range(1,len(h)):
    plt.text(-2.0-0.2,h[j]+0.005,r'$z_{%d}$' % (m[j]))
plt.show()

# refinement path figure
fig = plt.figure()
plt.xlabel(r'$k$')
plt.ylabel(r'$h$')
plt.grid(True)
box = [-0.005, 0.05, -0.02, 0.3]
plt.plot(box[0:2],[0.0, 0.0],'C0')
plt.plot([0.0, 0.0], box[2:],'C0')
plt.axis(box)
m = np.array([3, 4, 9, 19, 39],dtype=float)
h = 1.0 / (m + 1)
lamleft = (2.0/h**2) * (np.cos(m*np.pi*h) - 1.0)
k = 2.0 / (-lamleft)
plt.plot(k,h,'--o',color='C5',markersize=8)
plt.plot('refinement path'
plt.show()

