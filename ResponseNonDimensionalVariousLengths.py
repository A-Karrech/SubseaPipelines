from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as N


D=0.034	  #m
He = 2*D
L=8*D #m
A = 3.14*D**2/4
I = 3.14*D**4/64 #m^4
r = (I/A)**0.5
E = 48000 # N/m^2
EI = E*I
ks = 1000000 # N/m
qe = 2200*A*9.81
se = E*A*0.00001*50


H = He/r
S = se*L**2/EI;
Def = H**2/2;
q = qe*L**4/(r*EI);
k = (ks*L**4/EI)**0.25
eta = (S - Def)**0.5
eta0 = S**0.5
alp = (k**2/2 - (eta0/2)**2)**0.5 ;
bet = (k**2/2 + (eta0/2)**2)**0.5 ;

a = N.array([[1, 0, 1, 0, -1, 0], [0, 1, 0, eta, alp, -bet], [0, 0, -eta**2, 0, bet**2 - alp**2, 2*alp*bet], [0, 0, 0, eta**3, 3*bet**2*alp - alp**3,3*alp**2*bet - bet**3], [1, -1, N.cos(eta), -N.sin(eta), 0, 0], [0,1, eta*N.sin(eta), eta*N.cos(eta), 0, 0]])
b = N.array([q/k**4, 0,-q/eta**2,0,-H-q/eta**2/2,q/eta**2])
coefs = N.linalg.solve(a, b)


def BeamOnElasticResponse(x):
    """
    This function calculates the response of a beam on an elastic foundation
    """
    return q/k**4 + N.exp(-alp*x)*(coefs[4]*N.cos(bet*x) +coefs[5]*N.sin(bet*x))
 
def BeamOnFreeSpanResponse(x):
    """
    This function calculates the response of a free spanning beam
    """
            
    return q*x**2/(2*eta**2) + coefs[0] + coefs[1]*x + coefs[2]*N.cos(eta*x) + coefs[3]*N.sin(eta*x)
    
 
# np.linspace(start = 0, stop = 100, num = 5)
x = np.linspace(0,1,50)
y = BeamOnElasticResponse(x)
plt.plot(x,y, 'r',marker='*')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,y1, 'b', marker='*', label='E = 48 kPa')
plt.legend(loc='lower right') #plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('w')
 
 
E = 480000 # N/m^2
EI = E*I
H = He/r
S = se*L**2/EI;
Def = H**2/2;
q = qe*L**4/(r*EI);
k = (ks*L**4/EI)**0.25
eta = (S - Def)**0.5
eta0 = S**0.5
alp = (k**2/2 - (eta0/2)**2)**0.5 ;
bet = (k**2/2 + (eta0/2)**2)**0.5 ;

a = N.array([[1, 0, 1, 0, -1, 0], [0, 1, 0, eta, alp, -bet], [0, 0, -eta**2, 0, bet**2 - alp**2, 2*alp*bet], [0, 0, 0, eta**3, 3*bet**2*alp - alp**3,3*alp**2*bet - bet**3], [1, -1, N.cos(eta), -N.sin(eta), 0, 0], [0,1, eta*N.sin(eta), eta*N.cos(eta), 0, 0]])
b = N.array([q/k**4, 0,-q/eta**2,0,-H-q/eta**2/2,q/eta**2])
coefs = N.linalg.solve(a, b)
 

x = np.linspace(0,1,50)
y = BeamOnElasticResponse(x)
plt.plot(x,y, 'r',marker='.')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,y1, 'b', marker='.', label='E = 480 kPa')
plt.legend(loc='lower right') #plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('w')



E = 4800000 # N/m^2
EI = E*I
H = He/r
S = se*L**2/EI;
Def = H**2/2;
q = qe*L**4/(r*EI);
k = (ks*L**4/EI)**0.25
eta = (S - Def)**0.5
eta0 = S**0.5
alp = (k**2/2 - (eta0/2)**2)**0.5 ;
bet = (k**2/2 + (eta0/2)**2)**0.5 ;

a = N.array([[1, 0, 1, 0, -1, 0], [0, 1, 0, eta, alp, -bet], [0, 0, -eta**2, 0, bet**2 - alp**2, 2*alp*bet], [0, 0, 0, eta**3, 3*bet**2*alp - alp**3,3*alp**2*bet - bet**3], [1, -1, N.cos(eta), -N.sin(eta), 0, 0], [0,1, eta*N.sin(eta), eta*N.cos(eta), 0, 0]])
b = N.array([q/k**4, 0,-q/eta**2,0,-H-q/eta**2/2,q/eta**2])
coefs = N.linalg.solve(a, b)
 

x = np.linspace(0,1,50)
y = BeamOnElasticResponse(x)
plt.plot(x,y, 'r',marker='1')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,y1, 'b', marker='1', label='E = 4.8 MPa')
plt.legend(loc='lower right') #plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('w')
 
 
 
E = 48000000 # N/m^2
EI = E*I
H = He/r
S = se*L**2/EI;
Def = H**2/2;
q = qe*L**4/(r*EI);
k = (ks*L**4/EI)**0.25
eta = (S - Def)**0.5
eta0 = S**0.5
alp = (k**2/2 - (eta0/2)**2)**0.5 ;
bet = (k**2/2 + (eta0/2)**2)**0.5 ;

a = N.array([[1, 0, 1, 0, -1, 0], [0, 1, 0, eta, alp, -bet], [0, 0, -eta**2, 0, bet**2 - alp**2, 2*alp*bet], [0, 0, 0, eta**3, 3*bet**2*alp - alp**3,3*alp**2*bet - bet**3], [1, -1, N.cos(eta), -N.sin(eta), 0, 0], [0,1, eta*N.sin(eta), eta*N.cos(eta), 0, 0]])
b = N.array([q/k**4, 0,-q/eta**2,0,-H-q/eta**2/2,q/eta**2])
coefs = N.linalg.solve(a, b)
 

x = np.linspace(0,1,50)
y = BeamOnElasticResponse(x)
plt.plot(x,y, 'r')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,y1, 'b', label='E = 48 MPa')
plt.legend(loc='lower right') #plt.legend(loc='upper left')
plt.xlabel('x')
plt.ylabel('w')

plt.show()
 
 
 
'''
fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Deflection per washer (mm)')
ax1.set_ylabel('Force (kN)', color=color)
ax1.plot(x, y/1000, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Energy (kJ)', color=color)  # we already handled the x-label with ax1
ax2.plot(x, z*21, color=color)
ax2.tick_params(axis='y', labelcolor=color)

#both the xy (arrow tip) and xytext locations (text location) are in data coordinates. There are a variety of other coordinate systems one can choose -- you can specify the coordinate system of xy and xytext with one of the following strings for xycoords and textcoords (default is 'data')

ax1.annotate('Fully installed height', xy=(6.15, 0),  xycoords='data',
            xytext=(0.5, 0.5), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05, width = .5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            )
            
ax1.annotate('Valve open', xy=(12, 0),  xycoords='data',
            xytext=(0.84, 0.5), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', shrink=0.05, width = .5, headwidth = 5),
            horizontalalignment='right', verticalalignment='top',
            )

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


'''
