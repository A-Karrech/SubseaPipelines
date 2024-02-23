from pylab import *
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as N


D = 0.452      #m
thk = D/15. #DNVGL-ST-F101 (2017),wall-thickness-to-outer-diameter ratio (D t) of 15.
Di = D - thk
He = 0.75*D
L=75*D #m
A = 3.14*(D**2-Di**2)/4
I = 3.14*(D**4 - Di**4)/64 #m^4
r = (I/A)**0.5
Def = (D/r)**2;
E = 200e9 # N/m^2
EI = E*I
ks = 100000 # N/m
qe = 2200*A*9.81
se = E*A*0.00001*40


H = He/r
S = se*L**2/EI;
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
plt.plot(x,-y, 'r',marker='.')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,-y1, 'b', marker='.', label= r'$EI = 98.7 \times 10^6 N m^2$')
plt.legend(loc='upper right', prop = {'size' : 14}) #plt.legend(loc='upper right')
plt.xlabel('x', fontsize=14)
plt.ylabel('w', fontsize=14)
 
 
 
E = 150e9 # N/m^2
EI = E*I
se = E*A*0.00001*40
 
 
H = He/r
S = se*L**2/EI;
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
plt.plot(x,-y, 'r',marker='p')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,-y1, 'b', marker='p', label=r'$EI = 74.1 \times 10^6 N m^2$')
plt.legend(loc='upper right', prop = {'size' : 14}) #plt.legend(loc='upper right')
plt.xlabel('x', fontsize=14)
plt.ylabel('w', fontsize=14)
 

E = 100e9 # N/m^2
EI = E*I
se = E*A*0.00001*40
 
 
H = He/r
S = se*L**2/EI;
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
plt.plot(x,-y, 'r',marker='x')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,-y1, 'b', marker='x', label=r'$EI = 49.4 \times 10^6 N m^2$')
plt.legend(loc='upper right', prop = {'size' : 14}) #plt.legend(loc='upper right')
plt.xlabel('x', fontsize=14)
plt.ylabel('w', fontsize=14)
 
 
E = 50e9 # N/m^2
EI = E*I
se = E*A*0.00001*40
 
 
H = He/r
S = se*L**2/EI;
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
plt.plot(x,-y, 'r',marker='4')

x1 = np.linspace(-1,0,50)
y1 = BeamOnFreeSpanResponse(x1)
  
plt.plot(x1,-y1, 'b', marker='4', label=r'$EI = 24.7 \times 10^6 N m^2$')
plt.legend(loc='upper right', prop = {'size' : 14}) #plt.legend(loc='upper right')
plt.xlabel('x', fontsize=14)
plt.ylabel('w', fontsize=14)
  
plt.show()

 
