import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize as opt




def foo(x,y,p):
	A = ((p**2 + y**2)**0.5 - p)**0.5
	B = ((p**2 + y**2)**0.5 + p)**0.5
	Zet = ((x**4-y**2)/4)**0.25
	a = np.array([[1, 0, 1, 0, -1, 0], [0, A, 0, B, Zet, -Zet], [A**2, 0, -B**2, 0, 0, 2*Zet**2], [0, A**3, 0, -B**3, -2*Zet**3, -2*Zet**3], [np.cosh(A), -np.sinh(A), np.cos(B), -np.sin(B), 0, 0], [A*np.sinh(A), -A*np.cosh(A), -B*np.sin(B), -B*np.cos(B), 0, 0]])
	return np.linalg.det(a);

X= np.linspace(2*np.pi, 50*np.pi, 500)

D = 0.452      #m
thk = D/15. #DNVGL-ST-F101 (2017),wall-thickness-to-outer-diameter ratio (D t) of 15.
Di = D - thk
A = 3.14*(D**2-Di**2)/4
I = 3.14*(D**4 - Di**4)/64 #m^4
r = (I/A)**0.5
Def = (D/r)**2;

pp=-Def
p0=0
Z = 0*X;
Z0 = 0*X;

for i in range(len(X)):
	Z[i] = opt.brentq(lambda y: foo(X[i], y,pp), 1, 30)
	Z0[i] = opt.brentq(lambda y: foo(X[i], y,p0), 1, 30)
	#print(X[i],Z[i],Z0[i])


plt.plot(X, np.pi/Z**0.5, 'r--', X, np.pi/Z0**0.5, 'b-.')
plt.xlabel(r'$\kappa=\left(\frac{K_wL^4}{EI} \right)^{\frac{1}{4}}$', fontsize=14)
plt.ylabel(r'$\frac{L_e}{L}}$', fontsize=14)
plt.legend(['midplane stretching', 'Small strain'],prop = {'size' : 14})
plt.ylim((0.6,0.9))

plt.tight_layout()
plt.show()
