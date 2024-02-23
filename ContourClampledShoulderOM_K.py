import matplotlib.pyplot as plt
import numpy as np

#plt.style.use('_mpl-gallery-nogrid')

# make data


X, Y = np.meshgrid(np.linspace(0, 50*np.pi, 500), np.linspace(0, 68*np.pi, 500))
Z = 0*X;

D = 0.452      #m
thk = D/15. #DNVGL-ST-F101 (2017),wall-thickness-to-outer-diameter ratio (D t) of 15.
Di = D - thk
A = 3.14*(D**2-Di**2)/4
I = 3.14*(D**4 - Di**4)/64 #m^4
r = (I/A)**0.5
Def = (D/r)**2;

p=-Def
for i in range(len(Y)):
	for j in range(len(Y)):
		A = ((p**2 + Y[i,j]**2)**0.5 - p)**0.5
		B = ((p**2 + Y[i,j]**2)**0.5 + p)**0.5
		Zet = ((X[i,j]**4-Y[i,j]**2)/4)**0.25
		a = np.array([[1, 0, 1, 0, -1, 0], [0, A, 0, B, Zet, -Zet], [A**2, 0, -B**2, 0, 0, 2*Zet**2], [0, A**3, 0, -B**3, -2*Zet**3, -2*Zet**3], [np.cosh(A), -np.sinh(A), np.cos(B), -np.sin(B), 0, 0], [A*np.sinh(A), -A*np.cosh(A), -B*np.sin(B), -B*np.cos(B), 0, 0]])
		Z[i,j] =  np.linalg.det(a);

Z0 = 0*X;
p=0
for i in range(len(Y)):
	for j in range(len(Y)):
		A = ((p**2 + Y[i,j]**2)**0.5 - p)**0.5
		B = ((p**2 + Y[i,j]**2)**0.5 + p)**0.5
		Zet = ((X[i,j]**4-Y[i,j]**2)/4)**0.25
		a = np.array([[1, 0, 1, 0, -1, 0], [0, A, 0, B, Zet, -Zet], [A**2, 0, -B**2, 0, 0, 2*Zet**2], [0, A**3, 0, -B**3, -2*Zet**3, -2*Zet**3], [np.cosh(A), -np.sinh(A), np.cos(B), -np.sin(B), 0, 0], [A*np.sinh(A), -A*np.cosh(A), -B*np.sin(B), -B*np.cos(B), 0, 0]])
		Z0[i,j] =  np.linalg.det(a);
# plot
fig, ax = plt.subplots()

#ax.contour(X, Y, Z, levels=0)

cntr1 = ax.contour(X, Y, Z, levels=0)
cntr2 = ax.contour(X, Y, Z0, levels=0, linestyles="--", colors='C0')
h1,_ = cntr1.legend_elements()
h2,_ = cntr2.legend_elements()
ax.legend([h1[0], h2[0]], ['midplane stretching', 'Small strain'],prop = {'size' : 14})
ax.set_xlabel(r'$\kappa=\left(\frac{K_w L^4}{EI} \right)^{\frac{1}{4}}$', fontsize=14)
ax.set_ylabel(r'$\Omega = \omega\sqrt{ \frac{m_e L^4}{EI}}$', fontsize=14)

plt.tight_layout()
plt.show()
