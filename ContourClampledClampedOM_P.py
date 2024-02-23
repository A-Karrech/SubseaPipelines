import matplotlib.pyplot as plt
import numpy as np

#plt.style.use('_mpl-gallery-nogrid')

# make data



def f(x,y):
	A = ((x**2 + y**2)**0.5 - x)**0.5
	B = ((x**2 + y**2)**0.5 + x)**0.5
	return 2*A*B*(1-np.cosh(A)*np.cos(B))+(A**2-B**2)*np.sinh(A)*np.sin(B)
   
 

X, Y = np.meshgrid(np.linspace(-50*np.pi, 50*np.pi, 256), np.linspace(0, 100*np.pi, 256))

D = 0.452      #m
thk = D/15. #DNVGL-ST-F101 (2017),wall-thickness-to-outer-diameter ratio (D t) of 15.
Di = D - thk
A = 3.14*(D**2-Di**2)/4
I = 3.14*(D**4 - Di**4)/64 #m^4
r = (I/A)**0.5
Def = (D/r)**2;

p=-Def
Z = f(X+p,Y)
Z0 = f(X,Y)
#levels = np.linspace(np.min(Z), np.max(Z), 1)

# plot
fig, ax = plt.subplots()

cntr1 = ax.contour(X, Y, Z, levels=0)
cntr2 = ax.contour(X, Y, Z0, levels=0, linestyles="--", colors='C0')
h1,_ = cntr1.legend_elements()
h2,_ = cntr2.legend_elements()
ax.legend([h1[0], h2[0]], ['midplane stretching', 'Small strain'],prop = {'size' : 14})
ax.set_xlabel(r'$P= \frac{S_eL^2}{2EI}-\Lambda$', fontsize=14)
ax.set_ylabel(r'$\Omega = \omega\sqrt{ \frac{m_e L^4}{EI}}$', fontsize=14)

plt.tight_layout()
plt.show()
