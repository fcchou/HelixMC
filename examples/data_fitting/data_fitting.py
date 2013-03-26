from helixmc import fitfxn
import numpy as np
from scipy import odr, stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Fitting log(F) vs. Z-extension with Bouchiat's WLC model #
z_vs_f = np.loadtxt('Z_vs_F.out')
Z = z_vs_f[:,1]
F = z_vs_f[:,0]
logF = np.log10(F)

# Fitting with inextensible WLC
def fit_func1(z, Lp, L0):
    'Return log(F) with input z-extension, bending persistense, and contour length.'
    return np.log10( fitfxn.wlc_bouchiat(Lp, L0, z) )

L0_start = 13000
Lp_start = 500

logF1 = logF[:13] # Fit to medium-low force regime only.
Z1 = Z[:13]

[Lp, L0], pcov = curve_fit( fit_func1, Z1, logF1, p0 = [Lp_start, L0_start] )
z_fit = np.linspace(Z1[0], Z1[-1], 100)
logF_fit = [ fit_func1(z, Lp, L0) for z in z_fit ]

plt.figure()
plt.title( 'Inextensible WLC')
plt.plot(Z1/10000, logF1, 'ro') # Convert unit to micrometer
plt.plot(z_fit/10000, logF_fit, 'b-')
plt.figtext(0.2, 0.8, 'Lp = %.1f nm' % (Lp/10))
plt.xlabel(r'end-to-end ($\mu m$)')
plt.ylabel('log(F) (pN)')

# Fitting with extensible WLC
def fit_func2(params, x):
    'Implicit fitting function for extensible WLC.'
    z = x[0] #params = [Lp, L0, S], x = [Z, logF]
    F = 10 ** (x[1]) #Convert log(F) to F
    zero = fitfxn.wlc_bouchiat_impl( *params, z=z, F=F )
    return zero

params_start = [300,18000,2000] #Lp, L0, S
lsc_data = odr.Data( np.vstack( (Z,logF) ), y = 1 )
lsc_model = odr.Model( fit_func2, implicit=True )
lsc_odr = odr.ODR( lsc_data, lsc_model, params_start )
lsc_out = lsc_odr.run()
Lp, L0, S = lsc_out.beta

z_fit = np.linspace(Z[0], Z[-1], 100)
logF_fit = np.log10( fitfxn.f_wlc_bouchiat_impl(Lp, L0, S, z_fit) )

plt.figure()
plt.title( 'Extensible WLC')
plt.plot(Z/10000, logF, 'ro') # Convert unit to micrometer
plt.plot(z_fit/10000, logF_fit, 'b-')
plt.figtext(0.2, 0.80, 'Lp = %.1f nm' % (Lp/10))
plt.figtext(0.2, 0.75, 'S  = %.1f pN' % S)

# Fitting Effective Torsional Persisten vs. Force with Moroz-Nelson
# model
ceff_vs_f = np.loadtxt('Ceff_vs_F.out')
Ceff = ceff_vs_f[:,1]
F = ceff_vs_f[:,0]

def fit_func3(F, C) :
    return fitfxn.moroz_3rd(Lp, C, F)

[C], pcov = curve_fit( fit_func3, F, Ceff, p0 = [650] )
x = np.linspace(0.02, F[-1], 1000)
y = np.array(map( lambda x:fit_func3(x, C), x))
for i,j in enumerate(y):
    if j > Ceff[0]:
        break
x = x[(i-1):]
y = y[(i-1):]

plt.figure()
plt.plot(F, Ceff/10, 'ro')
plt.plot(x, y/10, 'b-')
plt.title('Moroz_Nelson')
plt.figtext(0.2, 0.8, 'C = %.1f nm' % (C/10))
plt.xlabel('Force (pN)')
plt.ylabel('Effective torsional persitence (nm)')


plt.show()
