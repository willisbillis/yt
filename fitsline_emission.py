from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from lmfit.models import GaussianModel

################################################################################
# load in fits data
hdulist = fits.open("gass_78_-48_1522043360.fits")
hdu = hdulist[0]
data = hdu.data
x,y = 110,30

hdu.header['CUNIT1'] = 'deg'
hdu.header['CUNIT2'] = 'deg'
hdu.header['CUNIT3'] = "m/s"
w = WCS(hdu.header,naxis=2)
x_dim = hdu.header['NAXIS1'] # x (deg), 187 pixels
y_dim = hdu.header['NAXIS2'] # y (deg), 251 pixels
z_dim = 575 #hdu.header['NAXIS3'] # z (vel), 1201 frames
dz = hdu.header['CDELT3']
z_min = hdu.header['CRVAL3']

z_max = z_min + z_dim*dz
vel = np.array([i/1000 for i in np.arange(z_min,z_max,dz)])
K = []
for z in xrange(0,z_dim):
    k = data[z,y,x]
    K.append(k)
K = np.array(K)
mod = GaussianModel()
pars = mod.guess(K, x=vel)
out  = mod.fit(K, pars, x=vel)
result = np.array(out.params)
sigma = result[0]
print "Velocity Dispersion of peak: " + "%.2f" % sigma + " km/s"

plt.plot(vel,K,label="Emission line")
plt.plot(vel,out.best_fit,label="Gaussian fit")
plt.legend()
plt.xlabel("Velocity (km/s)")
plt.ylabel("K")
plt.savefig("fits_emission({0},{1}).png".format(x,y),dpi=750)
