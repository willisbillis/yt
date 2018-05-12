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
x,y = 0,11

hdu.header['CUNIT1'] = 'deg'
hdu.header['CUNIT2'] = 'deg'
hdu.header['CUNIT3'] = "m/s"
w = WCS(hdu.header,naxis=2)
x_dim = hdu.header['NAXIS1'] # x (deg), 187 pixels
y_dim = hdu.header['NAXIS2'] # y (deg), 251 pixels
z_dim = 580 #hdu.header['NAXIS3'] # z (vel), 1201 frames
dz = hdu.header['CDELT3']
z_min = hdu.header['CRVAL3']

z_max = z_min + z_dim*dz
vel = np.array([i/1000 for i in np.arange(z_min,z_max,dz)])
sigma = 0
while True:
    if (y-10) > 142.0/186.0*x:
        sigma = 0
        return False
    else:
        K = []
        for z in xrange(0,z_dim):
            k = data[z,y,x]
            if np.isnan(k) == True:
                K.append(0)
            else:
                K.append(k)
        K = np.array(K)
        mod = GaussianModel()
        pars = mod.guess(K, x=vel)
        out  = mod.fit(K, pars, x=vel)
        result = np.array(out.params)
        sigma = result[0]
        center = result[1]
        height = result[4]

        if height > 0.07 and center < -35.0:
            if sigma < 30.0:
                vbin[x] = sigma
                return False
            else:
                case2.append((x,y))
                return False
        elif height > 0.07 and center > -35.0:
            continue
        elif height < 0.07 and center < -35.0:
            no_signal.append((x,y))
            return
        elif height < 0.07 and center > -35.0:
            if sigma < 30.0:
                case7.append((x,y))
                return False
            else:
                case8.append((x,y))
                return False
    z_dim = z_dim - 10
    vel = vel[:-10]
print "Velocity Dispersion of peak: " + "%.2f" % sigma + " km/s"

plt.plot(vel,K,label="Emission line")
plt.plot(vel,out.best_fit,label="Gaussian fit")
plt.legend()
plt.xlabel("Velocity (km/s)")
plt.ylabel("K")
plt.savefig("fits_emission({0},{1}).png".format(x,y),dpi=750)
