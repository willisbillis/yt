from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import numpy as np
from lmfit.models import GaussianModel
import time


################################################################################
# load in fits data

hdulist = fits.open("gass_78_-48_1522043360.fits")
hdu = hdulist[0]
data = hdu.data
x,y = 110,10
hdu.header['CUNIT1'] = 'deg'
hdu.header['CUNIT2'] = 'deg'
hdu.header['CUNIT3'] = "m/s"
w = WCS(hdu.header,naxis=2)
"""
x_dim = hdu.header['NAXIS1'] # x (deg), 187 pixels
y_dim = hdu.header['NAXIS2'] # y (deg), 251 pixels
z_dim = 500#hdu.header['NAXIS3'] # z (vel), 1201 frames
z_min = hdu.header['CRVAL3']
dz = hdu.header['CDELT3']
z_max = z_min + z_dim*dz
vel = np.array([i/1000 for i in np.arange(z_min,z_max,dz)])

time_0 = time.time()
disp_array = []
for y in xrange(0,y_dim):
    print y
    vbin = []
    for x in xrange(0,x_dim):
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
        sigma = result[0]/1000
        vbin.append(sigma)
    disp_array.append(vbin)

da = np.array(disp_array)
np.save("fits_da",da)
time_1 = time.time()
print "time elapsed: " + "%.2f" % ((time_1-time_0)/60) + " min"
"""

da = np.load("fits_da.npy")

fig = plt.figure()
ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8],wcs=w)
fig.add_axes(ax)
im = plt.imshow(da,origin = "lower",cmap='jet')
#plt.plot(110,10,'ro',markersize=1)

plt.xlabel("longitude")
plt.ylabel("latitude")
cbar = plt.colorbar(pad = 0)
cbar.ax.tick_params(labelsize=8)
cbar.set_label("Velocity Dispersion (km/s)",size=7)
plt.savefig("fits_vel_disp.png",dpi=750,bbox_inches='tight')
