from matplotlib import pyplot as plt
import matplotlib
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import astropy.units as u
import astropy.coordinates as ac
import cartopy.crs as ccrs
from astropy.coordinates import ICRS
import numpy as np

hdulist = fits.open("gass_78_-48_1522043360.fits")
hdu = hdulist[0]
data = hdu.data[538]

hdu.header['CUNIT1'] = 'deg'
hdu.header['CUNIT2'] = 'deg'
hdu.header['CUNIT3'] = "m/s"

w = WCS(hdu.header,naxis=2)
"""
w.printwcs()
NAXIS1 = hdu.header['NAXIS1']
NAXIS2 = hdu.header['NAXIS2']
x = np.arange(NAXIS1) * u.degree
y = np.arange(NAXIS2) * u.degree
X, Y = np.meshgrid(x, y)
ra_data, dec_data = w.wcs_pix2world(X, Y, 0)

w.wcs.crpix = [94.0, 719.25]
#w.wcs.cdelt = np.array([-0.1184149301180, 0.08])
w.wcs.crval = [77.5, 0]
w.wcs.ctype = ["GLON-CAR", "GLAT-CAR"]

origin_frame = ac.ICRS
rep = ac.SkyCoord(ra_data,dec_data,0, unit = 'deg', frame = 'fk5')
"""
fig = plt.figure()
ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8],wcs=w)

x,y = 110,10
pixx,pixy = w.wcs_world2pix(x,y,0)
print pixx,pixy

lon,lat = w.wcs_pix2world(pixx,pixy,0)
print lon, lat
fig.add_axes(ax)
plt.plot(lon,lat,'rx',markersize=5)
plt.plot(110,30,'bx',markersize=5)
ax.imshow(data, origin='lower')
ax.grid(False)
plt.pcolor(data, vmin=min(data.flatten()), vmax=max(data.flatten()), cmap=plt.cm.viridis)

icrs = ICRS()
overlay = ax.get_coords_overlay(icrs)
overlay.grid()

plt.xlabel('longitude')
plt.ylabel('latitude')
plt.savefig("GASS.png",dpi=750,bbox_inches='tight')
