from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization.wcsaxes import WCSAxes
import numpy as np
from lmfit.models import GaussianModel
from multiprocessing import Pool
import os,h5py,time

################################################################################
# load in fits data
hdulist = fits.open("gass_78_-48_1522043360.fits")
hdu = hdulist[0]
data = hdu.data

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

no_signal = []
case2 = []
case7 = []
case8 = []
def line_emission(y):
    z_dim = 580
    dz = hdu.header['CDELT3']
    z_min = hdu.header['CRVAL3']
    z_max = z_min + z_dim*dz
    vel = np.array([i/1000 for i in np.arange(z_min,z_max,dz)])
    x_dims = [i for i in xrange(0,x_dim)]
    vbin = [0 for i in x_dims]

    while len(x_dims) > 0:
        print y, len(x_dims)
        for x in x_dims:
            if (y-10) > 142.0/186.0*x:
                x_dims.remove(x)
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
                        vbin[x] = center
                        x_dims.remove(x)
                    else:
                        case2.append((x,y))
                        x_dims.remove(x)
                elif height > 0.07 and center > -35.0:
                    continue
                elif height < 0.07 and center < -35.0:
                    no_signal.append((x,y))
                    x_dims.remove(x)
                elif height < 0.07 and center > -35.0:
                    if sigma < 30.0:
                        case7.append((x,y))
                        x_dims.remove(x)
                    else:
                        case8.append((x,y))
                        x_dims.remove(x)
        z_dim = z_dim - 10
        vel = vel[:-10]

    extended_vbin = [y] + vbin
    return extended_vbin

if __name__ == '__main__':
    pool = Pool(processes=3)
    y_dim = 249
    chunksize = 3 # best if a multiple of processes for efficiency
    if y_dim % chunksize == 0:
        for a in np.arange(0,y_dim,chunksize):
            inputs = np.arange(a,a+chunksize)
            datachunk = np.array(pool.map(line_emission, inputs))
            if os.path.isfile('center.h5'):
                with h5py.File('center.h5', 'a') as hf:
                    hf['alldata'].resize((hf['alldata'].shape[0] + datachunk.shape[0]), axis=0)
                    hf['alldata'][-datachunk.shape[0]:] = datachunk
            else:
                with h5py.File("center.h5",'w') as hf:
                    hf.create_dataset('alldata',data=datachunk,maxshape=(None,188))
    else:
        print "choose different chunksize. remainder = {0}".format(y_dim % chunksize)

    # use data from hdf5 file to make plot
    with h5py.File('center.h5','r') as hf:
        data = hf['alldata'][:]
        data = data[data[:,0].argsort()] #reorder data by first entry of each element
        data = np.array([list(i[1:]) for i in data]) # take off ordering and
                                                     # returns np.array
        fig = plt.figure()
        ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8],wcs=w)
        fig.add_axes(ax)
        im = plt.imshow(data,origin = "lower",cmap='jet')

        plt.xlabel("longitude")
        plt.ylabel("latitude")
        cbar = plt.colorbar(pad = 0)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label("Velocity Dispersion (km/s)",size=7)
        plt.savefig("fits_centers.png",dpi=750,bbox_inches='tight')
