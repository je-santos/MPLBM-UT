import urllib, hdf5storage, scipy.ndimage, vedo # import libs

# download image from the Digital Rocks Portal
drp = 'https://www.digitalrocksportal.org/'
url = drp + 'projects/374/images/309544/download/'
urllib.request.urlretrieve( url, im_name:='374_01_00_256.mat' )

# read-in image
bin_im = hdf5.storage.loadmat(im_name)['bin']

# obtain the Euclidean distance of the pore-space
edist  = scipy.ndimage.morphology.distance_transform_edt(bin_im==0)

# remove upper corner
half_size = bin_im.shape[0]//2
edist[ -half_size:, -half_size:, :half_size] =  0
bin_im[-half_size:, -half_size:, :half_size] = -1

# plot using Vedo
plot  = vedo.Volume(edist).legosurface(vmin=1, vmax=5,  cmap='turbo')
plot += vedo.Volume(bin_im).legosurface(vmin=1, vmax=2).c('lightgray').opacity(0.05)
vedo.plotter(plot)
