"""
SCRIPT USED TO RUN TGSS ADR (150MHz) SURVEY FITS FILES

@author: elsonleung
"""
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
from matplotlib.patches import Ellipse


NameJ = 'PSRJ2021+3651.FITS'

psr = fits.open(NameJ)

data = np.genfromtxt('PSR_RA_DEC.txt', dtype = 'str', skip_header = 1)

#splitting the data array back into lists
for k in data:
    if k[0] == NameJ[3:-5]:
        RaJ = k[1]
        DecJ = k[2]
        Edot = k[3]

Ra = (float(RaJ[:2]) + float(RaJ[3:5])/60 + float(RaJ[6:])/3600)*15

print 'The Ra is', Ra

if DecJ[0] == '-':
    Dec = (-float(DecJ[1:3]) - float(DecJ[4:6])/60 - float(DecJ[7:])/3600)
else:
    Dec = (float(DecJ[1:3]) + float(DecJ[4:6])/60 + float(DecJ[7:])/3600)
print 'The Dec is', Dec

#import header
header = psr[0].header
image_data = psr[0].data
mywcs = wcs.WCS(header)

#convert the 4D array to 2D to make things simpler
image2D = image_data[0][0]
y_size, x_size = image2D.shape

#ra, dec, freq, stokes = [psr[0].header['crval1'], psr[0].header['crval2'], 0, 1]
ra, dec, freq, stokes = [Ra, Dec, psr[0].header['crval2'], 1]
(xpix,ypix, freq, stokes), = mywcs.wcs_world2pix([[ra,dec, freq, stokes]], 0)

#grab the xpixel and ypixel coordinate
print 'The corresponding xcoordinate and ycoordinate is', (xpix, ypix)

deg = 1.0
beamsize = 25.0
pix_per_degree = len(image_data[0][0][0])/deg #######there might be something wrong here######
pix_per_arcSec = pix_per_degree/3600.0
pix_per_beam = pix_per_arcSec * beamsize
print 'Amount of pixels per beam is', pix_per_beam

#create an array of all values
vals = image_data.flatten()

#calculate rms for sigma to grey scale    
rms = np.sqrt(sum(vals)**2/len(vals))
#print rms

                 
#shows original figure
original = plt.figure()
#a = original.add_subplot(2,1,1)
plt.imshow(image2D, cmap = 'gray', vmin = -1.*rms, vmax = 3*rms)
plt.colorbar()
original.show()
plt.title('%s' %NameJ[:-5])
#plt.plot([xpix],[ypix], 'yo', markerfacecolor = 'none',  markersize = 20, markeredgewidth = 0.1)
ax=plt.gca()
ax.invert_yaxis()
circle = plt.Circle((50.0, int(y_size)-50.0), 50.0, edgecolor = 'g', fc = 'None', lw = 0.5)

###################################
#### FOR INNER CIRCLE #############
###################################
#create a gride with values on them
x = np.arange(0, x_size)
y = np.arange(0, y_size)

xv, yv = np.meshgrid(x,y)

#print xv
#print yv

#distances from the source
dist = np.sqrt((xpix - xv)**2. + (ypix - yv)**2.)

#find coordinates in which are within the radius of beam
xIn_Circ, yIn_Circ = np.where(dist <= pix_per_beam*1.5)
Inner_CircInd = np.where(dist <= pix_per_beam*1.5)

Inner_CircVals = image2D[Inner_CircInd]
mu_inner = sum(Inner_CircVals)
#Source_flux = sum(Inner_CircVals)/pix_per_beam ######or is it the area/ number of pixels???#######

#check which pixels are plotted
plt.plot(yIn_Circ, xIn_Circ, 'r,')

"""
###################################
#### OUTER INNER CIRCLE ###########
###################################

xOut_Circ , yOut_Circ = np.where(dist <= 50.0)
Outer_CircInd = np.where(dist <= 50.0)

Outer_CircVals = image2D[Outer_CircInd]

#plt.plot(yOut_Circ, xOut_Circ, 'r,')

#find the background flux
BG_flux = (sum(Outer_CircVals)-sum(Inner_CircVals))/len(Outer_CircVals) #instead of area, divided by number of pixels
print 'The background flux is', BG_flux

if max(vals.flatten()) <= Source_flux*1.1 and max(vals.flatten()) >= Source_flux*0.9: #within 10% of the peak flux or integrated flux??
    print 'The source is resolved. The source flux is', Source_flux, 'and the peak value is', max(vals.flatten())
else:
    print 'The source is unresolved. The source flux is', Source_flux, 'but the peak value is', max(vals.flatten())
"""
###################################
#### OUTER INNER CIRCLE ###########
###################################
#distance to radius of outer circle
dist_out = np.sqrt((50. - xv)**2 + (int(y_size) - 50.0 - yv)**2)

#find coordinates in which are within the radius of the outer circle
xOut_Circ, yOut_Circ = np.where(dist_out <= 50.0)
Outer_CircInd = np.where(dist_out <= 50.0)
Outer_CircVals = image2D[Outer_CircInd]

#plotted just to check
plt.plot(yOut_Circ, xOut_Circ, 'r,')


"""
#find the background flux
BG_flux = (sum(Outer_CircVals)-sum(Inner_CircVals))/len(Outer_CircVals) #instead of area, divided by number of pixels
print 'The background flux is', BG_flux

#calibrated flux
caliFlux = sum(Inner_CircVals) - len(Inner_CircVals)*BG_flux #background flux is negative so I'm effectively ADDING on the flux???
#print caliFlux

#check if image is resolved or unresolved
if max(Inner_CircVals) <= caliFlux*1.1 and max(Inner_CircVals) >= caliFlux*0.9: #within 10% of the peak flux or integrated flux??
    print 'The source is resolved. The calibrated flux is', caliFlux, 'and the peak value is', max(vals.flatten())
else:
    print 'The source is unresolved. The calibrated flux is', caliFlux, 'but the peak value is', max(vals.flatten())
"""
  
#plot the circle
if float(DecJ[:3]) > 19.0:
    bmaj = 25.0
    bmin = 25.0
    ellipse = Ellipse(xy=(xpix, ypix), width = pix_per_beam, height = pix_per_beam , edgecolor='y', fc = 'None', lw = 0.1)
    areaEllip = pix_per_beam*pix_per_beam*np.pi
else:
    bmaj = abs(pix_per_beam/np.cos(np.radians(Dec-19.0)))
    bmin = 25.0
    ellipse = Ellipse(xy=(xpix, ypix), width = bmaj, height = bmin, edgecolor='y', fc = 'None', lw = 0.1)
    areaEllip = pix_per_beam*abs(pix_per_beam/np.cos(np.radians(Dec-19.0)))*np.pi

ax.add_patch(ellipse)
ax.add_artist(circle)

##################################
###### CALCULATIONS ##############
##################################
cdelt1 = psr[0].header['CDELT1']*3600
cdelt2 = psr[0].header['CDELT2']*3600
beam_area = 1.1331 * ((bmaj*bmin)/abs(cdelt1*cdelt2))

InnerCirc_beam = len(Inner_CircVals)/beam_area
OuterCirc_beam = len(Outer_CircVals)/beam_area

InnerFlux = mu_inner * InnerCirc_beam




#original.savefig('%s.png' %NameJ[:-5], dpi = 1500)










