"""
SCRIPT USED TO RUN TGSS ADR (150MHz) SURVEY FITS FILES

@author: elsonleung
"""
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
from matplotlib.patches import Ellipse


NameJ = 'PSRJ0218+4232.FITS'

psr = fits.open(NameJ)

data = np.genfromtxt('PSR_RA_DEC.txt', dtype = 'str', skip_header = 1)
beamsize = np.loadtxt(NameJ[:-5]+'.txt', usecols = (0, 10, 12)) #pulls out average beam sizes from catalogue

print 'The beamsized used is', beamsize[:,0][0], 'arcminutes away'

##average beam sizes
bmaj = np.median(beamsize[:,1])
bmin = np.median(beamsize[:,2])

#print bmaj,bmin

#splitting the data array back into lists
for k in data:
    if k[0] == NameJ[3:-5]:
        RaJ = k[1]
        DecJ = k[2]
        Edot = k[3]


Ra = (float(RaJ[:2]) + float(RaJ[3:5])/60 + float(RaJ[6:])/3600)*15

#print 'The Ra is', Ra

if DecJ[0] == '-':
    Dec = (-float(DecJ[1:3]) - float(DecJ[4:6])/60 - float(DecJ[7:])/3600)
else:
    Dec = (float(DecJ[1:3]) + float(DecJ[4:6])/60 + float(DecJ[7:])/3600)
#print 'The Dec is', Dec

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
circle = plt.Circle((xpix, ypix), 50.0, edgecolor = 'g', fc = 'None', lw = 0.5)



#########################################################################################################################
    
cdelt1 = psr[0].header['CDELT1']*3600 #pull this value from the header
cdelt2 = psr[0].header['CDELT2']*3600 #pull this value from the header
beam_area = 1.1331 * ((bmaj*bmin)/abs(cdelt1*cdelt2)) #the equation given to me    

########################## PLOT THE BEAM SIZE ###########################################################################

ellipse = Ellipse(xy=(xpix, ypix), width = bmaj/6.2*5, height = bmin/6.2*5, edgecolor='y', fc = 'None', lw = 0.5)

ax.add_patch(ellipse)
ax.add_artist(circle)
ax.legend([circle, ellipse], ['Region 2','Region 1'])


#########################################################################################################################
    
    
   
###################################
#### FOR INNER CIRCLE #############
###################################
#create a gride with values on them
x = np.arange(0, x_size)
y = np.arange(0, y_size)

xv, yv = np.meshgrid(x,y)

#print xv
#print yv

dx = xv - xpix
dy = yv - ypix

#distance to source
dist = (dx**2/(bmaj/6.2/2*5)**2) + (dy**2/ (bmin/6.2/2*5)**2) #remember to change it from arcsec to pixels

#find coordinates in which are within the radius of beam
xIn_Circ, yIn_Circ = np.where(dist <= 1.)
Inner_CircInd = np.where(dist <= 1.)

Inner_CircVals = image2D[Inner_CircInd]
npixInner = len(Inner_CircVals)
mu_inner = np.mean(Inner_CircVals)

#check which pixels are plotted
#plt.plot(yIn_Circ, xIn_Circ, 'r,')


###################################
#### OUTER INNER CIRCLE ###########
###################################

xOut_Circ , yOut_Circ = np.where(dist <= 50.0)
Outer_CircInd = np.where(dist <= 50.0)

Outer_CircVals = image2D[Outer_CircInd]
npixOuter = len(Outer_CircVals) - npixInner
mu_outer = (sum(Outer_CircVals)-sum(Inner_CircVals))/npixOuter
#plt.plot(yOut_Circ, xOut_Circ, 'g,')

###################################
#########CALCULATIONS##############
###################################

innersum = sum(Inner_CircVals)
innermean = mu_inner - abs(mu_outer)
innermax = max(Inner_CircVals) - abs(mu_outer)

NbeamsInner = npixInner/beam_area
total_flux = innermean * NbeamsInner


source_beam_area = 1.1331 * ((beamsize[:,1][0]*beamsize[:,2][0])/abs(cdelt1*cdelt2))
print bmaj, bmin
print beamsize[:,1][0], beamsize[:,2][0]
print cdelt1
print 'The beam area is', beam_area
print 'The source beam area is', source_beam_area


if beam_area < source_beam_area:
    if total_flux <= max(Inner_CircVals)*0.9: #condition so that its within 10%
        print "The calculation is wrong, you're a piece of shit and you know it"
    else:
        print "The source is resolved and the integrated flux is therefore", total_flux
else:
    print 'The source is unresolved and the flux is therefore', max(Inner_CircVals)


original.savefig('%s.png' %NameJ[:-5], dpi = 1500)






