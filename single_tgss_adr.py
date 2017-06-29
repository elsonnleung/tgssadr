"""
SCRIPT USED TO RUN TGSS ADR (150MHz) SURVEY FITS FILES

@author: elsonleung
"""
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy import wcs
from matplotlib.patches import Ellipse
from average import average


NameJ = 'PSRJ0218+4232.FITS'

psr = fits.open(NameJ)

data = np.genfromtxt('PSR_RA_DEC.txt', dtype = 'str', skip_header = 1)
beamsize = np.loadtxt(NameJ[:-5]+'.txt', usecols = (10, 12)) #pulls out average beam sizes from catalogue

#average beam sizes
bmaj = average(beamsize)[0]
bmin = average(beamsize)[1]


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

#distances from the source
dist = np.sqrt((xpix - xv)**2. + (ypix - yv)**2.)

#find coordinates in which are within the radius of beam
xIn_Circ, yIn_Circ = np.where(dist <= beam_area)
Inner_CircInd = np.where(dist <= beam_area)

Inner_CircVals = image2D[Inner_CircInd]
mu_inner = sum(Inner_CircVals)/len(Inner_CircVals)

#check which pixels are plotted
plt.plot(yIn_Circ, xIn_Circ, 'r,')


###################################
#### OUTER INNER CIRCLE ###########
###################################

xOut_Circ , yOut_Circ = np.where(dist <= 50.0)
Outer_CircInd = np.where(dist <= 50.0)

Outer_CircVals = image2D[Outer_CircInd]
mu_outer = sum(Outer_CircVals)/len(Outer_CircVals)

#plt.plot(yOut_Circ, xOut_Circ, 'g,')


##################################
###### CALCULATIONS ##############
##################################
mu_inner = sum(Inner_CircVals)/len(Inner_CircVals) #average in the inner circle
mu_outer = sum(Outer_CircVals)/len(Outer_CircVals) #average in the outer circle

                  
#plot the circle
if float(DecJ[:3]) > 19.0:
    ellipse = Ellipse(xy=(xpix, ypix), width = beam_area, height = beam_area , edgecolor='y', fc = 'None', lw = 0.1)
    areaEllip = bmaj*bmin*np.pi #these two lines, something wrong
else:
    ellipse = Ellipse(xy=(xpix, ypix), width = abs(beam_area/np.cos(np.radians(Dec-19.0))), height = beam_area, edgecolor='y', fc = 'None', lw = 0.1)
    areaEllip = beam_area*abs(beam_area/np.cos(np.radians(Dec-19.0)))*np.pi

ax.add_patch(ellipse)
ax.add_artist(circle)

InnerCirc_beam = len(Inner_CircVals)/beam_area #number of beams inside the INNER circle
OuterCirc_beam = len(Outer_CircVals)/beam_area #number of beams inside the OUTER circle

InnerFlux = mu_inner * InnerCirc_beam #inner beam flux: average of inner * number of beams inner
OuterFlux = mu_outer * OuterCirc_beam #outer beam flux: average of outer * number of beams outer

FluxDiff = OuterFlux - InnerFlux #difference in flux
BeamDiff = OuterCirc_beam - InnerCirc_beam #difference in beams

AnnulusBg = FluxDiff/ BeamDiff #the background value in the annulus

bg = AnnulusBg * InnerCirc_beam #the background: Annulus background * number of beams inner

SourceFlux = mu_inner * InnerCirc_beam - bg #source flux: average of inner * number of beams - background

print SourceFlux
print max(Inner_CircVals)


original.savefig('%s.png' %NameJ[:-5], dpi = 1500)










