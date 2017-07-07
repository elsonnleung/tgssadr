#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 14:42:31 2017

@author: elsonleung
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
#from astropy.utils.data import download_file #only used for testing
from matplotlib.colors import LogNorm #scale log like in ds9
import glob, os #this is for looping the pictures with same .FITS
from astropy import wcs
from matplotlib.patches import Ellipse

count = 0
resolved = 0
unresolved = 0
shitcounter = 0
resolvedList = []
unresolvedList= []
psrList = []
os.chdir('/Users/elsonleung/Google Drive/PWN research/PHY479 - PWN research/TGSS ADR')
for files in glob.glob('*.FITS'):
    psrList.append(files)
    
data = np.genfromtxt('PSR_RA_DEC.txt', dtype = 'str', skip_header = 1)
nameJ = data[:,0]

visiblepsr = np.genfromtxt('visiblepsr.txt', dtype = 'str')
vispsrList = []

for i in visiblepsr:
    vispsrList.append(i)

for psr in vispsrList:
    for k in data:
        if k[0] == psr[3:-5]: #do this to get rid of PSR and .txt
            RaJ = k[1] #actual Ra
            DecJ = k[2] #actual Dec
            Edot = k[3] #Edot
    print '---------------', psr[:-5], '----------------------'      

            
    Ra = (float(RaJ[:2]) + float(RaJ[3:5])/60 + float(RaJ[6:])/3600)*15 #convert RaJ to degrees
#    print Ra
    if DecJ[0] == '-': #convert DecJ to degrees
        Dec = (-float(DecJ[1:3]) - float(DecJ[4:6])/60 - float(DecJ[7:])/3600) #if its negative, need to subtract all
    else:
        Dec = (float(DecJ[1:3]) + float(DecJ[4:6])/60 + float(DecJ[7:])/3600) #if positive, add all
#    print Dec
    
    beamsize = np.loadtxt(psr[:-5]+'.txt', usecols = (0,10,12)) #extract the beam sizes
    bmaj = np.median(beamsize[:,1]) #take the median of beamsizes
    bmin = np.median(beamsize[:,2])

    pulsar = fits.open(psr)
    header = pulsar[0].header
    #grab image data
    image_data = pulsar[0].data
                         
    #convert the 4D array to 2D to make things simpler
    image2D = image_data[0][0]
    y_size, x_size = image2D.shape
    
    mywcs = wcs.WCS(header)
    #convert ra and dec to pixels on the image
    ra, dec, freq, stokes = [Ra, Dec, pulsar[0].header['crval2'], 1]
    (xpix,ypix, freq, stokes), = mywcs.wcs_world2pix([[ra,dec, freq, stokes]], 0)


#create an array of all values
    vals = image_data.flatten()

#calculate rms for sigma to grey scale    
    rms = np.sqrt(sum(vals)**2/len(vals))
#    print rms
                         
    #shows original figure
    original = plt.figure()
    #a = original.add_subplot(2,1,1)
    plt.imshow(image2D, cmap = 'gray', vmin = -0.1*rms, vmax = 0.5*rms)
    plt.colorbar()
    original.show()
    plt.title('%s' %psr[:-5])
    #plt.plot([xpix],[ypix], 'yo', markerfacecolor = 'none',  markersize = 20, markeredgewidth = 0.1)
    ax=plt.gca()
    ax.invert_yaxis()
    circle = plt.Circle((xpix,ypix), 50.0, edgecolor = 'g', fc = 'None', lw = 0.5)
    
    #########################################################################################################################
    
    cdelt1 = pulsar[0].header['CDELT1']*3600 #pull this value from the header
    cdelt2 = pulsar[0].header['CDELT2']*3600 #pull this value from the header
    beam_area = 1.1331 * ((bmaj*bmin)/abs(cdelt1*cdelt2)) #the equation given to me    
    
    ########################## PLOT THE BEAM SIZE ###########################################################################
    
    ellipse = Ellipse(xy=(xpix, ypix), width = bmaj/6.2*5, height = bmin/6.2*5, edgecolor='y', fc = 'None', lw = 0.5)
    
    ax.add_patch(ellipse)
    ax.add_artist(circle)
    ax.legend([circle, ellipse], ['Region 2','Region 1'])

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
    
    ###################################
    #########CALCULATIONS##############
    ###################################
    
    innersum = sum(Inner_CircVals)
    innermean = mu_inner - abs(mu_outer)
    innermax = max(Inner_CircVals) - abs(mu_outer)
    
    NbeamsInner = npixInner/beam_area
    total_flux = innermean * NbeamsInner
    
    
    source_beam_area = 1.1331 * ((beamsize[:,1][0]*beamsize[:,2][0])/abs(cdelt1*cdelt2))
#    print bmaj, bmin
#    print beamsize[:,0][0], beamsize[:,0][1]
#    print cdelt1
    print 'The beam area is', beam_area
    print 'The source beam area is', source_beam_area
    
    
    if beam_area < source_beam_area:
        if total_flux <= max(Inner_CircVals)*0.9: #condition to check if its within 10%
            print psr[:-5],"'s calculation is wrong, you're a piece of shit and you know it"
            shitcounter += 1
        else:
            print psr[:-5]," is resolved and the integrated flux is therefore", total_flux,'mJy'
            resolvedList.append(psr[:-5])
    else:
        print psr[:-5],' is unresolved and the flux is therefore', max(Inner_CircVals),'mJy'
        unresolvedList.append(psr[:-5])
        
    original.savefig('%s.png' %psr[:-5], dpi = 1500)
     
    count += 1 #just to keep track which file I'm on
#    print 'We are on image number', count
#    print 'The current pulsar is', psr[:-5]

print 'shitcounter:', shitcounter




    





