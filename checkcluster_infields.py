# SCRIPT TO CHECK WHETHER A CLUSTER WITH GIVEN (RA,DEC) POSITION IS OBSERVED BY LoTSS
# It adds the field ID, separation to the field center, the noise from the rms 6'' map of each field, the mean value, and the expected noise from the Declination to the catalog
#
# G. Di Gennaro
# Jan 2025

import warnings
warnings.filterwarnings('ignore')

import pickle
import pandas as pd
import numpy as np
import sys, os, glob
#import pyfits
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table


def sepn(r1,d1,r2,d2):
  """
  Calculate the separation between 2 sources, RA and Dec must be
  given in radians. Returns the separation in radians [TWS]
  """
  cos_sepn = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
  sepn = np.arccos(cos_sepn)
  return sepn   

def gaussian_response(rasep, decsep, FWHM=3.96):
  """
  Calculate the response of the LOFAR antennas (approximating them to Gaussian)
  given the distance of the target from the pointing centre
  """ 
  s = 1.0 * np.exp(-( (rasep**2.0)/(2.0*(FWHM/2.355)**2.0) + (decsep**2.0)/(2.0*(FWHM/2.355)**2.0) ))
  return s


def analytical_noise(dec, phi=52.90888889, A=62E-6):
  #elev = np.pi/2. + (53.* (np.pi/180.)) - (pdec * (np.pi/180.))
  noise = A * np.cos(dec  - (phi * (np.pi/180.)))**(-2)
  return noise


def map_noise(ID, ra, dec, width=1.0):
  """
  Map noise provided by T. Shimwell
  We calculate the expected noise associated to the cluster area from the cluster coordinates and averaged within 0.5deg
  """

  highresrms, hdr = fits.getdata('/local/work/g.digennaro/LoTSS-DR3/high-rmsmaps/'+str(ID)+'-mosaic.rms.fits', header=True)

  w = WCS(hdr, naxis=2)
  xmin, ymin = w.wcs_world2pix(ra+(width/2)/np.cos(dec*np.pi/180.), dec-width/2, 1)
  xmax, ymax = w.wcs_world2pix(ra-(width/2)/np.cos(dec*np.pi/180.), dec+width/2, 1)
  xmin, ymin = int(xmin), int(ymin)
  xmax, ymax = int(xmax), int(ymax)

  if xmin < 0 : xmin = 0
  if xmax < 0 : xmax = 0
  if ymin < 0 : ymin = 0
  if ymax < 0 : ymax = 0

  rms = np.sqrt(np.nanmean(highresrms[0,0,ymin:ymax,xmin:xmax]**2))
  
  return (rms)

def convertDR3status(status):
  if status == 1: plan = 'Progress'
  if status == 2: plan = 'Final'
  if status == 'None': plan = 'None'

  return (plan)


def match_LoTSSpointingID_to_catalogue(clustercatalogue, dcoord=2.2):
  arcsec2deg=1.0/3600
  arcmin2deg=1.0/60
  deg2rad=np.pi/180
  deg2arcsec = 1.0/arcsec2deg
  rad2deg=180.0/np.pi
  arcmin2rad=arcmin2deg*deg2rad
  arcsec2rad=arcsec2deg*deg2rad
  rad2arcmin=1.0/arcmin2rad
  rad2arcsec=1.0/arcsec2rad
  steradians2degsquared = (180.0/np.pi)**2.0
  degsquared2steradians = 1.0/steradians2degsquared

  # LoTSS pointing table
  print (clustercatalogue)

  file = open('fieldsdict.pkl','rb')
  data = pickle.load(file)
  ID      = np.array([data[i]['id'] for i in range(len(data))])
  status  = np.array([data[i]['status'] for i in range(len(data))])
  RAp     = np.array([data[i]['ra'] for i in range(len(data))])
  DECp    = np.array([data[i]['decl'] for i in range(len(data))])
  DR3stat = np.array([data[i]['dr3'] for i in range(len(data))])

  idy = np.where( (DR3stat == 1) | (DR3stat == 2) )[0] # FINAL DR3 AREA; DR3stat == 1 >> the field will be included ; DR3stat == 2 >> final map; not all in the list because the mosaic is done only when also the neighbour fields are present
  #idy = np.arange(0,len(RAp)) # ALL POINTINGS

  # cluster tables
  data = fits.open(clustercatalogue)[1].data
  namelist  = np.array(data['Name'])
  RA        = np.array(data['RAJ2000'])
  DEC       = np.array(data['DEJ2000'])

  DECcut = 0. # cut in declination because of LoTSS pointing
  idx = np.where( DEC >= DECcut )[0] 


  # cross match cluster list and LoTSS-DR3 pointings
  allmatchID, allstatusID, allnoiseID, noiseavgID, expnoiseID, allsepdeg = ['']*len(RA), ['']*len(RA), [-1]*len(RA), [-1]*len(RA), [-1]*len(RA), ['']*len(RA)
  for i in idx: # LoTSS-DR3 clusters
    matchID, statusID, noiseID, sepdeg = [], [], [], []
    for j in idy:
      if sepn(RA[i]*deg2rad, DEC[i]*deg2rad, RAp[j]*deg2rad, DECp[j]*deg2rad)*rad2deg <= 2.2 :
        #print (namelist[i])
        matchID   = np.append(matchID, ID[j])
        statusID  = np.append(statusID, convertDR3status(DR3stat[j]))
        sepdeg    = np.append(sepdeg, round(sepn(RA[i]*deg2rad, DEC[i]*deg2rad, RAp[j]*deg2rad, DECp[j]*deg2rad)*rad2deg, 1))

        try:
          noiseID   = np.append(noiseID, map_noise(str(ID[j]), RA[i], DEC[i])*1e6 )
        except:
          #noiseID   = np.append(noiseID, analytical_noise(DEC[i]*deg2rad)*1e6 ) # this needs to go out as soon as all the fields will be in the final mosaic 
          noiseID   = np.append(noiseID, np.nan ) # this needs to go out as soon as all the fields will be in the final mosaic 

    allmatchID[i]  = ' '.join(map(str,matchID))
    allstatusID[i] = ' '.join(map(str,statusID))
    expnoiseID[i]  = round(analytical_noise(DEC[i]*deg2rad)*1e6,0)
    allnoiseID[i]  = ' '.join(format(x, ".1f") for x in noiseID)
    noiseavgID[i]  = round(np.nanmean(noiseID), 1)
    allsepdeg[i]   = ' '.join(map(str,sepdeg))

  # write the table with the LoTSS pointing
  newclustercatalogue = clustercatalogue.replace('.fits','matched.fits')

  table = Table.read(clustercatalogue)
  table['POINTING_ID'] = allmatchID
  table['Status']      = allstatusID
  table['separation']  = allsepdeg
  table['noise_ID']    = allnoiseID
  table['noise']       = noiseavgID
  table['noise_DEC']   = expnoiseID

  table.write(newclustercatalogue, overwrite=True)
  
  table = Table.read(newclustercatalogue)
  try:
    ids = np.where( (table['noise'] != -1.) & (np.isfinite(table['M500'])) )[0]
  except:
    ids = np.where((table['noise'] != -1.))[0]
  newtable = table[ids]
  newtable.write(newclustercatalogue, overwrite=True)

if True:
  match_LoTSSpointingID_to_catalogue('./cluster_catalogues/ACT-DR5.fits')
  match_LoTSSpointingID_to_catalogue('./cluster_catalogues/MCXC2.fits')
  match_LoTSSpointingID_to_catalogue('./cluster_catalogues/PSZ2.fits')
  match_LoTSSpointingID_to_catalogue('./cluster_catalogues/eROSITA-GE.fits')

