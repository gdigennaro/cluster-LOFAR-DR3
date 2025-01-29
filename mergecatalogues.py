# SCRIPT TO CROSS-MATCH CLUSTERS IN DIFFERENT CATALOGUES
# 
# G. Di Gennaro
# October 2024

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import sys, os, glob
import itertools
#import pyfits
from astropy.io import fits
from astropy.table import Table

def sepn(r1,d1,r2,d2):
  """
  Calculate the separation between 2 sources, RA and Dec must be
  given in radians. Returns the separation in radians [TWS]
  """
  cos_sepn = np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)
  sepn = np.arccos(cos_sepn)
  return sepn   

#Define various angle conversion factors (multiply to undertake operation)
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


cataloglist = [
  './cluster_catalogues/PSZ2matched.fits',
  './cluster_catalogues/ACT-DR5matched.fits',
  './cluster_catalogues/MCXC2matched.fits',
  './cluster_catalogues/eROSITA-GEmatched.fits'
]

allcatalog = './cluster_catalogues/ALLmatched.fits'

if True:
  newname, newRA, newDEC, newz, newM500, newID, newdist, newsig = np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
  for catalog in cataloglist:
    print (catalog)
    data  = fits.open(catalog)[1].data
    name  = data['Name']
    RA    = data['RAJ2000']
    DEC   = data['DEJ2000']
    z     = data['z']
    P_ID  = data['POINTING_ID']
    dP    = data['separation']
    noise = data['noise']
    try:
      #Y500  = data['Y5R500']
      #M     = data['MSZ']
      #print (data['z'][ids])
      ids = []
      ids = np.where((data['MSZ'] != 0.))[0]
      #M[ids] = MassPSZ(z[ids], Y500[ids]*1E-3)/1e14
      M     = data['MSZ'][ids]
    except:
      pass
    try:
      ids = []
      ids = np.where((data['M500cC'] != 0.))[0]
      M = data['M500cC'][ids]
    except:
      pass
    try:
      ids = []
      ids = np.where((data['M500'] != 0.))[0]
      M = data['M500']
    except:
      pass  
    
    newname = np.append(newname, name[ids])
    newRA   = np.append(newRA, RA[ids])
    newDEC  = np.append(newDEC, DEC[ids])
    newz    = np.append(newz, z[ids])
    newM500 = np.append(newM500, M[ids])
    newID   = np.append(newID, P_ID[ids])
    newdist = np.append(newdist, dP[ids])
    newsig  = np.append(newsig, noise[ids])


  c1 = fits.Column(name='Name', array=newname, format='30A')
  c2 = fits.Column(name='z', array=newz, format='F8.5')
  c3 = fits.Column(name='RAJ2000', array=newRA, format='F8.4')
  c4 = fits.Column(name='DEJ2000', array=newDEC, format='F8.4')
  c5 = fits.Column(name='M500', array=newM500, format='F5.2')
  c6 = fits.Column(name='POINTING_ID', array=newID, format='40A')
  c7 = fits.Column(name='distance', array=newdist, format='40A')  
  c8 = fits.Column(name='noise', array=newsig, format='D')
  t = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8])
  t.writeto(allcatalog, overwrite=True)

  # add cross-match from other catalogs
  data  = fits.open(allcatalog)[1].data
  ras    = data['RAJ2000']
  decs   = data['DEJ2000']
  zs     = data['z']

  ALLaltname = ['']*len(data['Name'])

  for i in range(len(data['Name'])):
    altname, multindex = [], []
    for j in range(len(data['Name'])):
      if sepn(ras[i]*deg2rad, decs[i]*deg2rad, ras[j]*deg2rad, decs[j]*deg2rad)*rad2arcmin <= 1. and i != j and abs(zs[i]-zs[j]) <= 0.01 :
        altname   = np.append(altname, data['Name'][j])
        multindex = np.append(multindex, j)
    
    ALLaltname[i] = ', '.join(map(str,altname))

  table = Table.read(allcatalog)
  table['Alt.Name'] = ALLaltname
  table.write(allcatalog, overwrite=True)


#reftable = allcatalog
dataref  = fits.open(allcatalog)[1].data
nameref  = dataref['Name'] 
raref    = dataref['RAJ2000']
decref   = dataref['DEJ2000']
altname  = dataref['Alt.Name']

k = len(fits.open('./cluster_catalogues/PSZ2matched.fits')[1].data) # number of lines in PSZ2
multline = np.array([])
for i in range(1, len(cataloglist)):
  print (cataloglist[i])
  data  = fits.open(cataloglist[i])[1].data
  name  = data['Name']

  ids = np.where( [any(name[j] in alt for alt in altname[0:k]) for j in range(len(name))] )[0]
  #print (ids)
  
  multline = list(map(int, np.append(multline, ids+k)))
  k+=len(name)


#print (multline)
table = Table.read(allcatalog)
table.remove_rows(multline)
table.write(allcatalog.replace("matched","crossmatched"), overwrite=True)
