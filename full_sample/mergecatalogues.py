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
import argparse
from astropy.io import fits
from astropy.table import Table


def prepare_catalogues(clustercatalogue, DECmin, DECmax):  
  print (clustercatalogue, DECmin, DECmax)

  # cluster tables
  table = Table.read(clustercatalogue)
  namelist  = np.array(table['Name'])
  z         = np.array(table['z'])
  RA        = np.array(table['RAJ2000'])
  DEC       = np.array(table['DEJ2000'])
  M500      = np.array(table['M500'])

  #DECcut < 0. # cut in declination because of LoTSS pointing
  idx = np.where( (DEC >= DECmin) & (DEC <= DECmax) )[0] 

  # cross match cluster list and LoTSS-DR3 pointings
  R500kpc, R500amin = [-1]*len(RA), [-1]*len(RA)
  
  for i in idx: # LoTSS-DR3 clusters
    if (z[i] != -1) and (M500[i] > 0.):
      R500kpc[i]     = round(radius(M500[i], z[i], rho500=True)[0])
      R500amin[i]    = round(radius(M500[i], z[i], rho500=True)[1],3)
  
  newclustercatalogue = clustercatalogue.replace('.fits','matched.fits')
  # write the table with the LoTSS pointing
  table = Table.read(clustercatalogue)
  table['R500kpc']     = R500kpc
  table['R500amin']    = R500amin

  # add error on redshift; if null, use 0.01*(1+z)
  try:
    table['e_z'].fill_value = 0
    ids1 = np.where( (table['e_z'] == 0.) | (table['e_z'].filled() == 0.) )[0]
  except:
    ids1 = np.where( table['e_z'] == 0. )[0]
  table['e_z'][ids1] = 0.01 * (1+table['z'][ids1])

  table.write(newclustercatalogue, overwrite=True)
  
  #remove bad lines
  table = Table.read(newclustercatalogue)
  ids2 = np.where( (table['R500kpc'] != -1 ) )[0]
  newtable = table[ids2]
  newtable.write(newclustercatalogue, overwrite=True)


def radius(M500, z, rho500=True):
  from astropy.cosmology import FlatLambdaCDM
  cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
  
  rhoc = cosmo.critical_density(z)
  M500 *= 1.e14*(2.e33) #M500 in g
  R500 = ( (3.*M500)/(4.*np.pi* (500*rhoc.value)) )**(1./3.) #R500 in cm
  R500kpc  = R500/3.08e21 #R500 in kpc
  R500amin = (R500/3.08e21 / cosmo.kpc_proper_per_arcmin(z).value) #R500 in arcmin   

  if rho500:
    Rkpc = R500kpc ; Ramin = R500amin
  else:
    Rkpc = R500kpc/0.7 ; Ramin = R500amin/0.7

  return (Rkpc, Ramin)


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

parser = argparse.ArgumentParser(description='cross-match a catalogue with other(s)')
parser.add_argument('-c','--cataloglist', nargs='+', help='Catalog(s) to check for the crossmatch', required=True, type=str)
parser.add_argument('-o', '--outputname', help='output name of the cross-matched catalogue', default='crossmatched.fits', required=False, type=str)
parser.add_argument('--DEClims', nargs=2, help='limits in declination', required=True, type=float)
parser.add_argument('--outputcolumn', help='output name of the cross-matched column', default='crossmatched', required=False, type=str)
parser.add_argument('--onlyclusters', help='cut to M>1e14 Msun, so to have only clusters',action='store_true')
args = parser.parse_args()

cataloglist   = args.cataloglist
outputcatalog = args.outputname
declims       = args.DEClims

print (cataloglist)

for catalog in cataloglist:
  #print (catalog)
  if not os.path.exists(catalog.replace('.fits','matched.fits')):
    print ('calculate R500/R500amin and add it to the table')
    prepare_catalogues(catalog, DECmin=declims[0], DECmax=declims[1])


if not os.path.exists('all.fits'):
  newname, newRA, newDEC, newz, newez, newM500, newRkpc, newRamin = \
  np.array([]),  np.array([]), np.array([]),  np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
  
  for catalog in cataloglist:
    #print (catalog.replace('.fits','matched.fits'))
    data  = fits.open(catalog.replace('.fits','matched.fits'))[1].data
    name      = data['Name']
    RA        = data['RAJ2000']
    DEC       = data['DEJ2000']
    z         = data['z']
    e_z       = data['e_z']
    M500      = data['M500']
    R500kpc   = data['R500kpc']
    R500amin  = data['R500amin']

    newname = np.append(newname, name)
    newRA   = np.append(newRA, RA)
    newDEC  = np.append(newDEC, DEC)
    newz    = np.append(newz, z)
    newez   = np.append(newez, e_z)
    newM500 = np.append(newM500, M500)
    newRkpc = np.append(newRkpc, R500kpc)
    newRamin= np.append(newRamin, R500amin)

  newdata = Table()
  newdata['Name']        = np.array(newname, dtype=np.str)
  newdata['RAJ2000']     = np.array(newRA, dtype=np.float64)
  newdata['DEJ2000']     = np.array(newDEC, dtype=np.float64)
  newdata['z']           = np.array(newz, dtype=np.float64)
  newdata['e_z']         = np.array(newez, dtype=np.float64)
  newdata['M500']        = np.array(newM500, dtype=np.float64)
  newdata['R500kpc']     = np.array(newRkpc, dtype=np.float64)
  newdata['R500amin']    = np.array(newRamin, dtype=np.float64)

  newdata.write('all.fits', format='fits', overwrite=True)

  # add cross-match from other catalogs
  table = Table.read('all.fits')
  names  = table['Name']
  ras    = table['RAJ2000']
  decs   = table['DEJ2000']
  zs     = table['z']
  ezs    = table['e_z']
  ramin  = table['R500amin']  

  ALLaltname = ['']*len(table['Name'])
  for i in range(len(table['Name'])):
    altname, multindex = [], []
    for j in range(len(table['Name'])):
      if (i != j) and (abs(zs[i]-zs[j]) <= np.sqrt(ezs[i]**2 + ezs[j]**2)):
        d = sepn(ras[i]*deg2rad, decs[i]*deg2rad, ras[j]*deg2rad, decs[j]*deg2rad)*rad2arcmin
        if (d < 0.5*ramin[i] or d <  0.5*ramin[j]):
          altname   = np.append(altname, table['Name'][j])
          multindex = np.append(multindex, j)
    
    ALLaltname[i] = ', '.join(map(str,altname))

  table[args.outputcolumn] = ALLaltname
  table.write('all.fits', overwrite=True)

dataref  = fits.open('all.fits')[1].data
nameref  = dataref['Name'] 
raref    = dataref['RAJ2000']
decref   = dataref['DEJ2000']
altname  = dataref[args.outputcolumn]

k = len(fits.open(cataloglist[0].replace('.fits','matched.fits'))[1].data) # number of lines in the first catalog
multline = np.array([])
for i in range(1, len(cataloglist)):
  print (cataloglist[i].replace('.fits','matched.fits'))
  data  = fits.open(cataloglist[i].replace('.fits','matched.fits'))[1].data
  name  = data['Name']

  ids = np.where( [any(name[j] in alt for alt in altname[0:k]) for j in range(len(name))] )[0]  
  multline = list(map(int, np.append(multline, ids+k)))
  k+=len(name)

#print (multline)
#sys.exit()
table = Table.read('all.fits')
table.remove_rows(multline)
table.write(outputcatalog, overwrite=True)

if args.onlyclusters:
  table = Table.read(outputcatalog)
  ids = np.where( table['M500'] >= 1. )[0]
  newtable = table[ids]
  newtable.write(outputcatalog, overwrite=True)
