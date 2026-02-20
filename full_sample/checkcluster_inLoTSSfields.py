# SCRIPT TO CHECK WHETHER A CLUSTER WITH GIVEN (RA,DEC) POSITION IS OBSERVED BY LOTSS
# 
# G. Di Gennaro
# Sept 2024

import warnings
warnings.filterwarnings('ignore')

import pickle
import pandas as pd
import numpy as np
import sys, os, glob
#import pyfits
import argparse
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


def convertDR3status(status):
  if status == 1: plan = 'Progress'
  if status == 2: plan = 'Final'
  if status == 'None': plan = 'None'

  return (plan)


def radius(M500, z, rho500=True):
  from astropy.cosmology import FlatLambdaCDM
  cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
  
  rhoc = cosmo.critical_density(z)
  M500 *= 1.e14*(2.e33) #M500 in g
  R500 = ( (3.*M500)/(4.*np.pi* (500*rhoc.value)) )**(1./3.) #R500 in cm
  R500kpc  = R500/3.08e21 #R500 in kpc
  R500amin = round(R500/3.08e21 / cosmo.kpc_proper_per_arcmin(z).value,3) #R500 in arcmin   

  if rho500:
    Rkpc = R500kpc ; Ramin = R500amin
  else:
    Rkpc = R500kpc/0.7 ; Ramin = R500amin/0.7

  return (Rkpc, Ramin)
  


def match_LoTSSpointingID_to_catalogue(clustercatalogues, pointinglist ,dcoord=2.2, DECcut=0):
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

  
  for clustercatalogue in clustercatalogues:
    # LoTSS pointing table
    print (clustercatalogue)

    datalotss = Table.read(pointinglist)
    #print (datalotss['POINTING_ID']) ; sys.exit()
    ID    = np.array([datalotss[i]['id'].strip() for i in range(len(datalotss))])
    RAp   = np.array([datalotss[i]['RAJ2000'] for i in range(len(datalotss))])
    DECp  = np.array([datalotss[i]['DEJ2000'] for i in range(len(datalotss))])
    flag  = np.array([int(datalotss[i]['flag']) for i in range(len(datalotss))])
    idy   = np.arange(0,len(ID),1)
    
    # cluster tables
    table = Table.read(clustercatalogue)
    namelist  = np.array(table['Name'])
    z         = np.array(table['z'])
    RA        = np.array(table['RAJ2000'])
    DEC       = np.array(table['DEJ2000'])
    M500      = np.array(table['M500'])

    #DECcut = 0. # cut in declination because of LoTSS pointing
    idx = np.where( (DEC >= DEClims[0]) & (DEC <= DEClims[1])  )[0] 

    # cross match cluster list and LoTSS-DR3 pointings
    R500kpc, R500amin = [-1]*len(RA), [-1]*len(RA)
    
    for i in idx: # LoTSS-DR3 clusters
      #print (namelist[i])
      
      matchID, flagID, sepdeg = [], [], []
      for j in idy:
        if sepn(RA[i]*deg2rad, DEC[i]*deg2rad, float(RAp[j])*deg2rad, float(DECp[j])*deg2rad)*rad2deg <= 2.2 :
          matchID   = np.append(matchID, ID[j])
          flagID  = np.append(flagID, '{:d}'.format(flag[j]))
          sepdeg    = np.append(sepdeg, round(sepn(RA[i]*deg2rad, DEC[i]*deg2rad, float(RAp[j])*deg2rad, float(DECp[j])*deg2rad)*rad2deg, 2))

      allmatchID[i]  = ' '.join(map(str,matchID))
      allflagID[i]   = ' '.join(map(str,flagID))
      allsepdeg[i]   = ' '.join(map(str,sepdeg))
      
      if (z[i] != -1) and (M500[i] > 0.):
        R500kpc[i]     = round(radius(M500[i], z[i], rho500=True)[0])
        R500amin[i]    = round(radius(M500[i], z[i], rho500=True)[1],3)
      else:
        R500kpc[i]     = -1
        R500amin[i]    = -1
    
    newclustercatalogue = clustercatalogue.replace('.fits','matched.fits')
    # write the table with the LoTSS pointing
    table = Table.read(clustercatalogue)
    table['POINTING_ID'] = allmatchID
    table['FLAG']        = allflagID
    table['separation']  = allsepdeg
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
    table['POINTING_ID'].fill_value = 'N/A'
    ids2 = np.where( (table['z'] != -1.) & (table['POINTING_ID'].filled() != 'N/A') )[0]
    newtable = table[ids2]
    newtable.write(newclustercatalogue, overwrite=True)
  

parser = argparse.ArgumentParser(description='Run extraction and selfcalibration of clusters in LoTSS; you can give either a catalog (FITS format) or the cluster name')
parser.add_argument('-c','--cataloglist', nargs='+', help='Catalog to use from which extract clusters; examples=ACT-DR5.fits', required=True, type=str)
parser.add_argument('--DEClims', nargs=2, help='limits in declination', required=True, type=int)
parser.add_argument('--pointinglist', help='Catalog to use from which extract clusters', default='fieldsdict.fits', required=False, type=str)
args = parser.parse_args()

cataloglist  = args.cataloglist
pointinglist = args.pointinglist
match_LoTSSpointingID_to_catalogue(cataloglist, pointinglist, args.DECcut)
