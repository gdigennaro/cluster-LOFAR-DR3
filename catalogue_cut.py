# SCRIPT TO CHECK MAKE A SUBCATALOGUE GIVEN CLUSTER MASS/REDSHIFT/LoTSS NOISE
# 
# G. Di Gennaro
# Sept 2024

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import sys, os, glob
#import pyfits
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib import gridspec
from matplotlib import cm as color_map
from matplotlib.colorbar import Colorbar
from palettable.cartocolors.diverging import Fall_4 as colmap
from palettable.matplotlib import Magma_20 as colmap2


import matplotlib.lines as mlines
import matplotlib as mpl
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction' ]= 'in'
mpl.rcParams['errorbar.capsize'] = 3
mpl.rcParams['legend.handlelength'] = 1.5

def subcatalogues(clustercatalogue, Mcut, Ncut, zcut):
  # LoTSS pointing table
  print (clustercatalogue)
  data  = fits.open(clustercatalogue)[1].data
  noise = np.array(data['noise'])
  z     = np.array(data['z'])
  try:
    M  = np.array(data['MSZ'])
  except:
    pass
  try:
    M   = np.array(data['M500cC'])
  except:
    pass
  try:
    M   = np.array(data['M500'])
  except:
    pass  
    
  cuts = ''
  if Mcut: cuts += 'M'
  if zcut: cuts += 'z'
  if Ncut: cuts += 's'  

  if cuts == 'Mzs':
    ids = np.where( (M >= Mcut)  & (M != '') & (z >= zcut) & (noise <= Ncut) )[0] 
  elif cuts == 'Mz':
    ids = np.where( (M >= Mcut)  & (M != '') & (z >= zcut) )[0] 
  elif cuts == 'zs':
    ids = np.where( (z >= zcut) & (noise <= Ncut) )[0]
  elif cuts == 'M':
    ids = np.where( (M >= Mcut)  & (M != '') )[0]
  elif cuts == 'z':
    ids = np.where( (z >= zcut) )[0] 
  elif cuts == 's':
    ids = np.where( (noise <= Ncut) )[0]
  print (Mcut, Ncut, zcut, len(ids))

  if True:
    # write the subtable
    newclustercatalogue = clustercatalogue.replace('matched.fits','matched_'+cuts+'cut.fits')
    table = Table.read(clustercatalogue)
    newtable = table[ids]
    newtable.write(newclustercatalogue, overwrite=True)

  if False:
    # plot cut catalog 
    name = clustercatalogue.split("matched")[0].split('./cluster_catalogues/')[1]

    fig = plt.figure(figsize=(8.5,6))
    gs = gridspec.GridSpec(nrows=2, ncols=4, width_ratios=[1,0.15,0.01,0.04], height_ratios=[0.2,1.0], wspace=0., hspace=0.)
    ax1 = fig.add_subplot(gs[1,0]) # main
    ax2 = fig.add_subplot(gs[1,1]) # histogram vertical
    ax3 = fig.add_subplot(gs[0,0]) # histogram horizontal
    ax4 = fig.add_subplot(gs[1,3]) # colorbar

    ax1.tick_params(axis='both',which='both',labelsize=12,right=True, top=True)
    ax2.tick_params(axis='both',which='both',labelsize=12,right=True, top=True)
    ax3.tick_params(axis='both',which='both',labelsize=12,right=True, top=True)

    plt.suptitle(name+" clusters in LoTSS-DR3\n"+r"$M\geq%s\times10^{14}\,{\rm M_\odot}, ~ z\geq%s, ~ \sigma_{\rm rms}\leq%s\,{\rm\mu Jy\,beam^{-1}}$"%(Mcut,zcut,Ncut), fontsize=16)
    #plt.suptitle(name+" clusters in LoTSS-DR3", fontsize=16)
    #ax3.set_title(r"$M\geq%s\times10^{14}\,{\rm M_\odot}, ~ z\geq%s, ~ \sigma_{\rm rms}\leq%s\,{\rm\mu Jy\,beam^{-1}}$"%(Mcut,zcut,Ncut), fontsize=14)    
    #ax3.set_title(name+" clusters in LoTSS-DR3\n"+r"$M\geq%s\times10^{14}\,{\rm M_\odot}, ~ z\geq%s, ~ \sigma_{\rm rms}\leq%s\,{\rm\mu Jy\,beam^{-1}}$"%(Mcut,zcut,Ncut), fontsize=16)
    ax1.set_xlabel(r"redshift $[z]$", fontsize=14)
    ax1.set_ylabel(r"cluster mass $[M_{500}~\rm (\times10^{14}~M_{\odot})]$", fontsize=14)  
    ax1.grid(color='lightgrey', linestyle=':', which='both', linewidth=0.7)
    ax2.grid(color='lightgrey', linestyle=':', which='both', linewidth=0.7)
    ax3.grid(color='lightgrey', linestyle=':', which='both', linewidth=0.7)

    ax1.set_ylim(M[ids].min()-0.1*M[ids].min(),M[ids].max()+0.1*M[ids].max())
    ax2.set_ylim(M[ids].min()-0.1*M[ids].min(),M[ids].max()+0.1*M[ids].max())
    ax1.set_xlim(z[ids].min()-0.1,z[ids].max()+0.1)
    ax3.set_xlim(z[ids].min()-0.1,z[ids].max()+0.1)

    ax1.set_yscale('log')
    ax2.set_yscale('log')

    plt1 = ax1.scatter(z[ids], M[ids], s=50, c=noise[ids], cmap=colmap2.mpl_colormap, vmin=62, vmax=np.max(noise[ids]), edgecolors='grey', linewidth=0.3, alpha=0.8)
    #ax1.plot(z[ids], M[ids], 'o', color='k', alpha=0.6, label=name+" (%s clusters)"%str(len(z)))
    
    xbins = np.arange(0,z[ids].max()+0.1, 0.1)
    ybins = 10**np.linspace(np.log10(M[ids].min()), np.log10(M[ids].max()), 20)    
    ax2.hist(M[ids], bins=ybins, density=False, histtype='stepfilled',orientation='horizontal', color='grey', alpha=0.35)
    ax3.hist(z[ids], bins=xbins, density=False, histtype='stepfilled',orientation='vertical', color='grey',alpha=0.35)

    ax2.axhline(np.median(M[ids]), color='grey', linewidth=2, linestyle='--',zorder=1)
    ax3.axvline(np.median(z[ids]), color='grey', linewidth=2, linestyle='--',zorder=1)

    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax1.yaxis.set_minor_formatter(mtick.FormatStrFormatter('%d'))
    ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))
    ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%d'))

    ax2.tick_params(axis='x',labelsize=9)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax2.get_yminorticklabels(), visible=False)

    ax3.tick_params(axis='y',labelsize=9)
    plt.setp(ax3.get_xticklabels(), visible=False)

    cbar = Colorbar(mappable=plt1, ax=ax4)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(r'expected noise map $(\rm \sigma_{rms}~[\mu Jy\,beam^{-1}])$', fontsize=12)

    marker  = mlines.Line2D([], [], color='grey', marker='o', linestyle='None', label="%s clusters"%str(len(z[ids])))
    ax1.legend(handles=[marker], handletextpad=0.1, prop={'size': 10}, frameon=True)


    plt.tight_layout()
    plt.savefig('./images/'+name+'_LoTSS-DR3_mass_z_subsample.pdf')
    #plt.show()
    plt.close()
    #sys.exit()  

Mcut = 4.0
Ncut = 100.
zcut = 0.6

subcatalogues('./cluster_catalogues/ACT-DR5matched.fits', Mcut, Ncut, zcut)
subcatalogues('./cluster_catalogues/MCXC2matched.fits', Mcut, Ncut, zcut)
subcatalogues('./cluster_catalogues/PSZ2matched.fits', Mcut=4, Ncut=200, zcut=0.6)
subcatalogues('./cluster_catalogues/eROSITA-GEmatched.fits', Mcut, Ncut, zcut)
#subcatalogues('./cluster_catalogues/DESI-WHmatched.fits', Mcut, Ncut, zcut)

