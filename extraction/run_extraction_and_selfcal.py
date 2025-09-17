import argparse
import glob, os, sys
import numpy as np
import pickle
try:
    import bdsf as bdsm
except ImportError:
    import lofar.bdsm as bdsm
from crossmatch_utils import *
import pyregion
import time
start_time = time.time()

from surveys_db import *
from db_utils import *

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import wcs
from astropy.wcs import WCS
from auxcodes import separator

sys.path.remove('/opt/lofar/ddf-pipeline/utils')
sys.path.insert(1, "/local/work/g.digennaro/software/extraction-utils/ddf-pipeline/utils")
sys.path.insert(1, "/local/work/g.digennaro/software/extraction-utils/lotss-hba-survey")
sys.path.remove('/opt/lofar/ddf-pipeline/scripts')
sys.path.insert(1, "/local/work/g.digennaro/software/extraction-utils/ddf-pipeline/scripts")

rclone = os.environ['RCLONE_CONFIG_DIR'] = '/local/work/g.digennaro/software/'
#print (rclone) ; sys.exit()

from reprocessing_utils import *

"""
Script to run extraction and automatic selfcal from LoTSS poining.

to run inside latest flocs (type flocs-setup, see bash_aliases)
"""

#DATADIR = '/local/work/g.digennaro/LoTSS-DR3/'
DATADIR = './'

def upload_extract(cname,uploadfilename):
	rclonepath = os.environ['RCLONE_CONFIG_DIR']
	print('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy %s maca_sksp_disk_extract:%s/selfcal/'%(rclonepath,uploadfilename,cname))
	os.system('rclone  --multi-thread-streams 1 --config=%s/maca_sksp_disk_extract.conf copy %s maca_sksp_disk_extract:%s/selfcal/'%(rclonepath,uploadfilename,cname))


parser = argparse.ArgumentParser(description='Run extraction and selfcalibration of clusters in LoTSS; you can give either a catalog (FITS format) or the cluster name')
parser.add_argument('--doextraction', help='set it True if you want to extract clusters from archive', action='store_true')
parser.add_argument('--doselfcal', help='set it True if you want to run facetselfcal', action='store_true')
parser.add_argument('--stopselfcal', help='set it True if you want to stop the selfcal to cycle 0', action='store_true')
parser.add_argument('--restartselfcal', help='set it True if you want to stop the selfcal to cycle 0', action='store_true')
parser.add_argument('--docopy', help='set it True if you want to copy the selfcal output on Surfsara', action='store_true')
parser.add_argument('-i','--clustername', help='cluster name, if you want to extract a single cluster', default='', required=False, type=str)
parser.add_argument('--RA', help='cluster RA (in deg)', required=False, type=float)
parser.add_argument('--DEC', help='cluster DEC (in deg)', required=False, type=float)
#parser.add_argument('--docatalog', help='set it True if you want to run the script for a catalog of clusters', action='store_true')
parser.add_argument('-c','--catalog', help='Catalog to use from which extract clusters', required=False, type=str)

args = vars(parser.parse_args())

if args['catalog'] and args['clustername']:
  print ("Error: either give a single target name or a cluster catalog")
  sys.exit()

elif args['catalog'] and not args['clustername']:
  print ("use catalog:", args['catalog'])
  data = fits.open(args['catalog'])[1].data
  clusterlist = np.array(data['Name']) 
  ra         = np.array(data['RAJ2000'])
  dec        = np.array(data['DEJ2000'])

elif args['clustername'] and not args['catalog']:
  clusterlist = [args['clustername']]
  ra         = [args['RA']]
  dec        = [args['DEC']]

else:
  print ("Error: give a single target name or a cluster catalog.")
  sys.exit()

for i, cluster in enumerate(clusterlist):
  print (os.getcwd())
  name = cluster.replace(' ','')
  print ("CLUSTER:", name)
  

  if args['doextraction']:
    #cmd = 'python /local/work/g.digennaro/software/extraction-utils/ddf-pipeline/scripts/optimise-extract-region.py --ra %s --dec %s --size 0.4 --region_file ./%s/%s.extraction.reg'%(ra[i], dec[i],name,name)
    #print (cmd) 
    #os.system(cmd)
    
    if not glob.glob(DATADIR+name+"/"+"P???+??*.dysco.sub.shift.avg.weights.ms.archive*"):
      
      # this is for the final extraction  
      #cmd = 'python /local/work/g.digennaro/software/extraction-utils/extraction-utils/ddf-pipeline/scripts/run_extraction_pipeline.py %s' %name
      #print (cmd)
      #os.system(cmd)

      cmd = 'python /local/work/g.digennaro/software/extraction-utils/ddf-pipeline/scripts/extraction.py %s' %name
      print (cmd)
      os.system(cmd)
    
    else:
      print ("Extraction already done")
     

  else:
    print ("No extraction requested for %s"%name)
  
  if os.path.exists(DATADIR+name) and args['doselfcal'] and len(glob.glob(DATADIR+name+"/"+"P???+??*.dysco.sub.shift.avg.weights.ms.archive*")) > 0:
    if not os.path.exists(DATADIR+name+'/'+name+'.ds9.tar.gz'):      
      # RUN SELF CALIBRATION AND COPYING IN THE SURFSARA DATABASE
      CLUSTERDIR = DATADIR+name+'/'
      os.chdir(CLUSTERDIR)
      print (os.getcwd())
      
      # this is for the final selfcal  -- maybe not
      #cmd = 'python /local/work/g.digennaro/software/extraction-utils/extraction-utils/ddf-pipeline/scripts/run_selfcal_pipeline.py %s' %( str(name) )
      #print (cmd)
      #os.system(cmd)

      if False: # if args['docopy']:
        # update the database to give success
        selfcal_status = 'STARTED'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['selfcal_status'] = selfcal_status
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print('Updated status to STARTED for',name)

      cmd  = 'python /local/work/g.digennaro/software/extraction-utils/lofar_facet_selfcal/facetselfcal.py '
      cmd += '--helperscriptspath="/local/work/g.digennaro/software/extraction-utils/lofar_facet_selfcal" '
      cmd += '--helperscriptspathh5merge="/local/work/g.digennaro/software/extraction-utils/lofar_helpers" '
      cmd += '-b %s.ds9.reg '%name
      cmd += '--remove-flagged-from-startend --auto '
      if args['stopselfcal']:
        cmd += '--stop 1 ' #this is only to check for eventually bad antennas
      if args['restartselfcal']:
        cmd += '--start 0 --stop 10' #this is to restart the selfcal process (after checking the bad antennas)      
      cmd += '-i %s *dysco.sub.shift.avg.weights.ms.archive?'%name
      
      print (cmd)
      os.system(cmd)


      # copy to SURF
      if args['docopy']:
        print('Archiving the results to SURF')
        f = glob.glob('%s.ds9.tar.gz'%name) + glob.glob('%s_009.png'%name)

        for uploadfilename in f:
          upload_extract(name,uploadfilename)
        
        # update the database to give success
        selfcal_status = 'SDONE'
        sdb=SurveysDB()
        extractdict = sdb.get_reprocessing(name)
        extractdict['selfcal_status'] = selfcal_status
        sdb.db_set('reprocessing',extractdict)
        sdb.close()
        print('Updated status to SDONE for',name)
    
    
      os.chdir('/local/work/g.digennaro/LoTSS-DR3/')  

    else:
      print ("Selfacl already done")   

  else:
    print ("check extraction output") 
    
  timerun = round((time.time() - start_time)/3600,1)
  separator("%s hr"%(str(timerun)))



