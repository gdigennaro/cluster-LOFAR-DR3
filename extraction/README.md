The script `run_extraction_and_selfcal.py` works both on single clusters, with coordinates cluster name provided, and with a FITS table with a list of clusters; in this latter case, the table needs to have Name, RAJ2000, DEJ2000 columns). This is made of two independent steps, the extraction (where the pointing is downloaded and sources outside a given region are subtracted) and the selfcal is performed using facetselfcal v17.14 in auto mode (4 cycles of phase only and 6 cycles of amp+phase calibration). 
It needs to be ran inside `flocs v6.1.0` (https://public.spider.surfsara.nl/project/lofarvwf/fsweijen/containers/)

Examples of how to run the script, for a cluster list:

`python run_extraction_and_selfcal.py -c clusterlist.fits [–-doextraction] [--doselfcal]`

and for a single cluster:

`python run_extraction_and_selfcal.py -i clustername --RA clusterRA --DEC clusterDEC [--size 0.5] [–-doextraction] [--doselfcal]`
