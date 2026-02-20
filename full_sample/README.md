Here we report a metacatalogue of clusters within the LoTSS-DR3 footprint. 
These codes were used in Stuardi+2026 to build an ICM-selected metacatalogue of clusters

- `check_cluster_inLoTSSfields.py` associates all the clusters in a given sample to the LoTSS-DR3 pointing. CLUSTERCATALOG.fits becomes CLUSTERCATALOGmatched.fits; it uses `fieldsdict.fits` to get LoTSS-DR3 fields information (coordinates and quality based on Stuardi+26)
- `mergecatalogues.py` double checks multiple cluster catalogues for cluster duplicates; see 
- `subcatalogue.py` gives sub-catalogues based on cuts in redshift/mass/noise (multiple choice)
