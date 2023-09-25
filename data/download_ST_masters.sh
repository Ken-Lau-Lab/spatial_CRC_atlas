# download ST_master/ folder containing fully annotated AnnData objects
# for all Visium samples, and move them to the data/ST/ directory
wget <ST_master H5AD Dropbox URL>

mv ST_master/*.h5ad ST/  # move files to ST/ dir
rm -r ST_master/
