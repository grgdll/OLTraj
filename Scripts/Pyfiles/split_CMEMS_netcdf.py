# Script to split the CMEMS global files into daily sanpshots
import xarray as xr
import pandas as pd
import os
import glob
import argparse
# Debug
import ipdb


def main(year):
    # Get list of CMEMS files to process
    pathin = '../../Data/'
    #fname = 'dataset-duacs-rep-global-merged-allsat-phy-l4*.nc'
    fname = str(year)+'_???_DATASET-DUACS-REP-GLOBAL-MERGED-ALLSAT-PHY-L4.nc'
    flist = glob.glob(os.path.join(pathin,fname))
    flist = sorted(flist)
    # Define out path
    # Keep CMEMS velocities in Data and only trajectories in Processed!!!!
    # Create subfolder to distinguinsh original product and daily files
    pathout = os.path.join(pathin, 'Daily')

    # Create folder if it does not exists
    if not os.path.isdir(pathout):
        os.makedirs(pathout)
    # Read each CMEMS file
    for ifile in flist:
        ds = xr.open_dataset(ifile)
        # Process each time step
        for i,itime in enumerate(pd.to_datetime(ds.time.values)):
            print("i=",i, " itime=",itime)
            # Create _pathout subdirectory (by year)
            _pathout = os.path.join(pathout,itime.strftime('%Y'))
            # Create subdirectory if it does not exist
            if not os.path.isdir(_pathout):
                os.makedirs(_pathout)
            # Create out file name
            fnameout = '%s%s.nc' % (itime.strftime('%Y%m%d'),os.path.basename(ifile).split('_')[0][8:])
            fnameout = os.path.join(_pathout,fnameout)
            # Process only if file do not exists
            if not os.path.isfile(fnameout):
                print('Save %s' % os.path.basename(fnameout))
                # Get only specific timestep
                _ds = ds.isel(time=[i])
                # Remove conflicting attributes to avoid
                # AttributeError: NetCDF: String match to name in use
                try:
                    del _ds.attrs['_NCProperties']
                    #del _ds.time.attrs['CLASS']
                    #del _ds.time.attrs['NAME']
                except KeyError:
                    pass
                # Save to netcdf
                _ds.to_netcdf(fnameout)
            else:
                print('File %s already created!' % os.path.basename(fnameout))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', required=True, help='Year to process')
    args = parser.parse_args()
    #ftag = int(args.ftag)
    year = int(args.year)
    #main(ftag,year)
    main(year)
    #main()
