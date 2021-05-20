# For a given day combine bkw and fwd trajectory mat files
# into a single netcdf file

# run in parallel using this (note 2 arguments for xargs):
#echo {1998..2019..1}" "{1..12..1}| xargs -n2 | xargs -P7 -I{} sh -c 'V="{}"; python combine_traj.py --year=${V% *} --month=${V#* }'

import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import argparse
import os
import glob
# Custom
import lagrangian_plot_utils as lagr_utils
# Debug
import ipdb


# def combine_traj(day,model_res):
def combine_traj(year,month):
    # Define input path of .mat files 
    # (outdir is not used)
    pathin, outdir = lagr_utils.set_path()
    # Find files
    flist = glob.glob(os.path.join(pathin,'traj_aviso_%04d%02d*fwd.mat' % (year,month)))
    flist = sorted(flist)
    # Cycle through all files
    for ii, fwdfile in enumerate(flist):
        day = fwdfile.split('_fwd')[0][-8:]
        print('Process %s' %  day)
        #-----------------------
        # First fwd
        print('\tfwd')
        # Keep only filename
        fwdfile = os.path.basename(fwdfile)
        adv_aviso = lagr_utils.read_traj(fwdfile)
        # Particles positions with time
        lons = adv_aviso['lons']
        lats = adv_aviso['lats']
        # Final day
        day0 = adv_aviso['day0']
        dayf = adv_aviso['dayf']
        Nstep = adv_aviso['Nstep']
        dates = lagr_utils.get_dates(day0,dayf,Nstep)
        # Subsample dates to daily resolution
        time = np.asarray(dates[::4])
        # 2D matrices with initial particle position
        inilon = adv_aviso['long']
        inilat = adv_aviso['latg']
        lons2d = np.int16(np.ones((lons.shape[0],)+inilon.shape)*32767)
        lats2d = np.int16(np.ones((lons.shape[0],)+inilon.shape)*32767)
        # Find indices where inilon and inilat are not nans
        # (transpose to have F order: first index changing fastest, and the last index changing slowest)
        [ilon,ilat]= np.nonzero(~(np.isnan(inilon) & np.isnan(inilat)).T)
        # Assign values
        lons2d[:,ilat,ilon] = lons
        lats2d[:,ilat,ilon] = lats
        # Convert to mask array
        trajlon = np.ma.array(lons2d, mask=(lons2d==32767))
        trajlat = np.ma.array(lats2d, mask=(lats2d==32767))
        #-----------------------

        #-----------------------
        # Then bkw
        print('\tbwd')
        bkwfile = fwdfile.replace('fwd','bkw')
        adv_aviso = lagr_utils.read_traj(bkwfile)
        # Particles positions with time
        lons = adv_aviso['lons']
        lats = adv_aviso['lats']
        # Final day
        day0 = adv_aviso['day0']
        dayf = adv_aviso['dayf']
        Nstep = adv_aviso['Nstep']
        dates = lagr_utils.get_dates(day0,dayf,Nstep)
        # Subsample dates to daily resolution
        _time = np.asarray(dates[::4])
        # 2D matrices with initial particle position
        inilon = adv_aviso['long']
        inilat = adv_aviso['latg']
        lons2d = np.int16(np.ones((lons.shape[0],)+inilon.shape)*32767)
        lats2d = np.int16(np.ones((lons.shape[0],)+inilon.shape)*32767)
        # Find indices where inilon and inilat are not nans
        # (transpose to have F order: first index changing fastest, and the last index changing slowest)
        [ilon,ilat]= np.nonzero(~(np.isnan(inilon) & np.isnan(inilat)).T)
        # Assign values
        lons2d[:,ilat,ilon] = lons
        lats2d[:,ilat,ilon] = lats
        # Subsample time (every 4 time steps = 1 day!!!)
        _trajlon = np.ma.array(lons2d, mask=(lons2d==32767))
        _trajlat = np.ma.array(lats2d, mask=(lats2d==32767))
        # Sort by ascending time
        # Remove last point (initial position already in the fwd array)
        ibkw = np.argsort(_time)
        _time = _time[ibkw][:-1]
        _trajlon = _trajlon[ibkw,:,:][:-1,:,:]
        _trajlat = _trajlat[ibkw,:,:][:-1,:,:]
        #-----------------------

        #-----------------------
        # Combine the two
        print('\tCombine')
        time = np.concatenate((_time,time))
        trajlon = np.ma.concatenate((_trajlon,trajlon),axis=0)
        trajlat = np.ma.concatenate((_trajlat,trajlat),axis=0)
        # Need to unscale the variables before saving
        # (passign scale_factor and add_offset to to_netcdf automatically scales the variables to be saved)
        trajlon = trajlon*adv_aviso['lon_scale_factor']+adv_aviso['lon_offset']
        trajlat = trajlat*adv_aviso['lat_scale_factor']+adv_aviso['lat_offset']
        #-----------------------

        #-----------------------
        # Create xarray and save to netcdf
        # Find columns and rows (if any) with not all nan values
        inanlon = np.where(np.sum(~np.isnan(inilon),0) > 0)[0] 
        inanlat = np.where(np.sum(~np.isnan(inilat),1) > 0)[0]
        ds = xr.Dataset({'trajlon' : (['time','lat','lon'],trajlon[:,inanlat,:][:,:,inanlon]),
            'trajlat' : (['time','lat','lon'],trajlat[:,inanlat,:][:,:,inanlon])},
                         coords = {'lon':np.nanmean(inilon[np.ix_(inanlat,inanlon)],0),
                                   'lat':np.nanmean(inilat[np.ix_(inanlat,inanlon)],1),
                                   'time':time})
        # Save to netcdf
        # Create proper directory
        pathout = os.path.join(pathin,'%04d' % year)
        if not os.path.isdir(pathout):
            os.makedirs(pathout)
        ds.to_netcdf(os.path.join(pathout,'%s_oltraj_uv_global.nc' % (day)),
                encoding = {'trajlon' : {'scale_factor':adv_aviso['lon_scale_factor'],'add_offset':adv_aviso['lon_offset'],'dtype' : np.int16,'_FillValue' : 32767},
                    'trajlat' : {'scale_factor':adv_aviso['lat_scale_factor'],'add_offset':adv_aviso['lat_offset'], 'dtype' : np.int16,'_FillValue' : 32767},
                    'time' : {'dtype': np.double, 'units': f"days since {year:04d}-{month:02d}-{day[-2:]} 00:00:00", 'calendar': "proleptic_gregorian"}})
        #============================
        # Need to remove original mat files once netcdf is created
        #============================
        os.remove(os.path.join(pathin,fwdfile))
        os.remove(os.path.join(pathin,bkwfile))
        #-----------------------


def main():
   parser = argparse.ArgumentParser()
   parser.add_argument('--year', required=True, help="Year to convert")
   parser.add_argument('--month', required=True, help="Year to convert")
   args = parser.parse_args()
   year = int(args.year)
   month = int(args.month)
   combine_traj(year,month)


if __name__ == '__main__':
      main()
