#this is needed because the current code generate files with longitudes that are not sorted.
#if you run again the code, you will need to modify it just before the files are written

#run in parallel using:
#find /data/datasets/Projects/OLTraj/From_AVISO/OLTraj/Processed/Lagrangian_traj/????/ -name "[1-2]*nc" | xargs -n1 -P1 sh -c 'sleep 2; echo $0' | xargs -n1 -P4 echo python3 shift_lon.py $1
import xarray as xr
import sys
import os

fn = sys.argv[1]
print(fn[-28:])

ds = xr.open_dataset(fn)
ds = ds.sortby("lon")

fntmp = fn[:-28]+'tmp_'+fn[-28:]
#print(fntmp)
ds.to_netcdf(fntmp)

ds.close()

fnnew = fn 

os.rename(fntmp, fnnew)


# after correcting all files using the above script, I got this suggestion based on the NCO scripts here: https://sourceforge.net/p/nco/discussion/9829/thread/d023efcec1/?limit=25
# ncks --msa_usr_rdr -d lon,-179.875,-0.125 -d lon,0.0,180.0 in.nc out.nc
# it runs twice faster than this script!
