import numpy as np
import os
import datetime as dt
import imp
import scipy.io as sio
import gsw
import ipdb
#import AVISO_utils


def set_path():
   pathin=os.path.join('..','..','..','Processed','Lagrangian_traj')
   outdir=os.path.join('..','..','..','Figures','Lagrangian_traj')
   return pathin, outdir


def read_traj(infile):
   '''Read mat file with trajectory defined by infile
      (The file is in ../../../../Argo/SSH/Pyfiles/AVISO_utils.py
      returns a dictionary with various variables to be read
   '''
   pathin, outdir = set_path()
   #AVISO_utils=imp.load_source('AVISO_utils','./AVISO_utils.py')
   adv_aviso = sio.loadmat(os.path.join(pathin,infile),squeeze_me=True)
   return adv_aviso


def get_dates(day0,dayf,Nstep):
   ''' Returns an list with the dates corresponding to the various timesteps 
       (computation based on passages in ../Mfiles/MAIN_advect_aviso.m)
   '''
   tcnes = dt.datetime(1950,1,1,0,0)
   initime = tcnes + dt.timedelta(day0)
   dtime = (dayf-day0)/Nstep
   dates = [initime+dt.timedelta(dtime*i) for i in np.arange(Nstep)]
   return dates


def compute_area(lon,lat):
   # Lon and lat are 2D matrices !!!
   # Distance is computed in meters !!!
   dy = [gsw.distance(lon[:,i],lat[:,i]).squeeze() for i in np.arange(lon.shape[1])]
   dy = np.asarray(dy).T
   dx = [gsw.distance(lon[i,:],lat[i,:]).squeeze() for i in np.arange(lon.shape[0])]
   dx = np.asarray(dx)
   # Area is average between areas computed using lon distance from 
   # top or bottom lat as reference
   areatop = dx[1:,:]*dy[:,:-1]
   areabottom = dx[:-1,:]*dy[:,:-1]
   area = (areatop+areabottom)/2.
   return area
#=============================
#    # Use area of sphere element r**2sin(pi-lat) dlat dlon
''' Very little difference (dlon, dlat are small enough that cartesian area
    for rectangle works as well):
    Even with this method is important where to set the lat from where the 
    area is computed (either from bottom or from top)
'''
#    # Convert origin latitude from equator to north pole
#    lat = 90-lat
#    # convert angles in radians
#    lat = np.deg2rad(lat)
#    lon = np.deg2rad(lon)
#    dlat = np.diff(lat,axis=0)
#    dlon = np.diff(lon,axis=1)
#    area = gsw.earth_radius**2 * np.sin(lat[:-1,:-1]) * -dlat[:,:-1] * dlon[:-1,:]
#    area2 = gsw.earth_radius**2 * np.sin(lat[1:,:-1]) * -dlat[:,:-1] * dlon[:-1,:]
#=============================

###################
# To get area of single element simply divide area of cell by number of particles
# (the one below is more to check everything is ok with the computation)
##################


def get_lonlat_vertices(lon,lat,delta0):
   # Computes the position of the vertices of the volumes associated with
   # deployed particle
   # (returns 2D matrices that can be used with compute_area) 
   #
   # Longitude
   # 1) Add extra row
   lon = np.vstack((lon[0,:],lon))
   # 2) shift by delta0/2
   lon += delta0/2.
   # 3) add first column shifted back by delta0
   icol = lon[:,0].reshape(-1,1)-delta0
   lon = np.hstack((icol,lon))
   # Latitude
   # 1) add extra column
   lat = np.hstack((lat[:,0].reshape(-1,1),lat))
   # 2) shift by delta0/2
   lat += delta0/2.
   # 3) add first row shifted back by delta0
   irow = lat[0,:]-delta0
   lat = np.vstack((irow,lat))
   return lon, lat


def find_cell(lonpart,latpart,longrid,latgrid):
   # - lonpart and latpart are 1d vectors with positions of all particles
   # - longrid and latgrid are 1d vectors with position of cell vertices
   # - returns 2 1D vector with i and j indexes of the cell within 
   #   which each particle belongs

   # Initialize particle cell belonging to -1
   jcellpart = -np.ones(lonpart.shape)
   icellpart = -np.ones(latpart.shape)
   for i in np.arange(longrid.shape[0]-1):
      jcellpart[np.where((lonpart>=longrid[i]) & (lonpart<longrid[i+1]))[0]]=i
   for j in np.arange(latgrid.shape[0]-1):
      icellpart[np.where((latpart>=latgrid[j]) & (latpart<latgrid[j+1]))[0]]=j
   return icellpart, jcellpart
