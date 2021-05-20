function [glon,glat,Uw,Vw,lonw,latw,U,V]=getncUV_cmems(nday)

% Get date string
ncnes=datenum([1950 1 1]);
dv=datevec(nday+ncnes);
day=sprintf('%d%02d%02d',dv(1),dv(2),dv(3));

% Define file to be used
global DATADIR
fname = sprintf('%s.nc',day);
names2glob = sprintf('%s/%04d/%s',DATADIR,dv(1),fname);
fnames = glob(names2glob);

if (~isempty(fnames))
   if length(fnames)!=1
      print(['Error ',str(length(fnames)),' files found'])
      return;
   else
      fnameout=fnames{1};
   end
else
   disp('File not found.');
   Uw=[];
   Vw=[];
   lonw=[];
   latw=[];
   return;
end
ncfile=sprintf('%s',fnameout);

% Define u and v variable string
uvar = 'ugos';
vvar = 'vgos';

ncfile
NbLongitudes=ncread(ncfile,'longitude');
%=========================================
% Lamta code assumes lon to be -180 180
% Need to convert from 0 360 to -180 180
if any(NbLongitudes>180)
    NbLongitudes(NbLongitudes>180) = NbLongitudes(NbLongitudes>180)-360;
end
%=========================================
NbLatitudes=ncread(ncfile,'latitude');
U=ncread(ncfile,uvar)(:,:);
V=ncread(ncfile,vvar)(:,:);
%================================
% IMPORTANT!!!!
% Velocities on land need to be  set to 0!!!!
uvmask = isnan(U) & isnan(V);
U(uvmask)=0;
V(uvmask)=0;
%================================
% Conversion m/s to cm/s
U=U.*100;
V=V.*100;
% Create lon lat matrices
nlon=length(NbLongitudes);
nlat=length(NbLatitudes);

lon=reshape(NbLongitudes,nlon,1);
lat=reshape(NbLatitudes,nlat,1);
% Create mask for initial deployment
% (also used to derive deg velocities!!!)
[glon,glat]=meshgrid(lon,lat);

% Velocities in deg/s (Used in Lyap analysis)
RT=6371e5;
Ug=U./(RT.*cos(glat'./180*pi)).*180/pi;
Vg=V*180/pi./RT;

% Deployments on land points are masked
glon(uvmask')=NaN;
glat(uvmask')=NaN;

lonw=lon;
latw=lat;
Uw=Ug;
Vw=Vg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In case of global domain!!!
% Need to add one last column to make it completely cyclical!!!
[lonw,lonsrti]=sort(lonw);
Uw=Uw(lonsrti,:);
Vw=Vw(lonsrti,:);
Uw=Uw([1:nlon 1],:);
Vw=Vw([1:nlon 1],:);
lonw=lonw([1:nlon 1]);
lonw(end)=lonw(end)+360;

U=U(lonsrti,:);
V=V(lonsrti,:);
U=U([1:nlon 1],:);
V=V([1:nlon 1],:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
