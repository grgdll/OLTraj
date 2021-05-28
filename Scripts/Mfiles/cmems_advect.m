function [lons,lats] = cmems_advect(dayv,lonv,latv,delta0,uvmask,numdays,freeze)
tic()
global outpath

% Compute longitude and latitude advection
day0=datenum(dayv)-datenum([1950 1 1]);

dayf=day0+numdays;
%----------------------------
lonvv=lonv(1):delta0:lonv(end);
latvv=latv(1):delta0:latv(end);
% New grid (deployment)
[long,latg]=meshgrid(lonvv,latvv);
% Old grid (velocities)
[lonmask,latmask]=meshgrid(lonv,latv);
% New velocities used for masking
uvmaskg = interp2(lonmask,latmask,double(uvmask'),long,latg);
% Need to convert long from 0 360 to -180 180
if any(any(long>180))
    long(long>180) = long(long>180)-360;
end
% Mask points on land (velocity = 0)
long(uvmaskg>=1)=NaN;
latg(uvmaskg>=1)=NaN;

% Array of initial deployment points
pts=zeros(length(long(:))*2,1);
pts(1:2:end)=long(:);
pts(2:2:end)=latg(:);
%=============================================
% Remove points on land (~30% of the points)
% This saves ~10 sec in adv time and ~10 sec in save time
% 70 sec to load 30 days, 30 sec to advect all points and 30 sec to save all points
% Removing land points is 20/(70+30+30) = ~15% less time to compute one day
pts(isnan(pts))=[];
numpts=length(pts(1:2:end));
%=============================================

daysteps = 4;
tspan=([day0 dayf]'*60*60*24);
Nstep=abs(diff(tspan))/(60*60*24)*daysteps; % 4 t-steps per day

if freeze==1
   freezetime(tspan(1));
end

trj=RK4(tspan,pts,Nstep);
[lons,lats]=trj2pos(trj,Nstep,numpts);
latf=trj(2*Nstep:Nstep*2:end);
lonf=trj(Nstep:Nstep*2:end);

daystr=datestr(dayv,'yyyymmdd');

if numdays>0
   advtag = 'fwd';
else
   advtag = 'bkw';
end

disp('Advection')
toc()
disp(' ')

% Convert output variables:
% 1) Select only daily output
lons = lons(1:daysteps:end,:);
lats = lats(1:daysteps:end,:);
% 2) convert from double to 16 bit integers
n=16;
% Longitude scale_factor and offset
lonmin = -180;
lonmax = 181; % This to leave 32767 for the _FillValue
lon_scale_factor = (lonmax-lonmin) / (2 ** n - 1);
lon_offset = lonmin + 2 ** (n - 1) * lon_scale_factor; 
% Latitude scale_factor and offset
latmin = -90;
latmax = 91; % This to leave 32767 for the _FillValue
lat_scale_factor = (latmax-latmin) / (2 ** n - 1);
lat_offset = latmin + 2 ** (n - 1) * lat_scale_factor; 
% Convert from double to int16
lons = int16(round((lons - lon_offset) / lon_scale_factor));
lats = int16(round((lats - lat_offset) / lat_scale_factor));

tic()
if freeze==1
   % Save only the variables needed to create the netcdf file with the traj
   save('-mat-binary',[outpath,'/traj_aviso_',daystr,'_',advtag,'_frozen.mat'],'lons','lats','long','latg','day0','dayf','Nstep','lon_scale_factor','lon_offset','lat_scale_factor','lat_offset') % to save the file in matlab format
else
   % Save only the variables needed to create the netcdf file with the traj
   save('-mat-binary',[outpath,'/traj_aviso_',daystr,'_',advtag,'.mat'],'lons','lats','long','latg','day0','dayf','Nstep','lon_scale_factor','lon_offset','lat_scale_factor','lat_offset') % to save the file in matlab format
end
disp('Saving')
toc()
