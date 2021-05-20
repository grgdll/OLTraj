function [glon,glat,err]=cmems_load(day0,dayf)
% Load cmems velocity field from AVISO
% returns glon and glat used to initialize particle deployment

err=0;

% In case day0 and dayf are passed as date vectors
if(length(day0)>1)
   day0=datenum(day0)-datenum([1950 1 1]);
end
if(length(dayf)>1)
   dayf=datenum(dayf)-datenum([1950 1 1]);
end

% Test day to be loaded and retrieve values for few paramters
onedayinupd=datenum([2019 4 1])-datenum([1950 1 1]);

disp(' ')
disp('%--- Define the type of velocity field ---')
% Function in lamta.dev/lamta_all.cc which call select_sys from lamta.dev/field.h
% It defines the type of velocity field (4 => user defined Lut_field)
select(4); 

% Load test day
[glon,glat,Uw,Vw,lonw,latw]=getncUV_cmems(onedayinupd);

if(isempty(Uw))
	disp(sprintf('Error in finding velocity files for reference day %d (upd).',onedayinupd));
	err=2;
	return;
end

sz=size(Uw);

numdays=(dayf-day0)+1;

disp(' ')
disp('%--- Set parameters of LUT field --- ')
% set_par defined in lamta.dev/field.h for each field
set_par([lonw(1) lonw(end) sz(1) latw(1) latw(end) sz(2) day0*60*60*24 dayf*60*60*24 numdays]');
% No need to print_par: already in set_par !!!
%print_par();
disp(' ')

disp('%--- Set geometry of LUT field --- ')
% field_geometry(disttype,datatype,gridtype) defined in lamta.dev/lamta_all.cc
%
% disttype=1 (Euclidean), 2 (sphere).
% datatype=1 (deg./sec.), 2 (cm/sec.)
% gridtype=0 (flat), 1 (sphere, regular), 2 (sphere Mercator)
field_geometry(2,1,1); % new grid for AVISO netcdf files (Cartesian instead of Mercator) 
% No need to print_par: field geometry is not displayed !!!
% print_par();
disp(' ')

disp('%--- Load velocities ---')
ctr=0;
% Cycle through each day
for ctday=day0:1:dayf
	% Read global velocity field 
	[glon,glat,Uw,Vw,lonw,latw]=getncUV_cmems(ctday);

	if(isempty(Uw))
		disp(sprintf('Error in finding nrt velocity files (day %d)',ctday));
		err=1;
		return;
	end
	% Update the velocity field in the LUT field used
	% for the Lagrangian analysis
	disp(['% ',datestr(ctday+datenum([1950 1 1]))])
	LUT_frame_fill(Uw,Vw,ctr);
	ctr=ctr+1;
end
