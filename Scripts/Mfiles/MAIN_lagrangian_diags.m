%!/usr/bin/octave
%============================================================
% This script test particle advection for AVISO
% velocity fields
%
% F. Nencioli May 16 2019
%============================================================
function MAIN_lagrangian_diags(year,month)

   % Defines wether the velocity field is kept constant (1) or not (0)
   freeze=0;
   % Used to define numdays!!!!
   advdays = 30;
   % Used to define particle deployment resolution (hence howm ofter in space trajectories are computed)
   % in degrees (lat or lon)
   delta0 = 0.125; % 1/8 of degree

   pkg unload netcdf
   pkg load netcdf

   %=============================
   % Path to AVISO data and FSLE main directory
   global DATADIR=['../../Data/Daily/'];
   global LAMTADIR=['./Lamta.dev/'];
   global outpath = '../../Processed/Lagrangian_traj/';
   %=============================

   addpath(LAMTADIR)

   % String with month and year
   strym = sprintf('%04d%02d',year,month);
   % Cycle through all days of given month
   for iday = 1:eomday(year,month)
   % Comment above and uncomment below to process only first day of month for test!!!
   %for iday = 1:1
      % Iniday is a string in the format 'yyyymmdd'
      % pass directly dayprod
      iniday = sprintf('%s%02d',strym,iday);
      dayprod=datenum(iniday,'yyyymmdd');
      dayv=datevec(dayprod);

      %--------------------------
      % Load the velocity field
      day0=datenum(dayv)-datenum([1950 1 1]);

      % Advect for both forward and backward times
      for iadvect = -1:2:1
         numdays = iadvect*advdays
         dayf=day0+numdays;
         % inilon and inilat to be used for particle deployment!!!
         tic()
         [lonv,latv,uvmask,err] = cmems_load(min([day0,dayf]),max([day0,dayf]));
         disp('Loading')
         toc()
         disp(' ')
         % %--------------------------
         % % Remove rows at the north pole
         % disp('Remove N pole particles')
         % inilon=inilon(1:end-1,:);
         % inilat=inilat(1:end-1,:);
         % %--------------------------

         [lonf,latf] = cmems_advect(dayv,lonv,latv,delta0,uvmask,numdays,freeze);
      end % for iadvect = -1:2:1
   end % for iday = 1:31
   clear outpath
