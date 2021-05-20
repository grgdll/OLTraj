function [lons,lats,M]=trj2pos(trj,Nstep,numpts)
%function M=trj2pos(trj,Nstep,numpts)
%M(times,1=lon 2=lat,pt) e.g., all longitudes of particle 5: lon=M(:,1,5);
% all initial longitudes: M(1,1,:);
%all final latitudes: M(end,2,:);
M=reshape(trj,[Nstep,2,numpts]);

lons=reshape(M(:,1,:),Nstep,numpts);
lats=reshape(M(:,2,:),Nstep,numpts);
	
