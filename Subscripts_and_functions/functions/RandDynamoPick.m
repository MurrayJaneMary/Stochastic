function [timestep, slat, slong, dec, inc, F, plat, plong, VDM, VADM]...
    = RandDynamoPick(gh_full, Npicks, StartTime, EndTime, SiteLat, SiteLong)
%Function to make Npicks of randomly times field measurements in time range
%given and at specfied location


%rand('state',sum(100*clock))

max_degree = 10;%6
SizeTimeRange = EndTime-StartTime +1;

timestep  = ceil(rand(Npicks,1).*SizeTimeRange) + StartTime-1;
slat(1:Npicks, 1) = SiteLat; slong(1:Npicks, 1) = SiteLong; 
dec=zeros(Npicks,1);dec=[]; inc=zeros(Npicks,1); inc=[]; 
F=zeros(Npicks,1); F=[]; plat=zeros(Npicks,1); plat=[];
plong=zeros(Npicks,1); plong=[]; VDM=zeros(Npicks,1); VDM=[];
VADM=zeros(Npicks,1); VADM=[];

for ii = 1:Npicks

    
     
    [timestep(ii,1), slat(ii,1), slong(ii,1), dec(ii,1), inc(ii,1), ...
        F(ii,1), plat(ii,1), plong(ii,1), VDM(ii,1), VADM(ii,1)]...
                = QueryModel(gh_full,timestep(ii),SiteLat,SiteLong,max_degree);

end
