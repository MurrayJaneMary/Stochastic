function VGPScatterInput = MakeVGPScatterInput_NEW(AA);
%Makes VGPScatterInput = [VGPangdis, kp, ndir]
plat=AA(:,1); plong=AA(:,2);, Maglat=AA(:,3);, kdir=AA(:,4); ndir=AA(:,5);

N = size(plat,1);
FisherPoleOut = fisherfun([plong, plat]);
Plongm_temp = FisherPoleOut(1,1);Platm_temp = FisherPoleOut(1,2); 

for ii = 1:N
    
    angdisin(ii,:) =[plat(ii), plong(ii), Platm_temp, Plongm_temp]; %Creates matrix for calculating angular distances
    
end
VGPangdis = angdis(angdisin);

kp = kdir ./ (1/8 * (5+18*(sin(radians(Maglat))).^2 + 9*(sin(radians(Maglat))).^4));

VGPScatterInput = [VGPangdis kp ndir];