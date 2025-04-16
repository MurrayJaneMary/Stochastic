function [S S_no_ws Sprime_no_ws] = VGPScatter_ws_func(VGPScatterInput);
%Calculates VGPs and then the angular dist of each from the mean
%Input  = [VGPScatterInput] = [VGPangdis kp nd] from MakeVGPScatterInput(rawdata); %

VGPangdis=[];
angdisin=[];
polecalcin=[];
%A1=A
%fisher % obtains MeanDec, MeanInc, N, k, a95, Pallat
%A2=A



N = size(VGPScatterInput,1);

scatt = zeros(9999,1); scatt=[];
for i = 1:N 
   
  		scatt(i,1) = ((VGPScatterInput(i,1)).^2 - (81/sqrt(VGPScatterInput(i,2))).^2 ./ VGPScatterInput(i,3));
      
      
end

sumscatt =sum(scatt);
sumscatt_no_ws = sum((VGPScatterInput(:,1)).^2);

S_no_ws = sqrt((1/(N-1))*sumscatt_no_ws);


if sumscatt > 0
   S = sqrt((1/(N-1))*sumscatt);
   
else
   S = - sqrt((1/(N-1)) * -sumscatt);
   
end
