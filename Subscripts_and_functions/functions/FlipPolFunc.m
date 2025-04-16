function [newVGPlat, newVGPlong, newdec, newinc] = ...
    FlipPolFunc(VGPlat, VGPlong, olddec, oldinc)

%Function to flip polarity of VGPs that are southern hemisphere and
%associated directions as well if they are input

DirectionsExist = exist('olddec');
N = size(VGPlat,1);
ii = 1;
while ii <= N
    
    if VGPlat(ii,1) < 0
        
        %Needs flipping
        newVGPlong(ii,1) = VGPlong(ii,1) + 180;
                
        if newVGPlong(ii,1) >= 360
            newVGPlong(ii,1) = newVGPlong(ii,1) -360;
        end
        newVGPlat(ii,1) = VGPlat(ii,1) .* -1;
        
        if DirectionsExist == 1
            %direction too needs flipping
            newdec(ii,1) = olddec(ii,1) + 180;
            if newdec(ii,1) >= 360
                newdec(ii,1) = newdec(ii,1) -360;
            end
            newinc(ii,1) = oldinc(ii,1) .* -1;
        end
        
    else
        newVGPlat(ii,1) = VGPlat(ii,1);
        newVGPlong(ii,1) = VGPlong(ii,1);
    
        if DirectionsExist == 1
            newdec(ii,1) = olddec(ii,1);
            newinc(ii,1) = oldinc(ii,1);
        end
        
    end
    
    ii=ii+1;
end