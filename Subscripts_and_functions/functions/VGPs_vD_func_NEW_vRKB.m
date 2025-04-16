function [slatm, slongm, Decm, Incm, Maglatm, platm, plongm, Norig, S_nocut, Nnew, cutoff, low95, TrueRes, high95]...
    =VGPs_vD_func_NEW(slat, slong, dec,inc, ndir, kdir, Nboot)
%Calculates VGPs and then scatter and 95% confidence limits after applying Vandamme (1994) variable cutoff

%clearvars -except slat slong dec inc ndir kdir Nboot

rand('state',sum(100*clock))

%Find mean direction

% Make VGPs
[plat, plong] = polecalc([dec, inc, slat, slong]);


Norig = size(plat, 1);

Maglat(1:Norig,1) = 90 - angdis([slat, slong, plat, plong]);
% MedianMaglat = median(Maglat);
if any(isnan(Maglat))
%     keyboard;
    Maglat = Maglat(~isnan(Maglat));
    plat = plat(~isnan(Maglat));
    plong = plong(~isnan(Maglat));
    kdir = kdir(~isnan(Maglat));
    ndir = ndir(~isnan(Maglat));
    dec = dec(~isnan(Maglat));
    inc = inc(~isnan(Maglat));
    slat = slat(~isnan(Maglat));
    slong = slong(~isnan(Maglat));
end
MedianMaglat = median(Maglat);

%Calculate S_nocut
% Make VGPScatterInput  = [VGPangdis kp ndir]
AA = [plat, plong, Maglat, kdir, ndir];
% for some reason, if plat is opposite sign, generates NaNs to remove those rows
platSign = mode(sign(AA(:,1)));
slat = slat(sign(AA(:,1))==platSign);
slong = slong(sign(AA(:,1))==platSign);
dec = dec(sign(AA(:,1))==platSign);
inc = inc(sign(AA(:,1))==platSign);
AA = AA(sign(AA(:,1))==platSign,:);

[VGPScatterInput] = MakeVGPScatterInput_NEW(AA);%msgbox('done'); pause;
if any(any(isnan(VGPScatterInput)))
    keyboard;
end
[S_nocut ] = VGPScatter_ws_func([VGPScatterInput]);

%***********variablecutoff.m
S= S_nocut;
stop = 0;
p=1;
try
    combo = [VGPScatterInput, AA, dec, inc, slat, slong];
catch
    keyboard;
end

if any(any(isnan(combo)))
    keyboard;
end
% tic;
while stop == 0    
%     tic;
    if S < 0
        cutoff = 5;
    else 
        cutoff = 1.8 * S + 5;
    end
    
    [VGPScatterInput] = MakeVGPScatterInput_NEW(AA);
    N=size(AA,1);
    Maxangdis = max(VGPScatterInput(:,1));

    %Make a combo array with [VGPangdis, kp, ndir, plat, plong, Maglat,
    %kdir, ndir, dec, inc, slat, slong]
    if Maxangdis > cutoff
        %gr=9999
%         VGPScatterInput =[];
%         AA=[];
%         rejected = combo(combo(:,1)>cutoff,:);
        combo(combo(:,1)>cutoff,:) = [];
%         while Maxangdis > cutoff
%             sorted=sortrows(combo,1); %sort by angular distance so that can take all below the maximum
% %             combo = [];
%             combo = sorted(1:N-1,:);
%             
% %             rejected(p,:) = sorted(N,:); %records the rejected site
%             
%             N = N-1; p = p+1;
%             %return
%             Maxangdis = max(combo(:,1));
%         end

        AA = combo(:, 4:end);
        [VGPScatterInput]= MakeVGPScatterInput_NEW(AA); %[VGPangdis kp nd]
        [S] = VGPScatter_ws_func(VGPScatterInput);
        combo = [VGPScatterInput, AA];
        if any(isnan(S))
            keyboard;
        end
    else
        stop = 1;
    end
% toc;    
end
% toc;
%end variable cutoff

% Maxangdis = max(combo(:,1));
TrueVGPScatterInput = VGPScatterInput;
Nnew = size(TrueVGPScatterInput,1);
TrueRes = S;
% if any(isnan(TrueRes))
%     keyboard;
% end
%TrueRes_no_ws = S_no_ws;

%Calculate mean direction, pole, and site location from surviving data
FisherDirOut = fisherfun([combo(:,9),  combo(:,10)]);
Decm = FisherDirOut(1,1); Incm = FisherDirOut(1,2);

FisherPoleOut = fisherfun([combo(:,5),  combo(:,4)]);
platm = FisherPoleOut(1,2); plongm = FisherPoleOut(1,1);

FisherSiteOut = fisherfun([combo(:,12),  combo(:,11)]);
slatm = FisherSiteOut(1,2); slongm = FisherPoleOut(1,1);

%Provide mean magnetic latitude
Maglatm = 90 - angdis([slatm, slongm, platm, plongm]);

%Run bootstraps
if Nboot == 0; low95 = 0; high95 = 0;
else
    realCheck = 0; % check to see if confidence values are real, not sure why needed
    while realCheck<2
        for jj=1:Nboot %determines the number of random picks for the bootstrap
            VGPScatterInput(:,:,jj) = TrueVGPScatterInput(...
                        randsample(1:size(TrueVGPScatterInput,1),size(TrueVGPScatterInput,1),'true'),:);
        end        
    %     randidx = randi(1,size(TrueVGPScatterInput,1),[size(TrueVGPScatterInput,1),Nboot]);
    %     VGPScatterInput = TrueVGPScatterInput(
        S = VGPScatter_ws_func(VGPScatterInput);
    %         randscatt(jj,1)=S;
        randscatt = S;

        conf95 = prctile(randscatt,[2.5,97.5]);
        low95 = conf95(1);
        high95 = conf95(2);
        
        realCheck = isreal(low95) + isreal(high95);
        if realCheck==2
            break;
        else
            disp('imag conf');
        end
    end
end


%[slatm, slongm, Decm, Incm, Maglatm, Platm, Plongm, Norig, Nnew, cutoff, low95, TrueRes, high95]
%headings=['Dec ' 'Inc ' 'No ' 'N ' 'Cutoff ' 'Dmax ' 'R ' 'k ' 'a95 ' 'Pallat ' 'Sl ' 'S ' 'Su ']
%Results = [Decm Incm Norig Nnew cutoff Maglatm platm plongm low95 TrueRes high95]