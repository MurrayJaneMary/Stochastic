function [ out_val, resid ] = ModelGfast_v2(lats, S)
%ModelGfast Uses fmincon to return a, b and resid
%   if optimization toolbox present, uses fmincon. 
%   otherwise uses fminsearch. 
%   expects lats, S, weights (optional) and initial guess [a,b] (optional)
    
    % check if optimization toolbox is present
    if license('test','optimization_toolbox')
        isOTP = 1;
    else
        isOTP = 0; 
    end
    
    % check if weights were passed, if not set to 1
    % currently not working (?)
    %if (~exist('weights','var'))
     %   weights = ones([size(lats,1),1]);
    %end
   
    % check if initial guess for a and b were provided, otherwise default
    % to 25 and 0.5
%     if (~exist('guess','var'))
%         guess = [25,0.5];
%     end
%     weights = ones(size(S)); 
    % define function for square diff between model G guess and data (S)
    % g(2,1) are a, b; D = data (S); L = lats; n = weights
%     modelG_func = @(g,D,L,n)sum(((D-sqrt(g(1)^2+(g(2)*L).^2)).^2)*n);
    modelG_func = @(g,D,L) sum(((D-sqrt(g(1)^2+(g(2).*L).^2)).^2));

    % display options
    options = optimset('Display','none'); % no display
%     options = optimset('Display','final'); % default
%     options = optimset('Display','iter'); % each iter
    
    % minization search
    if isOTP
        [out_val,fval] = fmincon(@(g) modelG_func(g,S,lats),[25,0.5],[],[],[],[],[1,0],[90,1],[],options);
    else
        [out_val,fval] = fminsearch(@(g) modelG_func(g,S,lats),[25,0.5]);
    end
    
    a = out_val(1);
    b = out_val(2);
    resid = fval;
end

