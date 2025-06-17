function [y_median Nx Nx_Sufficient] = RunningMedianThreshold(x, y, threshold)

%Function to return median of a variable y as a function of increasing x
% where x and y are both vectors of same length
%Set threshold as a fraction (e.g. 0.05) whereby Nx_Sufficient is enough
%steps for y_median to be within this range of the global y_median

y_median_global = median(y,2);
N=size(x,1);
% y_median  = zeros(N);y_median  =[];
% Nx = zeros(N);Nx  =[];
% dev  = zeros(N);dev  =[];
% Nx_Sufficient = zeros(N);Nx_Sufficient  =[];



for i = 1:N

    y_median(:,i) = median(y(:,1:i),2);

    for j = 1:size(y,1)


        Nx(j,i) = i;
        dev(j,i) = abs((y_median(j,i) - y_median_global(j))./ y_median_global(j));
        if dev(j,i) > threshold
            Nx_Sufficient(j,i) = 0;
        else
            Nx_Sufficient(j,i) = 1;
        end

    end
end