function  ang_sum  = angdissum( x, vgplat, vgplon )
%angdissum Returns sum of angular distances between test pole and vpgs
%   x - test pole, x(1) = lat, x(2) = lon
    test_lats = ones(size(vgplat))*x(1);
    test_lons = ones(size(vgplat))*x(2);
    ang_dist = angdis([test_lats,test_lons,vgplat,vgplon]);
    ang_sum = sum(ang_dist);

end

