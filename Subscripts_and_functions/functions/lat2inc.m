% calculates inclination from latitude based on pmag equation

function  inc = lat2inc(lat)
l=radians(lat);
I= atan(2*tan(l));
inc= degrees(I);