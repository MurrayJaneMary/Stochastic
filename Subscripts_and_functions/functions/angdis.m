% Calculates ang dist between pole positions. Input A=[lat1 long1 lat2 long2]

function b = angdis(A)
%A=rawdata;

lat1r = radians(A(:,1));
long1r = radians(A(:,2));
lat2r = radians(A(:,3));
long2r = radians(A(:,4));

% angular distance =   arccos(Sin(lat1) * Sin(lat2) + Cos(lat1) * Cos(lat2) * Cos(lon1-lon2)) 
b = degrees(acos(sin(lat1r) .* sin(lat2r) + cos(lat1r) .* cos(lat2r) .* cos(long1r-long2r)));

