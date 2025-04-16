%Calculates fisher stats from a series of Decs and Incs
%Input A = [Dec Inc] in degrees
%Output B = [MeanDec MeanInc N k a95 Pallat];

function  B = fisherfun(A)

Decr=radians(A(:,1));
Incr=radians(A(:,2));

l = cos(Decr) .* cos(Incr);
m = sin(Decr) .* cos(Incr);
n = sin(Incr);

R = sqrt((sum(l))^2+(sum(m))^2+(sum(n))^2);

X=1/R * sum(l);
Y=1/R * sum(m);
Z=1/R * sum(n);

MeanDec=degrees(atan(Y/X));

if X < 0
    MeanDec=MeanDec+180;
else
    if Y < 0
        MeanDec=MeanDec+360;
    end
end

MeanInc=degrees(asin(Z));
N=size(A,1);

if N > 1
    k=(N-1)/(N-R);
    a95 = degrees(acos(1-((N-R)/R)*((1/0.05)^(1/(N-1))-1)));
    %a95_2 = 140 / sqrt(k * N);
else
    k = 0;
    a95 = 0;

end
Pallat= degrees(atan((tan(radians(MeanInc)))/2));

B = [MeanDec MeanInc N k a95 Pallat];



