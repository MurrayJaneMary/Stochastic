function    B=g2fxFunc_vRKB(gh,lmax,nlat,elong,alt)
%****************************************************
%  Direct conversion of igrf11_f to MATLAB

%   OUTPUT
%     x     = north component (nT) if isv = 0, nT/year if isv = 1
%     y     = east component (nT) if isv = 0, nT/year if isv = 1
%     z     = vertical component (nT) if isv = 0, nT/year if isv = 1
%     f     = total intensity (nT) if isv = 0, rubbish if isv = 1
%     B     = [x; y; z]
%

%************************************************************************************
%global gh
%%%%%%%%%%%Variable Changes:
isv=0;
itype=2;
colat=90-nlat;
%****************************************************
cl=NaN*ones(1,13); 
sl=NaN*ones(1,13);  
x     = 0;
y     = 0;
z     = 0;
   
nmx   = lmax;%changed from 10 in original
% fprintf('lmax: %d\n',lmax)

% nmx = 9;
nc    = nmx*(nmx+2);

kmx   = (nmx+1)*(nmx+2)/2;

r     = alt+6371.2;
one   = colat*0.017453292;
ct    = cos(one);
st    = sin(one);
one   = elong*0.017453292;
cl(1) = cos(one);
sl(1) = sin(one);
cd    = 1.0;
sd    = 0.0;
l     = 1;
m     = 1;
n     = 0;

if itype==1  %Skip if itype==2
    %   Conversion from geodetic to geocentric coordinates (using the WGS84 spheroid)
    a2    = 40680631.6;
    b2    = 40408296.0;
    one   = a2*st*st;
    two   = b2*ct*ct;
    three = one + two;
    rho   = sqrt(three);
    r     = sqrt(alt*(alt + 2.0*rho) + (a2*one + b2*two)/three);
    cd    = (alt + rho)/r;
    sd    = (a2 - b2)/rho*ct*st/r;
    one   = ct;
    ct    = ct*cd -  st*sd;
    st    = st*cd + one*sd;
end
ratio = 6371.2/r;  
rr    = ratio*ratio;
p=NaN*ones(1,kmx);
q=NaN*ones(1,kmx);

%     computation of Schmidt quasi-normal coefficients p and x(=q)
p(1)  = 1.0;
p(3)  = st;
q(1)  = 0.0;
q(3)  = ct;
for k=2:kmx
    if n<m
        m     = 0;
        n     = n + 1;
        rr    = rr*ratio;
        fn    = n;
        gn    = n - 1;
    end
    fm    = m;
    if m==n && k~=3
        one   = sqrt(1.0 - 0.5/fm);
        j     = k - n - 1;
        p(k)  = one*st*p(j);
        q(k)  = one*(st*q(j) + ct*p(j));
        cl(m) = cl(m-1)*cl(1) - sl(m-1)*sl(1);
        sl(m) = sl(m-1)*cl(1) + cl(m-1)*sl(1);
    end
    if m~=n && k~=3
        gmm    = m*m;
        one   = sqrt(fn*fn - gmm);
        two   = sqrt(gn*gn - gmm)/one;
        three = (fn + gn)/one;
        i     = k - n;
        j     = i - n + 1;
        p(k)  = three*ct*p(i) - two*p(j);
        q(k)  = three*(ct*q(i) - st*p(i)) - two*q(j);
    end
    %Synthesis of x, y, and z
    lm    =  l;
%     fprintf('lmax: %d\tnmx: %d\tkmx: %d\tlm: %d\n',lmax, nmx, kmx, lm)
    one   = (gh(lm))*rr;
    %DEBUG 
%     [nyear,ncoef]=decode_coeff_pointer(lm);
%     fprintf('nyear %6i %5i %5i \n',lm,nyear,ncoef); 
    if m~=0
        two   = (gh(lm+1))*rr;
        three = one*cl(m) + two*sl(m);
        x     = x + three*q(k);
        z     = z - (fn + 1.0)*three*p(k);
        if st~=0
            y     = y + (one*sl(m) - two*cl(m))*fm*p(k)/st;
        else
            y     = y + (one*sl(m) - two*cl(m))*q(k)*ct;
        end
        l     = l + 2;
    else
        x     = x + one*q(k);
        z     = z - (fn + 1.0)*one*p(k);
        l     = l + 1;
    end
    m     = m + 1;
end %k loop

%     conversion to coordinate system specified by itype
one   = x;
x     = x*cd +   z*sd;
z     = z*cd - one*sd;
B     =[x; y; z];
%f     = sqrt(x*x + y*y + z*z);
% return