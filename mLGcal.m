function lbdG=mLGcal(Tutc) 
mjd = mJulday(Tutc);
mday = mjd - mJulday([2000;1;1;12;0;0]);
lbdG = mod(280.4606184 + 360.9856122863*mday,360); 
lbdG = lbdG*pi/180;


function MJD = mJulday(Tutc)
Ye=Tutc(1);Mo=Tutc(2);Da=Tutc(3);Hr=Tutc(4);Mi=Tutc(5);Se=Tutc(6);

J = Da - 32075 + fix(1461.*(Ye + 4800 + fix((Mo-14)/12))./4)...
       +fix(367.*(Mo - 2 - fix((Mo -14)/12)*12)./12)...
       -fix(3.*((Ye + 4900 + fix((Mo-14)/12))./100)./4);
JD = J-0.5+Hr./24.0+Mi./1440.0+Se./86400.0;
MJD= JD-2400000.5;