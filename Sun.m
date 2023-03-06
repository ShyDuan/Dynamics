function Si=Sun(Epoc,time)
tJ2015=mJulday(Epoc)-mJulday([2015;1;1;12;0;0]);
   dtd = time/86400+tJ2015+5479.0;   
    dtc = dtd/36525.0;	
    Sun_e=0.016709;
    Sun_i=23.4393*pi/180;
    Sun_w=(282.9373+0.32* dtc)*pi/180;
    Sun_M=(357.5291+0.9856*dtd)*pi/180;
    if  Sun_M >= 2*pi
         Sun_M = mod(Sun_M, 2*pi);
    elseif Sun_M<0
       Sun_M = mod(Sun_M, 2*pi)+2*pi;
    end
    Sun_f = Sun_M+2*Sun_e*sin(Sun_M);
    Si = [cos(Sun_w+Sun_f);...
          sin(Sun_w+Sun_f)*cos(Sun_i);...
          sin(Sun_w+Sun_f)*sin(Sun_i)];    

function MJD = mJulday(Tutc)
Ye=Tutc(1);Mo=Tutc(2);Da=Tutc(3);Hr=Tutc(4);Mi=Tutc(5);Se=Tutc(6);

J = Da - 32075 + fix(1461.*(Ye + 4800 + fix((Mo-14)/12))./4)...
    + fix(367.*(Mo-2-fix((Mo-14)/12)*12)./12)...
    - fix(3.*((Ye+4900+fix((Mo-14)/12))./100)./4);
JD = J-0.5+Hr./24.0+Mi./1440.0+Se./86400.0;
MJD=JD-2400000.5;

