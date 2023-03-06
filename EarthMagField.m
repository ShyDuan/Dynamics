function Bo = EarthMagField(time,Ele)
global lbdG0
u=Ele.omg+Ele.f;r=Ele.a*(1-Ele.e^2)/(1+Ele.e*cos(Ele.f)); we=7.2921159e-5;  

L1=atan(tan(u)*cos(Ele.i));
if cos(u)<0;L1=pi+atan(tan(u)*cos(Ele.i));end;
L=L1+Ele.OMG-(lbdG0+we*time);
dlta=atan(sin(L1)*tan(Ele.i));

% start
sita=pi/2-dlta;

P01=cos(sita);  P02=sin(sita);
P12=sqrt(3)*cos(sita);
P22=sqrt(3)/2*sin(sita);
P03=2.5*cos(sita)*((cos(sita))^2-0.6);  P13=5*sqrt(6)/4*((cos(sita))^2-0.2);
P23=sqrt(15)/4*sin(2*sita);  P33=sqrt(10)/4*(sin(sita))^2;
P04=35/8*((cos(sita))^2-6/7*(cos(sita))^2+3/35);
dP01=-sin(sita);  dP11=cos(sita);
dP02=-1.5*sin(2*sita);  dP12=sqrt(3)*cos(2*sita);
dP22=sqrt(3)/2*sin(2*sita);
dP03=7.5*sin(sita)*(0.2-(cos(sita))^2);  dP13=5*sqrt(6)/4*cos(sita)*(0.8-3*(sin(sita))^2);
dP23=sqrt(15)/2*sin(sita)*(3*(cos(sita))^2-1);    dP33=3*sqrt(10)/8*sin(sita)*sin(2*sita);
dP04=35/4*(-cos(sita)*cos(sita)+3/7)*sin(2*sita);

Re=6378.137;
g01=-29431.7; g11=-1482.9; %nT,2016
g02=-2453.8;   g12=3009.6;  g22=1678.8;
g03=1354.1;    g13=-2357.8; g23=1224.9;  g33=571.9;
g04=906.9;   %g14=813.9;
h11=4770.5;
h12=-2873;h22=-656;

Bs1=-(Re/r)^3*(g01*dP01+(g11*cos(L)+h11*sin(L))*dP11);
Bf1=-(Re/r)^3*(-g11*sin(L)+h11*cos(L));
Br1=2*(Re/r)^3*(g01*P01+(g11*cos(L)+h11*sin(L))*sin(sita));

Bs2=-(Re/r)^4*(g02*dP02+(g12*cos(L)+h12*sin(L))*dP12+(g22*cos(2*L)+h22*sin(2*L))*dP22);
Bf2=-(Re/r)^4*((-g12*sin(L)+h12*cos(L))*P12+(-g22*sin(2*L)+h22*cos(2*L))*2*P22);
Br2=3*(Re/r)^4*(g02*P02+(g12*cos(L)+h12*sin(L))*P12*sin(sita))+(g22*cos(2*L)+h22*sin(2*L))*P22*sin(sita)

Bs3=-(Re/r)^5*(g03*dP03+g13*cos(L)*dP13+g23*cos(2*L)*dP23+g33*cos(3*L)*dP33);
Bf3=-(Re/r)^5*(-g13*sin(L)*P13*1-g23*sin(2*L)*P23*2-g33*sin(3*L)*P33*3);
Br3=4*(Re/r)^5*(g03*P03+g13*P03+g13*cos(L)*P13*sin(sita)+g23*cos(2*L)*P23*sin(sita)+g33*cos(3*L))*P33*sin(sita);

Bs=Bs1+Bs2+Bs3-(Re/r)^6*g04*dP04;
Bf=Bf1+Bf2+Bf3;
Br=Br1+Br2+Br3+5*(Re/r)^6*g04*P04;
B=sqrt(Bs^2+Bf^2+Br^2);

if Ele.i>=pi/2
    cosks=-sin(Ele.i)*cos(L1);  sinks=sqrt(1-cosks*cosks);
    Box=Bs*cosks-Bf*sinks;
    Boy=-Bs*sinks-Bf*cosks;
    Boz=-Br;
else
    cosks=sin(Ele.i)*cos(L1);  sinks=sqrt(1-cosks*cosks);
    Box=-Bs*cosks+Bf*sinks;
    Boy= Bs*sinks+Bf*cosks;
    Boz=-Br;
end
Bo=[Box;Boy;Boz]*1e-9;  %Unit nT->T