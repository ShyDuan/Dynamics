 function RV=Ele2RV(Ele1)
d2r=pi/180;mui=398600.5;    %(km^3/s^2)
a=Ele1(1);e=Ele1(2);i=Ele1(3)*d2r;OMEGA=Ele1(4)*d2r;omega=Ele1(5)*d2r;M=Ele1(6)*d2r;
E=M;
for counter1=1:100
    E=M+e*sin(E);
end
cosE=cos(E);sinE=sin(E);
rcosf=a*(cosE-e);
rsinf=a*sqrt(1-e^2)*sinE;
r=a*(1-e*cosE);
if rsinf~=0
    sinf=rsinf/r;
    f=asin(sinf);
    if rcosf<0
        if  rsinf>=0
            f=pi-f;
        else
            f=-1*pi-f;
        end
    end
else if rcosf>0
        f=0;
    else
        f=pi;
    end
end
 
u=f+omega;
%R
r=a*(1-e*cosE);
rx=r*(cos(u)*cos(OMEGA)-sin(u)*sin(OMEGA)*cos(i));
ry=r*(cos(u)*sin(OMEGA)+sin(u)*cos(OMEGA)*cos(i));
rz=r*sin(i)*sin(u);
%V
p=a*(1-e^2);
kv=sqrt(mui/p);
vx=kv*(-cos(OMEGA)*sin(u)-sin(OMEGA)*cos(i)*cos(u)-cos(OMEGA)*e*sin(omega)-sin(OMEGA)*cos(i)*e*cos(omega));
vy=kv*(-sin(OMEGA)*sin(u)+cos(OMEGA)*cos(i)*cos(u)-sin(OMEGA)*e*sin(omega)+cos(OMEGA)*cos(i)*e*cos(omega));
vz=kv*(sin(i)*cos(u)+e*cos(omega)*sin(i));
R=[rx,ry,rz];
V=[vx,vy,vz];
RV=[R';V'];
