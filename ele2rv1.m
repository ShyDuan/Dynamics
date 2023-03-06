
function RV=ele2rv(Ele1)
mui=398600.5;
a=Ele1(1);e=Ele1(2);i=Ele1(3);OMEGA=Ele1(4);omega=Ele1(5);f=Ele1(6);
u=f+omega;
%R
p=a*(1-e^2);
r=p/(1+e*cos(f));
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


