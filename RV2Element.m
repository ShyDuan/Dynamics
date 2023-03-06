function Ele=RV2Element(OrbRV)
mui=398600.5;
x=OrbRV(1);y=OrbRV(2);z=OrbRV(3);vx=OrbRV(4);vy=OrbRV(5);vz=OrbRV(6);
r=(x^2+y^2+z^2)^0.5;
V=(vx^2+vy^2+vz^2)^0.5;
a=1/(2/r-V^2/mui);
RdotV=x*vx+y*vy+z*vz;
ex=((V^2-mui/r)*x-RdotV*vx)/mui;
ey=((V^2-mui/r)*y-RdotV*vy)/mui;
ez=((V^2-mui/r)*z-RdotV*vz)/mui;
e=(ex^2+ey^2+ez^2)^0.5;
hx=y*vz-z*vy;
hy=z*vx-x*vz;
hz=x*vy-y*vx;
h=(hx^2+hy^2+hz^2)^0.5;
Sn=(hx^2+hy^2)^0.5;
i=acos(hz/h);
if abs(i)<1e-7
  OMG=0;
  u=atan2(y,x);
else
  if hx>=0.0
    OMG=acos(-hy/Sn);
  else
    OMG=-acos(-hy/Sn);
  end
  sinu=z/(r*sin(i));
  cosu=(x*cos(OMG)+y*sin(OMG))/r;
  u=atan2(sinu,cosu);
end
E=atan2(RdotV/(mui*a)^0.5,1-r/a);
Esf=(1-e^2)^0.5*sin(E);
Ecf=cos(E)-e;
f=atan2(Esf,Ecf);
omg=u-f; 
omg=mod(omg,2*pi);
M=E-e*sin(E);
lbd=omg+M;
lbd=mod(lbd,2*pi);
Ele.a=a; 
Ele.e=e;
Ele.i=i;
Ele.OMG=OMG;
Ele.omg=omg;
Ele.M=M;
Ele.f=f;
Ele.lbd=lbd;
Ele.w0=sqrt(mui/(a*(1-e^2))^3)*(1+e*cos(f))^2;
Ele.Coi=[-cos(OMG)*sin(u)-sin(OMG)*cos(i)*cos(u) -sin(OMG)*sin(u)+cos(OMG)*cos(i)*cos(u) sin(i)*cos(u); ...
        -sin(OMG)*sin(i) cos(OMG)*sin(i) -cos(i); ...
        -cos(OMG)*cos(u)+sin(OMG)*cos(i)*sin(u) -sin(OMG)*cos(u)-cos(OMG)*cos(i)*sin(u) -sin(i)*sin(u)];
