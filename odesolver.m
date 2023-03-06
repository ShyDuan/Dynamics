function yn1 = odesolver(DYN,tn,yn,Tspan,h,input,para)

fyt=str2func(DYN);
for i=1:Tspan/h
    K1=fyt(tn,yn,input,para); 
    K2=fyt(tn+h/2,yn+0.5*h*K1,input,para);
    K3=fyt(tn+h/2,yn+0.5*h*K2,input,para);
    K4=fyt(tn+h,yn+h*K3,input,para);
    yn1=yn+h*(K1+2*K2+2*K3+K4)/6;
    tn=tn+h;
    yn=yn1;  
end


function dRV=OrbDyn(t,rv,a,para)  %a=F/m
mui=398600.5; Re=6378.137;
J2=0.00108263;  J3=-2.541000725156921e-6;   J4=-1.6185636e-6;
rx=rv(1);ry=rv(2);rz=rv(3);
vx=rv(4);vy=rv(5);vz=rv(6);
ax=a(1);ay=a(2);az=a(3);
r0=sqrt(rx^2+ry^2+rz^2);
v0=sqrt(vx^2+vy^2+vz^2);
drx=vx;
dry=vy;
drz=vz;
dvx=-(mui*rx/r0^3)*(1+1.5*J2*(Re/r0)^2*(1-5*(rz/r0)^2)+2.5*J3*(Re/r0)^3*(3*rz/r0-7*(rz/r0)^3)-0.625*J4*(Re/r0)^4*(3-42*(rz/r0)^2+63*(rz/r0)^4))+ax;
dvy=-(mui*ry/r0^3)*(1+1.5*J2*(Re/r0)^2*(1-5*(rz/r0)^2)+2.5*J3*(Re/r0)^3*(3*rz/r0-7*(rz/r0)^3)-0.625*J4*(Re/r0)^4*(3-42*(rz/r0)^2+63*(rz/r0)^4))+ay;
dvz=-(mui*rz/r0^3)*(1+1.5*J2*(Re/r0)^2*(3-5*(rz/r0)^2)-0.625*J4*(Re/r0)^4*(15-70*(rz/r0)^2+63*(rz/r0)^4))-2.5*J3*(Re/r0)^3*(6*(rz/r0)^2-7*(rz/r0)^4-0.6)*(mui/r0^2)+az;

dRV=[drx;dry;drz;dvx;dvy;dvz];

function dwda = AttDyn(t,x,Tc,para)
fai=x(1)*pi/180;sita=x(2)*pi/180; psi=x(3)*pi/180;
wx=x(4)*pi/180;wy=x(5)*pi/180;wz=x(6)*pi/180; w=[wx;wy;wz];

J=para.Jsat;
Hw=para.Hw;
w0=para.w0;

dw= inv(J)*(-cross(w,J*w+Hw)+Tc);
dwx=dw(1);dwy=dw(2);dwz=dw(3);

Cbo =[cos(sita)*cos(psi)                                       cos(sita)*sin(psi)                               -sin(sita);...
     -cos(fai)*sin(psi)+sin(fai)*sin(sita)*cos(psi)    cos(fai)*cos(psi)+sin(fai)*sin(sita)*sin(psi)   sin(fai)*cos(sita);...
      sin(fai)*sin(psi)+cos(fai)*sin(sita)*cos(psi)   -sin(fai)*cos(psi)+cos(fai)*sin(sita)*sin(psi)   cos(fai)*cos(sita)];
wr = [wx;wy;wz] + Cbo*[0;w0;0];

    dfai = wr(1)+tan(sita)*(wr(2)*sin(fai)+wr(3)*cos(fai));
    dsita = wr(2)*cos(fai)-wr(3)*sin(fai);
    dpsi=(wr(2)*sin(fai)+wr(3)*cos(fai))/cos(sita);
    
    dwda=[dfai;dsita;dpsi;dwx;dwy;dwz]*180/pi;
    
function  dHwheel= dHw(t,Hwi,Twi,para)
Twfi = 1e-4*sign(Hwi)+1e-3*Hwi;
dHwheel=+Twi-Twfi;     








