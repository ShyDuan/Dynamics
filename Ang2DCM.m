function A = Ang2DCM(bi,bj,bk,i,j,k)
if i==1; Ai=Rx(bi);end
if i==2; Ai=Ry(bi);end
if i==3; Ai=Rz(bi);end
if j==1; Aj=Rx(bj);end
if j==2; Aj=Ry(bj);end
if j==3; Aj=Rz(bj);end
if k==1; Ak=Rx(bk);end
if k==2; Ak=Ry(bk);end
if k==3; Ak=Rz(bk);end
A=Ak*Aj*Ai;

function A=Rx(x)
A=[1   0    0;
   0  cos(x) sin(x);
   0 -sin(x) cos(x)];

function A=Ry(x)
A=[cos(x) 0  -sin(x);
    0     1     0;
   sin(x) 0  cos(x)];

function A=Rz(x)
A=[cos(x)  sin(x) 0;
  -sin(x)  cos(x) 0;
     0        0   1];