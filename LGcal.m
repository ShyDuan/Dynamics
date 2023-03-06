
function lbdG = LGcal(OrbiniEpoc)
Y=OrbiniEpoc(1);
M=OrbiniEpoc(2);
D=OrbiniEpoc(3);
H=OrbiniEpoc(4);
Mi=OrbiniEpoc(5);
S=OrbiniEpoc(6);

n = datenum(Y,M,D,H,Mi,S); %
POD=rem(n,1);  
JD=367*Y-floor(7*(Y+floor((M+9)/12))/4) ...
   -floor(3*(floor((Y+(M-9)/7)/100)+1)/4) ...
   +floor(275*M/9)+D+1721028.5+POD;  % 
d = JD -  2451545.0; 
lbdG = mod(280.4606184 + 360.9856122863*d,360); 
lbdG = lbdG*pi/180;