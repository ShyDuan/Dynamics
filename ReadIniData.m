 
function  [OrbEpoc,OrbRV,Hw0,SatAtt0]=ReadIniData
fdata=fopen('SimIniData.m');
if (fdata==-1)
    display('Cannot Open The File, Check!')
    return
end
fgets(fdata);
fgets(fdata);
OrbEpoc=fscanf(fdata,'%g',6);
fgets(fdata);
fgets(fdata);
OrbEle=fscanf(fdata,'%g',6);
OrbRV=Ele2RV(OrbEle);

fgets(fdata);
fgets(fdata);
Hw0=fscanf(fdata,'%g',4);

fgets(fdata);
fgets(fdata);
SatAtt0=fscanf(fdata,'%g',6);


