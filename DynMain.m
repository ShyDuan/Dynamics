function [RV]=DynMain(Twm,MT,Tthrust)
global lbdG0;
d2r=pi/180;
para.Jsat=1e-6*diag([82891;  111375;  120296*1]);%�������
para.Cw=[1 0 0  0.57735; 
         0 1 0  0.57735;
         0 0 1  0.57735];%���Ӱ�װ��
Hwb=[0.01;0.01; 0.01; -0.01732];%����Ŀ��Ƕ���
CDssb=[1  0  0
       0  -1  0
       0  0  -1]; %̫����װ��
result=[];
%������������ļ�
 [OrbEpoc,OrbRV,Hw0,SatAtt0]=ReadIniData; 
 tnow=0;h=0.001;Acc=[0;0;0];Bo=[1;0;0]; lbdG0 = LGcal(OrbEpoc); 
Hw=Hw0;
while tnow<100
    %=======���Ͽ��ƣ��򻯣�=============
    Kp=[0.008; 0.011;0.012];  
    Kd=[0.037;0.049;0.053];
    Ki=[0.00018;0.00024;0.00026];
    Up=SatAtt0(1:3)*d2r;
    Ud=SatAtt0(4:6)*d2r;    
    Tc=Kp.*Up+Kd.*Ud;
    Twm=pinv(para.Cw)*Tc; 
    
    Cbo=Ang2DCM(SatAtt0(3)*d2r,SatAtt0(2)*d2r,SatAtt0(1)*d2r,3,2,1);
    Bb=Cbo*Bo;
    Kpm=1e-4;  
    dH=para.Cw*(Hw0-Hwb);
    MT = Kpm*cross(dH,Bb)/norm(Bb)^2;
    Tmt=cross(MT,Bb);
    %=====================
    
    %����ִ������������Twm�ʹ�������MT
    
    % �������
    RV=odesolver('OrbDyn',tnow,OrbRV,0.25,h,Acc,[]);   OrbRV=RV;
    ELE=RV2Element(RV);
    Si=Sun(OrbEpoc,tnow);
    Bo=EarthMagField(tnow,ELE);
    bShadow=EarthShadow(Si,RV(1:3));
    
    %��̬����������
    Tthrust=[0;0;0]; 
    
    %���Ӷ���ѧ
    for i=1:4
        if abs(Hw(i))>=0.015; Twm(i)=0;end
        Hw(i)=odesolver('dHw',tnow,Hw0(i),0.25,h,Twm(i),[]);    
    end
    Hw0=Hw;
    
    %���������������Ĵž�MT    
     
    para.Hw=para.Cw*Hw;
    para.w0=ELE.w0;
    Tc=-para.Cw*Twm + Tmt + Tthrust;
    % ��̬����ѧ
    SatAtt=odesolver('AttDyn',tnow,SatAtt0,0.25,h,Tc,para);
    SatAtt0=SatAtt;
    
    
    % ���������
    GnssNoise = [0.1;0.1;0.1;0.001;0.001;0.001]*randn(1);
    StNoise = [5;5;70]/3600*randn(1); %deg
    GANoise=0.09/60/(0.25)^0.5*randn(1);%deg
    GAb=[10;50;20]/3600;
    MMnoise=300e-9*randn(1);
    DssNoise=[0.5;0.5].*rand(2,1);%deg
    
    GNSSout=GNSSmodel(RV,tnow,GnssNoise);
    STout=STmodel(SatAtt(1:3),ELE.Coi,StNoise);
    GAout=GAmodel(SatAtt(4:6),GANoise,GAb);
    DSSout=DSSmodel(Cbo*ELE.Coi*Si,bShadow,CDssb,DssNoise);
    MMout=MMmodel(Bb,MMnoise);    
      tnow=tnow+0.25, %RV
    result=[result,[tnow;SatAtt;Hw0]];
end

figure(1)
plot(result(1,:),result(2:7,:));grid on
figure(2)
plot(result(1,:),result(8:11,:));grid on



function GNSSout=GNSSmodel(RV,tgps,noise)
GNSSout.rv=RV+noise;
GNSSout.Tgps=tgps;

function GAout=GAmodel(wbi,noise,b)
GAout=wbi+ b+noise; 

function STout=STmodel(ang,Coi,noise)
Csb=eye(3); %������װ��,�ݶ�
Cbo=Ang2DCM(ang(3)*pi/180,ang(2)*pi/180,ang(1)*pi/180,3,2,1);
Css=eye(3)+[0 noise(3) -noise(2);-noise(3) 0 noise(1);noise(2) -noise(1) 0]*pi/180;
STout=C2Q(Css*Csb*Cbo*Coi);

function DSSout=DSSmodel(Sb,bShadow,CDssb,noise)
Ss=CDssb*Sb;   
S61=   atan(Ss(2)/Ss(3));  %beta 
S62=  +atan(Ss(1)/Ss(3));  %alf
bDssValid = 3;
if S61>60*pi/180 || S61<-60*pi/180 || S62>60*pi/180 || S62<-60*pi/180 || bShadow==1
    S61=0;S62=0; bDssValid=0;
end
DSSout=[[S61;S62]*180/pi+noise;bDssValid]; % deg


function MMout=MMmodel(Bb,noise)
MMout=Bb+noise;

