%digital controlled effect of performance on the inverter
clear all; close all; clc

s=tf('s');

R=0.62;
L=1.22*(10^-3);
C=20*(10^-6);
wn=sqrt(1/(L*C));
kxi=R/(2*L*wn);
td=982.24*(10^(-6));

kk=[0.1:0.1:2]; %1x20
ts=kk*td;

G_fcns=wn^2/(s^2+(2*kxi*wn)*s+wn^2);
T_fcns=feedback(G_fcns,1);

%varies sample period
for pp=1:length(ts)
    G_fcnz{1,pp}=c2d(G_fcns,ts(pp),'zoh');
    T_fcnz{1,pp}=feedback(G_fcnz{1,pp},1);
    [numz{1,pp},denz{1,pp}]=tfdata(T_fcnz{1,pp});
    [numw{1,pp},denw{1,pp}]=z2wtrans(numz{1,pp},denz{1,pp},ts(pp));
    T_fcnw{1,pp}=tf(numw{1,pp},denw{1,pp});
    G_fcnw{1,pp}=d2c(G_fcnz{1,pp},'tustin');
    
end

figure(11)
bode(G_fcns,G_fcnz{1,1},G_fcnz{1,2},G_fcnz{1,3},G_fcnz{1,4},G_fcnz{1,5})
legend('continuous','0.1T_d','0.2T_d','0.3T_d','0.4T_d','0.5T_d','Location','SouthWest')
figure(12)
bode(G_fcnz{1,5},G_fcnz{1,6},G_fcnz{1,7},G_fcnz{1,8},G_fcnz{1,9},G_fcnz{1,10})
legend('0.5T_d','0.6T_d','0.7T_d','0.8T_d','0.9T_d','1.0T_d','Location','SouthWest')
figure(13)
bode(G_fcnz{1,10},G_fcnz{1,11},G_fcnz{1,12},G_fcnz{1,13},G_fcnz{1,14},G_fcnz{1,15})
legend('1.0T_d','1.1T_d','1.2T_d','1.3T_d','1.4T_d','1.5T_d','Location','SouthWest')
figure(14)
bode(G_fcnz{1,15},G_fcnz{1,16},G_fcnz{1,17},G_fcnz{1,18},G_fcnz{1,19},G_fcnz{1,20})
legend('1.5T_d','1.6T_d','1.7T_d','1.8T_d','1.9T_d','2.0T_d','Location','SouthWest')

figure(21)
bode(G_fcns,G_fcnw{1,1},G_fcnw{1,2},G_fcnw{1,3},G_fcnw{1,4},G_fcnw{1,5})
legend('continuous','0.1T_d','0.2T_d','0.3T_d','0.4T_d','0.5T_d','Location','SouthEast')
figure(22)
bode(G_fcnw{1,5},G_fcnw{1,6},G_fcnw{1,7},G_fcnw{1,8},G_fcnw{1,9},G_fcnw{1,10})
legend('0.5T_d','0.6T_d','0.7T_d','0.8T_d','0.9T_d','1.0T_d','Location','NorthEast')
figure(23)
bode(G_fcnw{1,10},G_fcnw{1,11},G_fcnw{1,12},G_fcnw{1,13},G_fcnw{1,14},G_fcnw{1,15})
legend('1.0T_d','1.1T_d','1.2T_d','1.3T_d','1.4T_d','1.5T_d','Location','NorthEast')
figure(24)
bode(G_fcnw{1,15},G_fcnw{1,16},G_fcnw{1,17},G_fcnw{1,18},G_fcnw{1,19},G_fcnw{1,20})
legend('1.5T_d','1.6T_d','1.7T_d','1.8T_d','1.9T_d','2.0T_d','Location','NorthEast')

[Gms,Pms,Wgs,Wps]=margin(G_fcns);

for pp=1:length(ts)
    [Gm(pp),Pm(pp),Wg(pp),Wp(pp)]=margin(G_fcnw{1,pp});
end

qwe=0:0.1:2;
Gm=[Gms Gm];
Pm=[Pms Pm];
figure(3)
plot(qwe,Gm,'o');


figure(4)
plot(qwe,Pm,'o');

