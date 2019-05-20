clear all; close all; clc;

s=tf('s');
z=tf('z');

R=0.62;
RL=10;
L=1.22*(10^-3);
C=20*(10^-6);
Utp=500; %magnitude of high freqency triangle carrier wave

wn(1)=sqrt(1/(L*C)); %no load
wn(2)=sqrt((R+RL)/(L*C*RL));
kxi(1)=(R*C)/(2*L*C*wn(1)); %no load
kxi(2)=(R*RL*C+L)/(2*L*C*RL*wn(2));
k_L=[1 RL/(R+RL)];

k=4000/4096;
ku=8*4;

tc=100*10^-6;
tt=[1:0.5:25];
ts=tt*tc;

tt2=[1:0.5:25];
ts2=tt2*tc;

kc=0.1:0.1:5;
kc2=0.01:0.01:10;

for hh=1:length(ts)
    
    opn_fns1{1,hh}=(wn(2)^2)/(s^2+2*kxi(2)*wn(2)*s+wn(2)^2);
    opn_fnz1{1,hh}=c2d(opn_fns1{1,hh},ts(hh),'zoh');
    opn_fnz_stpdy1{1,hh}=c2d(opn_fns1{1,hh},ts(hh),'zoh')*(1/z);
        
    r11=rlocus(opn_fnz1{1,hh},kc);
    rho11=abs(r11)-1;
    [m,n]=size(rho11);
    zero_index111=find(rho11(1,1:(end-1)).*rho11(1,(2:end))<0);
    zero_index112=find(rho11(2,1:(end-1)).*rho11(2,(2:end))<0);
    zero_index11=min([zero_index111 zero_index112]);
    stb_gain11(1,hh)=zero_index11(1,1)*0.1;
    
    r12=rlocus(opn_fnz_stpdy1{1,hh},kc);
    rho12=abs(r12)-1;
    [m,n]=size(rho12);
    zero_index121=find(rho12(1,1:(end-1)).*rho12(1,(2:end))<0);
    zero_index122=find(rho12(2,1:(end-1)).*rho12(2,(2:end))<0);
    zero_index12=min([zero_index121 zero_index122]);
    stb_gain12(1,hh)=zero_index12(1,1)*0.1;

%     S11{1,hh}=allmargin(opn_fnz1{1,hh});
%     S12{1,hh}=allmargin(opn_fnz_stpdy1{1,hh});
%     Gm11(hh)=S11{1,hh}.GainMargin(1,1); %gain margin
%     Gm12(hh)=S12{1,hh}.GainMargin(1,1);
    clear r11 r12 rho11 rho12 zero_index11 zero_index_12
end

for ff=1:length(ts2)
    
    opn_fns2{1,ff}=(wn(1)^2)/(s^2+2*kxi(1)*wn(1)*s+wn(1)^2);
    opn_fnz2{1,ff}=c2d(opn_fns2{1,ff},ts2(ff),'zoh');
    opn_fnz_stpdy2{1,ff}=c2d(opn_fns2{1,ff},ts2(ff),'zoh')*(1/z);
   
    r21=rlocus(opn_fnz2{1,ff},kc2);
    rho21=abs(r21)-1;
    [m,n]=size(rho21);
    zero_index211=find(rho21(1,1:(end-1)).*rho21(1,(2:end))<0);
    zero_index212=find(rho21(2,1:(end-1)).*rho21(2,(2:end))<0);
    zero_index21=min([zero_index211 zero_index212]);
    stb_gain21(1,ff)=zero_index21(1,1)*0.01;
    
    r22=rlocus(opn_fnz_stpdy2{1,ff},kc2);
    rho22=abs(r22)-1;
    [m,n]=size(rho22);
    zero_index221=find(rho22(1,1:(end-1)).*rho22(1,(2:end))<0);
    zero_index222=find(rho22(2,1:(end-1)).*rho22(2,(2:end))<0);
    zero_index22=min([zero_index221 zero_index222]);
    stb_gain22(1,ff)=zero_index22(1,1)*0.01;
    

%     S21{1,ff}=allmargin(opn_fnz2{1,ff});
%     S22{1,ff}=allmargin(opn_fnz_stpdy2{1,ff});
%     Gm21(ff)=S21{1,ff}.GainMargin(1,1);
%     Gm22(ff)=S22{1,ff}.GainMargin(1,1);
    clear r21 r22 rho21 rho22 zero_index21 zero_index22
end

figure(1)
plot(ts,stb_gain11,'b:.'); hold on;
title('Stability of the one-step-delay and no delay system. R_L = 10\Omega.')
xlabel('T_s(sec)'); ylabel('K_c');
plot(ts,stb_gain12,'r:*');
axis([10^-4 2.5*10^-3 0 3])
legend('R_L = 10\Omega no delay','R_L = 10\Omega one-step delay' )

figure(2)
plot(ts2,stb_gain21,'b:.'); hold on;
title('Stability of the one-step-delay and no delay system. R_L = \infty.')
xlabel('T_s(sec)'); ylabel('K_c');
plot(ts2,stb_gain22,'r:*');
axis([10^-4 2.5*10^-3 0 9])
legend('R_L = \infty, no delay','R_L = \infty, one-step delay')
