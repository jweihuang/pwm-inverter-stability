clear all; close all; clc;

s=tf('s');
z=tf('z');

R=0.62;
L=1.22*(10^-3);
C=20*(10^-6);

wn=[3400 6400 9400];
kxi=0.0397;
td=[(2*pi/3400) (2*pi/6400) (2*pi/9400)];

tt1=[0.01:0.01:4/(2*pi*10^3/3400)];
ts1=tt1.*td(1);

tt2=[0.01:0.01:4/(2*pi*10^3/6400)];
ts2=tt2.*td(2);

tt3=[0.01:0.01:4/(2*pi*10^3/9400)];
ts3=tt3.*td(3);

wid=0.1;
gain=0:wid:10;

for hh=1:length(ts1)
    
    opn_fns1{1,hh}=(wn(1)^2)/(s^2+2*kxi*wn(1)*s+wn(1)^2);
    opn_fnz1{1,hh}=c2d(opn_fns1{1,hh},ts1(hh),'zoh');
    opn_fnz_stpdy1{1,hh}=c2d(opn_fns1{1,hh},ts1(hh),'zoh')*(1/z);
    clo_fnz1{1,hh}=feedback(opn_fnz1{1,hh},1);
    clo_fnz_stpdy1{1,hh}=feedback(opn_fnz_stpdy1{1,hh},1);
%     opn_fnw1{1,hh}=d2c(opn_fnz1{1,hh}, 'tustin'); %bilinear transform w-domain
%     opn_fnw_stpdy1{1,hh}=d2c(opn_fnz_stpdy1{1,hh}, 'tustin');
    
    r11=rlocus(opn_fnz1{1,hh},gain);
    rho11=abs(r11)-1;
    [m,n]=size(rho11);
    zero_index111=find(rho11(1,1:(end-1)).*rho11(1,(2:end))<0);
    zero_index112=find(rho11(2,1:(end-1)).*rho11(2,(2:end))<0);
    zero_index11=min([zero_index111 zero_index112]);
    stb_gain11(1,hh)=zero_index11(1,1)*wid;
    
    r12=rlocus(opn_fnz_stpdy1{1,hh},gain);
    rho12=abs(r12)-1;
    [m,n]=size(rho12);
    zero_index121=find(rho12(1,1:(end-1)).*rho12(1,(2:end))<0);
    zero_index122=find(rho12(2,1:(end-1)).*rho12(2,(2:end))<0);
    zero_index12=min([zero_index121 zero_index122]);
    stb_gain12(1,hh)=zero_index12(1,1)*wid;
    
%     [Gm11(hh),Pm11(hh),Wg11(hh),Wp11(hh)]=margin(opn_fnz1{1,hh});
%     [Gm12(hh),Pm12(hh),Wg12(hh),Wp12(hh)]=margin(opn_fnz_stpdy1{1,hh});
    S11{1,hh}=allmargin(opn_fnz1{1,hh});
    S12{1,hh}=allmargin(opn_fnz_stpdy1{1,hh});
    Gm11(hh)=S11{1,hh}.GainMargin(1,1); %gain margin
    Gm12(hh)=S12{1,hh}.GainMargin(1,1);
    clear r11 r12 rho11 rho12 zero_index11 zero_index_12
end

for ff=1:length(ts2)
    
    opn_fns2{1,ff}=(wn(2)^2)/(s^2+2*kxi*wn(2)*s+wn(2)^2);
    opn_fnz2{1,ff}=c2d(opn_fns2{1,ff},ts2(ff),'zoh');
    opn_fnz_stpdy2{1,ff}=c2d(opn_fns2{1,ff},ts2(ff),'zoh')*(1/z);
    clo_fnz2{1,ff}=feedback(opn_fnz2{1,ff},1);
    clo_fnz_stpdy2{1,ff}=feedback(opn_fnz_stpdy2{1,ff},1);
%     opn_fnw2{1,ff}=d2c(opn_fnz2{1,ff}, 'tustin'); %bilinear transform w-domain
%     opn_fnw_stpdy2{1,ff}=d2c(opn_fnz_stpdy2{1,ff}, 'tustin');
    
    r21=rlocus(opn_fnz2{1,ff},gain);
    rho21=abs(r21)-1;
    [m,n]=size(rho21);
    zero_index211=find(rho21(1,1:(end-1)).*rho21(1,(2:end))<0);
    zero_index212=find(rho21(2,1:(end-1)).*rho21(2,(2:end))<0);
    zero_index21=min([zero_index211 zero_index212]);
    stb_gain21(1,ff)=zero_index21(1,1)*wid;
    
    r22=rlocus(opn_fnz_stpdy2{1,ff},gain);
    rho22=abs(r22)-1;
    [m,n]=size(rho22);
    zero_index221=find(rho22(1,1:(end-1)).*rho22(1,(2:end))<0);
    zero_index222=find(rho22(2,1:(end-1)).*rho22(2,(2:end))<0);
    zero_index22=min([zero_index221 zero_index222]);
    stb_gain22(1,ff)=zero_index22(1,1)*wid;
    
%     [Gm21(ff),Pm21(ff),Wg21(ff),Wp21(ff)]=margin(opn_fnz2{1,ff});
%     [Gm22(ff),Pm22(ff),Wg22(ff),Wp22(ff)]=margin(opn_fnz_stpdy2{1,ff});
    S21{1,ff}=allmargin(opn_fnz2{1,ff});
    S22{1,ff}=allmargin(opn_fnz_stpdy2{1,ff});
    Gm21(ff)=S21{1,ff}.GainMargin(1,1);
    Gm22(ff)=S22{1,ff}.GainMargin(1,1);
    clear r21 r22 rho21 rho22 zero_index21 zero_index22
end


for ss=1:length(ts3)
    
    opn_fns3{1,ss}=(wn(3)^2)/(s^2+2*kxi*wn(3)*s+wn(3)^2);
    opn_fnz3{1,ss}=c2d(opn_fns3{1,ss},ts3(ss),'zoh');
    opn_fnz_stpdy3{1,ss}=c2d(opn_fns3{1,ss},ts3(ss),'zoh')*(1/z);
    clo_fnz3{1,ss}=feedback(opn_fnz3{1,ss},1);
    clo_fnz_stpdy3{1,ss}=feedback(opn_fnz_stpdy3{1,ss},1);
%     opn_fnw3{1,ss}=d2c(opn_fnz3{1,ss}, 'tustin'); %bilinear transform w-domain
%     opn_fnw_stpdy3{1,ss}=d2c(opn_fnz_stpdy3{1,ss}, 'tustin');
    
    r31=rlocus(opn_fnz3{1,ss},gain);
    rho31=abs(r31)-1;
    [m,n]=size(rho31);
    zero_index311=find(rho31(1,1:(end-1)).*rho31(1,(2:end))<0);
    zero_index312=find(rho31(2,1:(end-1)).*rho31(2,(2:end))<0);
    zero_index31=min([zero_index311 zero_index312]);
    stb_gain31(1,ss)=zero_index31(1,1)*wid;
    
    r32=rlocus(opn_fnz_stpdy3{1,ss},gain);
    rho32=abs(r32)-1;
    [m,n]=size(rho32);
    zero_index321=find(rho32(1,1:(end-1)).*rho32(1,(2:end))<0);
    zero_index322=find(rho32(2,1:(end-1)).*rho32(2,(2:end))<0);
    zero_index32=min([zero_index321 zero_index322]);
    stb_gain32(1,ss)=zero_index32(1,1)*wid;
    
%     [Gm31(ss),Pm31(ss),Wg31(ss),Wp31(ss)]=margin(opn_fnz3{1,ss});
%     [Gm32(ss),Pm32(ss),Wg32(ss),Wp32(ss)]=margin(opn_fnz_stpdy3{1,ss});
    S31{1,ss}=allmargin(opn_fnz3{1,ss});
    S32{1,ss}=allmargin(opn_fnz_stpdy3{1,ss});
    Gm31(ss)=S31{1,ss}.GainMargin(1,1);
    Gm32(ss)=S32{1,ss}.GainMargin(1,1);

    clear r31 r32 rho31 rho32 zero_index31 zero_index_32
end

figure(10)
plot(ts1,stb_gain11,'b:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\omega_n = 3400 rad/s, 6400 rad/s, 9400 rad/s.')
xlabel('T_s(sec)'); ylabel('K');
plot(ts1,stb_gain12,'b:*'); hold on;
plot(ts2,stb_gain21,'g:.'); hold on;
plot(ts2,stb_gain22,'g:*'); hold on;
plot(ts3,stb_gain31,'r:.'); hold on;
plot(ts3,stb_gain32,'r:*')
axis([0 0.004 0 10])
legend('\omega_n = 3400 rad/s no step delay','\omega_n = 3400 rad/s one-step delay', '\omega_n = 6400 rad/s no step delay','\omega_n = 6400 rad/s one-step delay','\omega_n = 9400 rad/s no step delay','\omega_n = 9400 rad/s one-step delay' )

figure(11)
plot(ts1,stb_gain11,'b:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\omega_n = 3400 rad/s.')
xlabel('T_s(sec)'); ylabel('K');
plot(ts1,stb_gain12,'b:*'); 
axis([0 0.004 0 10])
legend('\omega_n = 3400 rad/s no step delay','\omega_n = 3400 rad/s one-step delay' )

figure(12)
plot(ts2,stb_gain21,'g:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\omega_n = 6400 rad/s.')
xlabel('T_s(sec)'); ylabel('K');
plot(ts2,stb_gain22,'g:*');
axis([0 0.004 0 10])
legend('\omega_n = 6400 rad/s no step delay','\omega_n = 6400 rad/s one-step delay' )

figure(13)
plot(ts3,stb_gain31,'r:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\omega_n = 9400 rad/s.')
xlabel('T_s(sec)'); ylabel('K');
plot(ts3,stb_gain32,'r:*')
axis([0 0.004 0 10])
legend('\omega_n = 9400 rad/s no step delay','\omega_n = 9400 rad/s one-step delay' )

figure(20)
plot(ts1,Gm11,'b:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\omega_n = 3400 rad/s, 6400 rad/s, 9400 rad/s.')
xlabel('T_s(sec)'); ylabel('absolute magnitude');
plot(ts1,Gm12,'b:*'); hold on;
plot(ts2,Gm21,'g:.'); hold on;
plot(ts2,Gm22,'g:*'); hold on;
plot(ts3,Gm31,'r:.'); hold on;
plot(ts3,Gm32,'r:*')
axis([0 0.004 0 10])
legend('\omega_n = 3400 rad/s no step delay','\omega_n = 3400 rad/s one-step delay', '\omega_n = 6400 rad/s no step delay','\omega_n = 6400 rad/s one-step delay','\omega_n = 9400 rad/s no step delay','\omega_n = 9400 rad/s one-step delay' )

Gm11_dB=20*log10(Gm11);
Gm12_dB=20*log10(Gm12);
Gm21_dB=20*log10(Gm21);
Gm22_dB=20*log10(Gm22);
Gm31_dB=20*log10(Gm31);
Gm32_dB=20*log10(Gm32);

Gm_dB_min=min([Gm11_dB Gm12_dB Gm21_dB Gm22_dB Gm31_dB Gm32_dB]);
Gm_dB_max=max([Gm11_dB Gm12_dB Gm21_dB Gm22_dB Gm31_dB Gm32_dB]);

figure(21)
plot(ts1,Gm11_dB,'b:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\omega_n = 3400 rad/s, 6400 rad/s, 9400 rad/s.')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts1,Gm12_dB,'b:*'); hold on;
plot(ts2,Gm21_dB,'g:.'); hold on;
plot(ts2,Gm22_dB,'g:*'); hold on;
plot(ts3,Gm31_dB,'r:.'); hold on;
plot(ts3,Gm32_dB,'r:*')
axis([0 0.004 Gm_dB_min Gm_dB_max])
% legend('\omega_n = 3400 rad/s no step delay','\omega_n = 3400 rad/s one-step delay', '\omega_n = 6400 rad/s no step delay','\omega_n = 6400 rad/s one-step delay','\omega_n = 9400 rad/s no step delay','\omega_n = 9400 rad/s one-step delay','Location', 'SouthEast' )

figure(22)
plot(ts1,Gm11_dB,'b:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\omega_n = 3400 rad/s')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts1,Gm12_dB,'b:*');
axis([0 0.004 Gm_dB_min Gm_dB_max])
legend('\omega_n = 3400 rad/s no step delay','\omega_n = 3400 rad/s one-step delay','Location', 'SouthEast')

figure(23)
plot(ts2,Gm21_dB,'g:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\omega_n = 6400 rad/s')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts2,Gm22_dB,'g:*'); 
axis([0 0.004 Gm_dB_min Gm_dB_max])
legend('\omega_n = 6400 rad/s no step delay','\omega_n = 6400 rad/s one-step delay','Location', 'SouthEast')

figure(24)
plot(ts3,Gm31_dB,'r:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\omega_n = 9400 rad/s.')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts3,Gm32_dB,'r:*')
axis([0 0.004 Gm_dB_min Gm_dB_max])
legend('\omega_n = 9400 rad/s no step delay','\omega_n = 9400 rad/s one-step delay','Location', 'SouthEast')
