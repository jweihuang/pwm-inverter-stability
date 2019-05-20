clear all; close all; clc;

s=tf('s');
z=tf('z');

R=0.62;
L=1.22*(10^-3);
C=20*(10^-6);

wn=6400;
kxi=[0.0397 0.3397 0.1397];
td=0.001;
tt=[0.05:0.01:2.5];
ts=tt.*td(1);

wid=0.1;
gain=0:wid:10;

for hh=1:length(ts)
    
    opn_fns1{1,hh}=(wn^2)/(s^2+2*kxi(1)*wn*s+wn^2);
    opn_fnz1{1,hh}=c2d(opn_fns1{1,hh},ts(hh),'zoh');
    opn_fnz_stpdy1{1,hh}=c2d(opn_fns1{1,hh},ts(hh),'zoh')*(1/z);
    clo_fnz1{1,hh}=feedback(opn_fnz1{1,hh},1);
    clo_fnz_stpdy1{1,hh}=feedback(opn_fnz_stpdy1{1,hh},1);
    
    r11=rlocus(opn_fnz1{1,hh},gain);
    rho11=abs(r11)-1;
    [m,n]=size(rho11);
    zero_index111=find(rho11(1,1:(end-1)).*rho11(1,(2:end))<0);
    zero_index112=find(rho11(2,1:(end-1)).*rho11(2,(2:end))<0);
    zero_index11=min([zero_index111 zero_index112]);
    stb_gain11(1,hh)=zero_index11(1,1)*wid;
    opn_fnw1{1,hh}=d2c(opn_fnz1{1,hh}, 'tustin'); %bilinear transform w-domain
    opn_fnw_stpdy1{1,hh}=d2c(opn_fnz_stpdy1{1,hh}, 'tustin');
    
    r12=rlocus(opn_fnz_stpdy1{1,hh},gain);
    rho12=(abs(r12)-1);
    [m,n]=size(rho12);
    zero_index121=find(rho12(1,1:(end-1)).*rho12(1,(2:end))<0);
    zero_index122=find(rho12(2,1:(end-1)).*rho12(2,(2:end))<0);
    zero_index12=min([zero_index121 zero_index122]);
    stb_gain12(1,hh)=zero_index12(1,1)*wid;
    
    S11{1,hh}=allmargin(opn_fnz1{1,hh});
    S12{1,hh}=allmargin(opn_fnz_stpdy1{1,hh});
    Gm11(hh)=S11{1,hh}.GainMargin(1,1); %gain margin
    Gm12(hh)=S12{1,hh}.GainMargin(1,1);
%     Pm11(hh)=max(S11{1,hh}.PhaseMargin); %phase margin
%     Pm12(hh)=max(S12{1,hh}.PhaseMargin);
%     Wg11(hh)=max(S11{1,hh}.PMFrequency); %gain crossover freqency
%     Wg12(hh)=max(S12{1,hh}.PMFrequency);
    Wp11(hh)=S11{1,hh}.GMFrequency(1,1); %phase crossover frequency
    Wp12(hh)=S12{1,hh}.GMFrequency(1,1);
    clear r11 r12 rho11 rho12 zero_index11 zero_index_12
end

for ff=1:length(ts)
    
    opn_fns2{1,ff}=(wn^2)/(s^2+2*kxi(2)*wn*s+wn^2);
    opn_fnz2{1,ff}=c2d(opn_fns2{1,ff},ts(ff),'zoh');
    opn_fnz_stpdy2{1,ff}=c2d(opn_fns2{1,ff},ts(ff),'zoh')*(1/z);
    clo_fnz2{1,ff}=feedback(opn_fnz2{1,ff},1);
    clo_fnz_stpdy2{1,ff}=feedback(opn_fnz_stpdy2{1,ff},1);
    
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
    
    S21{1,ff}=allmargin(opn_fnz2{1,ff});
    S22{1,ff}=allmargin(opn_fnz_stpdy2{1,ff});
    Gm21(ff)=S21{1,ff}.GainMargin(1,1);
    Gm22(ff)=S22{1,ff}.GainMargin(1,1);
    Wp21(ff)=S21{1,ff}.GMFrequency(1,1);
    Wp22(ff)=S22{1,ff}.GMFrequency(1,1);
    clear r21 r22 rho21 rho22 zero_index21 zero_index22
end


for ss=1:length(ts)
    
    opn_fns3{1,ss}=(wn^2)/(s^2+2*kxi(3)*wn*s+wn^2);
    opn_fnz3{1,ss}=c2d(opn_fns3{1,ss},ts(ss),'zoh');
    opn_fnz_stpdy3{1,ss}=c2d(opn_fns3{1,ss},ts(ss),'zoh')*(1/z);
    clo_fnz3{1,ss}=feedback(opn_fnz3{1,ss},1);
    clo_fnz_stpdy3{1,ss}=feedback(opn_fnz_stpdy3{1,ss},1);
    
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
    
    S31{1,ss}=allmargin(opn_fnz3{1,ss});
    S32{1,ss}=allmargin(opn_fnz_stpdy3{1,ss});
    Gm31(ss)=S31{1,ss}.GainMargin(1,1);
    Gm32(ss)=S32{1,ss}.GainMargin(1,1);
    Wp31(ss)=S31{1,ss}.GMFrequency(1,1);
    Wp32(ss)=S32{1,ss}.GMFrequency(1,1);
    clear r31 r32 rho31 rho32 zero_index31 zero_index32
end

figure(10)

plot(ts,stb_gain11,'b:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\xi = 0.0397, 0.3397, 0.1397')
xlabel('T_s(sec)'); ylabel('K');
plot(ts,stb_gain12,'b:x'); hold on;
plot(ts,stb_gain21,'g:.'); hold on;
plot(ts,stb_gain22,'g:x'); hold on;
plot(ts,stb_gain31,'r:.'); hold on;
plot(ts,stb_gain32,'r:x'); 
axis([0 0.0025 0 10])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay', '\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay','\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay')

figure(11)

plot(ts,stb_gain11,'b:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\xi = 0.0397')
xlabel('T_s(sec)'); ylabel('K');
plot(ts,stb_gain12,'b:x'); hold on;
axis([0 0.0025 0 10])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay')

figure(12)
plot(ts,stb_gain21,'g:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\xi = 0.3397')
xlabel('T_s(sec)'); ylabel('K');
plot(ts,stb_gain22,'g:x'); hold on;
axis([0 0.0025 0 10])
legend('\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay')

figure(13)
plot(ts,stb_gain31,'r:.'); hold on;
title('comparisons stability of the one-step-delay with no delay.\xi = 0.1397')
xlabel('T_s(sec)'); ylabel('K');
plot(ts,stb_gain32,'r:x'); 
axis([0 0.0025 0 10])
legend('\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay')

figure(21)
plot(ts,Gm11,'b:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\xi = 0.0397, 0.3397, 0.1397')
xlabel('T_s(sec)'); ylabel('absolute magnitude');
plot(ts,Gm12,'b:x'); hold on;
plot(ts,Gm21,'g:.'); hold on;
plot(ts,Gm22,'g:x'); hold on;
plot(ts,Gm31,'r:.'); hold on;
plot(ts,Gm32,'r:x')
axis([0 0.0025 0 10])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay', '\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay','\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay')

Gm11_dB=20*log10(Gm11);
Gm12_dB=20*log10(Gm12);
Gm21_dB=20*log10(Gm21);
Gm22_dB=20*log10(Gm22);
Gm31_dB=20*log10(Gm31);
Gm32_dB=20*log10(Gm32);

Gm_dB_min=min([Gm11_dB Gm12_dB Gm21_dB Gm22_dB Gm31_dB Gm32_dB]);
Gm_dB_max=max([Gm11_dB Gm12_dB Gm21_dB Gm22_dB Gm31_dB Gm32_dB]);

figure(22)
plot(ts,Gm11_dB,'b:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\xi = 0.0397, 0.3397, 0.1397')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts,Gm12_dB,'b:x'); hold on;
plot(ts,Gm21_dB,'g:.'); hold on;
plot(ts,Gm22_dB,'g:x'); hold on;
plot(ts,Gm31_dB,'r:.'); hold on;
plot(ts,Gm32_dB,'r:x')
axis([0 0.0025 Gm_dB_min Gm_dB_max])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay', '\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay','\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay','Location', 'SouthEast')

figure(23)
plot(ts,Gm11_dB,'b:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\xi = 0.0397')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts,Gm12_dB,'b:x'); 
axis([0 0.0025 Gm_dB_min Gm_dB_max])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay','Location', 'SouthEast')

figure(24)
plot(ts,Gm21_dB,'g:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\xi = 0.3397')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts,Gm22_dB,'g:x'); 
axis([0 0.0025 Gm_dB_min Gm_dB_max])
legend('\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay','Location', 'SouthEast')

figure(25)
plot(ts,Gm31_dB,'r:.'); hold on;
title('Gain margin of the one-step-delay with no delay.\xi = 0.1397')
xlabel('T_s(sec)'); ylabel('magnitude(dB)');
plot(ts,Gm32_dB,'r:x')
axis([0 0.0025 Gm_dB_min Gm_dB_max])
legend('\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay','Location', 'SouthEast')

pm_min=min([Pm11 Pm12 Pm21 Pm22 Pm31 Pm32]);
pm_max=max([Pm11 Pm12 Pm21 Pm22 Pm31 Pm32]);

figure(3)
plot(ts,Pm11,'b.'); hold on;
title('Phase margin of the one-step-delay with no delay.\xi = 0.0397, 0.3397, 0.1397')
xlabel('T_s(sec)'); ylabel('Phase(degree)');
plot(ts,Pm12,'bx'); hold on;
plot(ts,Pm21,'g.'); hold on;
plot(ts,Pm22,'gx'); hold on;
plot(ts,Pm31,'r.'); hold on;
plot(ts,Pm32,'rx')
axis([0 0.0025 pm_min pm_max])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay', '\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay','\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay')

Wg_min=min([Wg11 Wg12 Wg21 Wg22 Wg31 Wg32]);
Wg_max=max([Wg11 Wg12 Wg21 Wg22 Wg31 Wg32]);

figure(4)
plot(ts,Wg11,'b.'); hold on;
title('Gain crossover frequency of the one-step-delay with no delay.\xi = 0.0397, 0.3397, 0.1397')
xlabel('T_s(sec)'); ylabel('rad/s');
plot(ts,Wg12,'bx'); hold on;
plot(ts,Wg21,'g.'); hold on;
plot(ts,Wg22,'gx'); hold on;
plot(ts,Wg31,'r.'); hold on;
plot(ts,Wg32,'rx')
axis([0 0.0025 Wg_min Wg_max])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay', '\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay','\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay')

Wp_min=min([Wp11 Wp12 Wp21 Wp22 Wp31 Wp32]);
Wp_max=max([Wp11 Wp12 Wp21 Wp22 Wp31 Wp32]);

figure(5)
plot(ts,Wp11,'b.'); hold on;
title('Phase crossover frequency of the one-step-delay with no delay.\xi = 0.0397, 0.3397, 0.1397')
xlabel('T_s(sec)'); ylabel('rad/s');
plot(ts,Wp12,'bx'); hold on;
plot(ts,Wp21,'g.'); hold on;
plot(ts,Wp22,'gx'); hold on;
plot(ts,Wp31,'r.'); hold on;
plot(ts,Wp32,'rx')
axis([0 0.0025 Wp_min Wp_max])
legend('\xi = 0.0397 no step delay','\xi = 0.0397 one-step delay', '\xi = 0.3397 no step delay','\xi = 0.3397 one-step delay','\xi = 0.1397 no step delay','\xi = 0.1397 one-step delay')

figure(81)
margin(opn_fnz1{1,29});
figure(82)
margin(opn_fnz1{1,30});
figure(83)
margin(opn_fnz1{1,64});
figure(84)
margin(opn_fnz1{1,65});

figure(71)
margin(opn_fnw1{1,29});
figure(72)
margin(opn_fnw1{1,30});
figure(73)
margin(opn_fnw1{1,64});
figure(74)
margin(opn_fnw1{1,65});


figure(91)
margin(opn_fnz1{1,124});
figure(92)
margin(opn_fnz_stpdy1{1,126});
figure(93)
margin(opn_fnz1{1,128});
figure(94)
margin(opn_fnz_stpdy1{1,130});