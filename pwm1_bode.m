clear all; close all; clc

td=982.24*(10^(-6));
tt=[0.1:0.1:2]; %1x20
ts=tt*td;
 w=0:10:10^5;
 
for kk=1:length(w)
    fn(kk)=(1.688*exp(1j*w(kk)*ts(15))+1.161)/(exp(1j*w(kk)*ts(15))^2+1.375*exp(1j*w(kk)*ts(15))+0.473);
    mag(kk)=abs(fn(kk));
    phe(kk)=angle(fn(kk));
end

figure(1)
subplot(2,1,1)
semilogx(w,20*log10(mag));
grid on
ylabel('Magnitude (dB)')
title('Bode Diagram')
subplot(2,1,2)
semilogx(w,phe*(180/pi))
grid on
ylabel('Phase (deg))')
xlabel('Frequency (rad/s)')

