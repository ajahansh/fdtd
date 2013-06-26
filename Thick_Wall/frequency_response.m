%% CALCULATE TRANSFER FUNCTION
clear all;close all;clc
fp=[68,3*68,3^2*68];
k_max=2085;                       % Optimized for 500 Hz
m=3;                              % k_max is set for this m.
f_kmax=200;                       % Previous section
k=k_max*fp/f_kmax;                % Linear relationship with frequency
f=zeros(100*length(fp),1);
r1=zeros(100*length(fp),1);
r2=zeros(100*length(fp),1);
r3=zeros(100*length(fp),1);
fprintf(sprintf('Thick Wall Calculation...\n'))
for it=1:length(fp)               % Scan ricker wavelet
    wall=thick_wall(0,k(it),m,fp(it),0);
    f (1+(it-1)*100:it*100)=wall.f;
    r1(1+(it-1)*100:it*100)=wall.r1;
    r2(1+(it-1)*100:it*100)=wall.r2;
    r3(1+(it-1)*100:it*100)=wall.r3;
end
fprintf(sprintf('Free Space Calculation...\n'))
for it=1:length(fp)               % Divide transfer function by free wall
    free=thick_wall(0,k(it),m,fp(it),1);
    r1(1+(it-1)*100:it*100)=r1(1+(it-1)*100:it*100)./free.r1;
    r2(1+(it-1)*100:it*100)=r2(1+(it-1)*100:it*100)./free.r1;
    r3(1+(it-1)*100:it*100)=r3(1+(it-1)*100:it*100)./free.r1;
end
save('trans.mat','r1','r2','r3','f')
%% PLOT TRANSFER FUNCTION
load('trans.mat','r1','r2','r3','f')
kd=2*pi*f/340; %d=1m
semilogy(kd,r1,'r')
hold on
semilogy(kd,r2,'g')
semilogy(kd,r3,'b')
xlabel('Kd')
title('Frequency response of thin wall divided by free space')
legend('Recorder 1','Recorder 2','Recorder 3','Location', ...
       'NorthEast')
getframe();
saveas(gcf,'thick_wall','pdf')
matlab2tikz('transfer_thick_wall.tikz', 'height', '\figureheight', 'width', ...
            '\figurewidth','showInfo',false);