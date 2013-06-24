%% FIND BEST KMAX FOR 500 HZ
% $$$ clear all;close all;clc
% $$$ delete('frequency-response.txt')
% $$$ diary('frequency-response.txt')
% $$$ disp('Finding best K_max in 500 Hz to have least reflection')
% $$$ tic
% $$$ f_kmax=500;                      % Kmax is found for this frequency
% $$$ m=3;
% $$$ options=optimset('MaxFunEvals',30,'GradObj','on','Display','iter');
% $$$ % In this part m=3. Initial value for k_max is 20e3.
% $$$ [k_max,return_args]=fminsearch(@(k)SIT_SIP_PML(0,k,m,f_kmax,1),...
% $$$                                20e3,options);
% $$$ disp(['Optimal k_max=',num2str(k_max),' at ', num2str(f_kmax), ' Hz'])
% $$$ time=toc;
% $$$ disp(['Optimization for k_max took: ',num2str(time), ' s'])
% $$$ disp(['Remaining Energy=',num2str(return_args)])
% $$$ diary off
%% CALCULATE TRANSFER FUNCTION
clear all;close all;clc
fp=[68,3*68,3^2*68,3^3*68,3^4*68];
k_max=5637;                       % Previous section at 500 Hz with wall.
m=3;                              % k_max is set for this m.
f_kmax=500;                       % Previous section
k=k_max*fp/f_kmax;                % Linear relationship with frequency
f=zeros(100*length(fp),1);
r1=zeros(100*length(fp),1);
r2=zeros(100*length(fp),1);
r3=zeros(100*length(fp),1);
fprintf(sprintf('Thin Wall Calculation...\n'))
for it=1:length(fp)               % Scan ricker wavelet
    wall=thin_wall(0,k(it),m,fp(it),0);
    f (1+(it-1)*100:it*100)=wall.f;
    r1(1+(it-1)*100:it*100)=wall.r1;
    r2(1+(it-1)*100:it*100)=wall.r2;
    r3(1+(it-1)*100:it*100)=wall.r3;
end
fprintf(sprintf('Free Space Calculation...\n'))
for it=1:length(fp)               % Divide transfer function by free wall
    free=thin_wall(0,k(it),m,fp(it),1);
    r1(1+(it-1)*100:it*100)=r1(1+(it-1)*100:it*100)./free.r1;
    r2(1+(it-1)*100:it*100)=r2(1+(it-1)*100:it*100)./free.r1;
    r3(1+(it-1)*100:it*100)=r3(1+(it-1)*100:it*100)./free.r1;
end
save('trans.mat','r1','r2','r3','f')
%% PLOT TRANSFER FUNCTION
load('trans.mat','r1','r2','r3','f')
kd=2*pi*f/340; %d=1m
figure 
hold on
plot(kd,r1,'r')
plot(kd,r2,'g')
plot(kd,r3,'b')
xlabel('Kd')
title('Frequency response of thin wall divided by free space')
legend('Recorder 1','Recorder 2','Recorder 3','Location', ...
       'NorthEast')
getframe();
saveas(gcf,'thin_wall','pdf')
%matlab2tikz('transfer_thin_wall.tikz', 'height', '\figureheight', 'width', ...
%            '\figurewidth','showInfo',false);