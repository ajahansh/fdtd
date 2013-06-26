%% CALCULATE TRANSFER FUNCTION
clear all;close all;clc
fp=[68,3*68,3^2*68];
k_max=2085;                       % Previous section at 500 Hz with wall.
m=3;                              % k_max is set for this m.
f_kmax=200;                       % Previous section
k=k_max*fp/f_kmax;                % Linear relationship with frequency
N=[3,5,7,9];                        % Number of discretization in fine grid
f=zeros(1,length(fp)*100);
r1=zeros(length(fp)*100,length(N));
r2=r1;
r3=r2;
fprintf(sprintf('Wedge Calculation...\n'))
for it_n=1:length(N)
    fprintf(sprintf('N=%d\n',N(it_n)))
    for it=1:length(fp)               % Scan ricker wavelet
        fine=wedge_fine(0,k(it),m,fp(it),N(it_n));
        r1(1+(it-1)*100:it*100,it_n)=fine.r1;
        r2(1+(it-1)*100:it*100,it_n)=fine.r2;
        r3(1+(it-1)*100:it*100,it_n)=fine.r3;
    end
end
fprintf(sprintf('Free Space Calculation...\n'))
field='free';
for it=1:length(fp)               % Run in free space and divide by it
    free=free_space(0,k(it),m,fp(it),1);
    f(1+(it-1)*100:it*100)=free.f;
    for it_n=1:length(N)
        r1(1+(it-1)*100:it*100,it_n)=r1(1+(it-1)*100:it*100,it_n)./free.r1;
        r2(1+(it-1)*100:it*100,it_n)=r2(1+(it-1)*100:it*100,it_n)./free.r1;
        r3(1+(it-1)*100:it*100,it_n)=r3(1+(it-1)*100:it*100,it_n)./free.r1;
    end
end
save('trans.mat','r1','r2','r3','f')
%% PLOT TRANSFER FUNCTION WITH RESPECT TO N
load('trans.mat','r1','r2','r3','f')
kd=2*pi*f/340; %d=1m
h1=semilogy(kd,r1(:,1),'r');
hold on
h2=semilogy(kd,r1(:,2),'g');
h3=semilogy(kd,r1(:,3),'b');
h4=semilogy(kd,r1(:,4),'k');
xlabel('Kd')
legend([h1,h2,h3,h4],'N=3','N=5','N=7','N=9')
matlab2tikz('r1.tikz', 'height', '\figureheight', 'width', ...
            '\figurewidth','showInfo',false);
export_fig('r1','-pdf','-transparent')
hold off
h1=semilogy(kd,r2(:,1),'r');
hold on
h2=semilogy(kd,r2(:,2),'g');
h3=semilogy(kd,r2(:,3),'b');
h4=semilogy(kd,r2(:,4),'k');
xlabel('Kd')
legend([h1,h2,h3,h4],'N=3','N=5','N=7','N=9')
matlab2tikz('r2.tikz', 'height', '\figureheight', 'width', ...
            '\figurewidth','showInfo',false);
export_fig('r2','-pdf','-transparent')
hold off
h1=semilogy(kd,r3(:,1),'r');
hold on
h2=semilogy(kd,r3(:,2),'g');
h3=semilogy(kd,r3(:,3),'b');
h4=semilogy(kd,r3(:,4),'k');
xlabel('Kd')
legend([h1,h2,h3,h4],'N=3','N=5','N=7','N=9')
matlab2tikz('r3.tikz', 'height', '\figureheight', 'width', ...
            '\figurewidth','showInfo',false);
export_fig('r3','-pdf','-transparent')