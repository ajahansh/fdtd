clear all
clc
close all
global nu;
nu=linspace(1,2700,1000);
c=340;
k=2*pi*nu/c;
%% Source Normal Recorder Normal - Recorder 1
a=sqrt(13)/2;
b=a;
theta_r=2*pi-atan(2/3);
theta_s=atan(2/3);

data1=Fresnel_Diffraction(a, b, theta_r, theta_s);
%% Source Reflection Recorder Normal - Recorder 1
a=sqrt(29)/2;
b=sqrt(13)/2;
theta_r=2*pi-atan(2/3);
theta_s=atan(2/5);

data2=Fresnel_Diffraction(a, b, theta_r, theta_s);
%% Source Reflection Recorder Reflection - Recorder 1
a=sqrt(29)/2;
b=a;
theta_r=2*pi-atan(2/5);
theta_s=atan(2/5);

data3=Fresnel_Diffraction(a,b,theta_r,theta_s);
%% Source normal Recorder Reflection - Recorder 1
a=sqrt(13)/2;
b=sqrt(29)/2;
theta_r=2*pi-atan(2/5);
theta_s=atan(2/3);

data4=Fresnel_Diffraction(a, b, theta_r, theta_s);
semilogy(k,abs(data1-data2+data3-data4),'r')
hold on
%% Source Normal Recorder Normal - Recorder 2
a=sqrt(13)/2;
b=sqrt(25)/2;
theta_r=2*pi-atan(4/3);
theta_s=atan(2/3);

data1=Fresnel_Diffraction(a, b, theta_r, theta_s);
%% Source Reflection Recorder Normal - Recorder 2
a=sqrt(29)/2;
b=sqrt(25)/2;
theta_r=2*pi-atan(4/3);
theta_s=atan(2/5);

data2=Fresnel_Diffraction(a, b, theta_r, theta_s);
%% Source Reflection Recorder Reflection - Recorder 2
a=sqrt(29)/2;
b=sqrt(41)/2;
theta_r=2*pi-atan(4/5);
theta_s=atan(2/5);

data3=Fresnel_Diffraction(a,b,theta_r,theta_s);
%% Source normal Recorder Reflection - Recorder 2
a=sqrt(13)/2;
b=sqrt(41)/2;
theta_r=2*pi-atan(4/5);
theta_s=atan(2/3);

data4=Fresnel_Diffraction(a, b, theta_r, theta_s);
semilogy(k,abs(data1-data2+data3-data4),'g')
%% Source Normal Recorder Normal - Recorder 3
a=sqrt(13)/2;
b=sqrt(45)/2;
theta_r=2*pi-atan(2);
theta_s=atan(2/3);

data1=Fresnel_Diffraction(a, b, theta_r, theta_s);
%% Source Reflection Recorder Normal - Recorder 3
a=sqrt(29)/2;
b=sqrt(45)/2;
theta_r=2*pi-atan(2);
theta_s=atan(2/5);

data2=Fresnel_Diffraction(a, b, theta_r, theta_s);
%% Source Reflection Recorder Reflection - Recorder 3
a=sqrt(29)/2;
b=sqrt(61)/2;
theta_r=2*pi-atan(6/5);
theta_s=atan(2/5);

data3=Fresnel_Diffraction(a,b,theta_r,theta_s);
%% Source normal Recorder Reflection - Recorder 3
a=sqrt(13)/2;
b=sqrt(61)/2;
theta_r=2*pi-atan(6/5);
theta_s=atan(2/3);

data4=Fresnel_Diffraction(a, b, theta_r, theta_s);
semilogy(k,abs(data1-data2+data3-data4),'b')

xlabel('Kd')
ylabel('Absolute Fresnel')
legend('Recorder 1', 'Recorder 2', 'Recorder 3')
saveas(gcf,'fresnel_diffraction','pdf')
matlab2tikz('fresnel_diffraction.tikz', 'height', '\figureheight', 'width', ...
            '\figurewidth','showInfo',false);
