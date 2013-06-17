clear all
clc
close all
global nu;
nu=0:918;
c = 340;
k = 2*pi*nu/c;
%% Source Reflection Recorder Normal
a = sqrt(13)/2;
b = a;
theta_r = 2*pi - atan(2/3);
theta_s = atan(2/3);

data1 = Fresnel_Diffraction(a, b, theta_r, theta_s);

%% Source Reflection Recorder Normal
a = sqrt(29)/2;
b = sqrt(13)/2;
theta_r = 2*pi - atan(2/3);
theta_s = atan(2/5);

data2 = Fresnel_Diffraction(a, b, theta_r, theta_s);

%% Source Reflection Recorder Reflection
a = sqrt(29)/2;
b = a;
theta_r = 2*pi - atan(2/5);
theta_s = atan(2/5);

data3 = Fresnel_Diffraction(a, b, theta_r, theta_s);

%% Source normal Recorder Reflection
a = sqrt(13)/2;
b = sqrt(29)/2;
theta_r = 2*pi - atan(2/5);
theta_s = atan(2/3);

data4 = Fresnel_Diffraction(a, b, theta_r, theta_s);
plot(k,abs(data1-data2+data3-data4))
xlabel('Kd')
ylabel('Absolute Fresnel')
matlab2tikz('fresnel_diffraction_r1.tikz', 'height', '\figureheight', 'width', ...
            '\figurewidth','showInfo',false);