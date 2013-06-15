function data = Fresnel_Diffraction(a, b, theta_r, theta_s)

global nu;
c = 340;
k = 2*pi*nu/c;

L = a*b/(a+b);

alpha_p = theta_r + theta_s;
alpha_m = theta_r - theta_s;

n = 2;

N_p = 1;
N_m = 0;

A_p_alpha_p = 2*cos((2*n*pi*N_p - alpha_p)/2)^2; 
A_p_alpha_m = 2*cos((2*n*pi*N_p - alpha_m)/2)^2; 
A_m_alpha_p = 2*cos((2*n*pi*N_m - alpha_p)/2)^2; 
A_m_alpha_m = 2*cos((2*n*pi*N_m - alpha_m)/2)^2; 

Constant = exp((a+b)*k*1j + 1j*pi/4)/sqrt(L)./sqrt(2*pi*k)/2/n;

data = Constant.*(cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_m) + ...
                  cot((pi - alpha_p)/2/n)*fcs(k*L*A_m_alpha_p) + ...
                  cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_p) + ...
                  cot((pi + alpha_m)/2/n)*fcs(k*L*A_p_alpha_m));


end