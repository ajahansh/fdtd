d = 1;
a = sqrt(13)/2*d;
b = sqrt(13)/2*d;

nu = 0:1000;
c = 340;
k = 2*pi.*nu./c;

L = a*b/(a+b);

theta_r = 2*pi - atan(1/3);
theta_s = atan(1/3);
alpha_p = theta_r + theta_s;
alpha_m = theta_r - theta_s;

n = 1;

N_p = 1;
N_m = 0;

A_p_alpha_p = 2*cos((2*2*pi*N_p - alpha_p)/2)^2; 
A_p_alpha_m = 2*cos((2*2*pi*N_p - alpha_m)/2)^2; 
A_m_alpha_p = 2*cos((2*2*pi*N_m - alpha_p)/2)^2; 
A_m_alpha_m = 2*cos((2*2*pi*N_m - alpha_m)/2)^2; 


Constant = exp((a+b)*k*1j + 1j*3*pi/4)/sqrt(L)./sqrt(2*pi*k)/2/n;

Field = Constant.*(cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_m) + ...
                  cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_p) + ...
                  cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_p) + ...
                  cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_m));

 %% Source Reflection Recorder Normal
d = 1;
a = sqrt(29)/2*d;
b = sqrt(13)/2*d;

nu = 0:1000;
c = 340;
k = 2*pi.*nu./c;

L = a*b/(a+b);

theta_r = 2*pi - atan(1/3);
theta_s = atan(2/5);
alpha_p = theta_r + theta_s;
alpha_m = theta_r - theta_s;

n = 1;

N_p = 1;
N_m = 0;

A_p_alpha_p = 2*cos((2*2*pi*N_p - alpha_p)/2)^2; 
A_p_alpha_m = 2*cos((2*2*pi*N_p - alpha_m)/2)^2; 
A_m_alpha_p = 2*cos((2*2*pi*N_m - alpha_p)/2)^2; 
A_m_alpha_m = 2*cos((2*2*pi*N_m - alpha_m)/2)^2; 


Constant = exp((a+b)*k*1j + 1j*3*pi/4)/sqrt(L)./sqrt(2*pi*k)/2/n;

Field = Constant.*(cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_m) + ...
                   cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_p) + ...
                   cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_p) + ...
                   cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_m)) + Field;

  %% Source Normal Recorder Reflection
d = 1;
a = sqrt(13)/2*d;
b = sqrt(29)/2*d;

nu = 0:1000;
c = 340;
k = 2*pi.*nu./c;

L = a*b/(a+b);

theta_r = 2*pi - atan(2/5);
theta_s = atan(1/3);
alpha_p = theta_r + theta_s;
alpha_m = theta_r - theta_s;

n = 1;

N_p = 1;
N_m = 0;

A_p_alpha_p = 2*cos((2*2*pi*N_p - alpha_p)/2)^2; 
A_p_alpha_m = 2*cos((2*2*pi*N_p - alpha_m)/2)^2; 
A_m_alpha_p = 2*cos((2*2*pi*N_m - alpha_p)/2)^2; 
A_m_alpha_m = 2*cos((2*2*pi*N_m - alpha_m)/2)^2; 


Constant = exp((a+b)*k*1j + 1j*3*pi/4)/sqrt(L)./sqrt(2*pi*k)/2/n;

Field = Constant.*(cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_m) + ...
                   cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_p) + ...
                   cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_p) + ...
                   cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_m)) + Field;
 %% Source Reflection Recorder Reflection
d = 1;
a = sqrt(29)/2*d;
b = sqrt(29)/2*d;

nu = 0:1000;
c = 340;
k = 2*pi.*nu./c;

L = a*b/(a+b);

theta_r = 2*pi - atan(2/5);
theta_s = atan(2/5);
alpha_p = theta_r + theta_s;
alpha_m = theta_r - theta_s;

n = 1;

N_p = 1;
N_m = 0;

A_p_alpha_p = 2*cos((2*2*pi*N_p - alpha_p)/2)^2; 
A_p_alpha_m = 2*cos((2*2*pi*N_p - alpha_m)/2)^2; 
A_m_alpha_p = 2*cos((2*2*pi*N_m - alpha_p)/2)^2; 
A_m_alpha_m = 2*cos((2*2*pi*N_m - alpha_m)/2)^2; 


Constant = exp((a+b)*k*1j + 1j*3*pi/4)/sqrt(L)./sqrt(2*pi*k)/2/n;

Field = Constant.*(cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_m) + ...
                   cot((pi - alpha_m)/2/n)*fcs(k*L*A_m_alpha_p) + ...
                   cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_p) + ...
                   cot((pi + alpha_p)/2/n)*fcs(k*L*A_p_alpha_m)) + Field
               
               
plot(2*pi*nu/c, abs(Field));              