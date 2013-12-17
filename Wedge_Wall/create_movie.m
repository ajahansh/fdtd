%% CALCULATE TRANSFER FUNCTION
clear all;close all;clc
fp=[68,3*68,3^2*68];
k_max=2085;                       % Previous section at 500 Hz with wall.
m=3;                              % k_max is set for this m.
f_kmax=200;                       % Previous section
k=k_max*fp/f_kmax;                % Linear relationship with frequency
N=[3,5,7,9];                              % Number of discretization in fine grid
fprintf(sprintf('Creating Movies...\n'))
for it_n=1:length(N)
    for it=1:length(fp)
        fprintf(sprintf('f=%d, N=%d\n\n',fp(it),N(it_n)))
        wedge_fine(0,k(it),m,fp(it),N(it_n),1);
    end
end

