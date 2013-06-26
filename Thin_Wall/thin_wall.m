function return_args=thin_wall(ifshow,k_max,m,fp,if_free)
if nargin<1
    close all
    clc
    ifshow=1;
    fp=68*9;
    k_max=2085*fp/200;
    m=3;
    if_free=1;
end
%% INITIALIZATION OF GRID AND SIMULATING PARAMETERS
c=340;          % Speed of sound
lp=c/fp;        % Associated wavelength with peak frequency
md=2;           % dr=md/fp
d=1.0;          % d distance in exercise (in m)
dx=c/fp/20;     % Spatial discretisation step, at least 20 samples
                % per wavelength
nd=ceil(d/dx);  % number of grids in d.
if mod(nd,2)==1 % nd should be even, because nd/2 is needed.
    nd=nd+1;
end
dx=d/nd;             % Recalculate dx
dy=dx;               % Same grid in x,y
np=lp/dx;            % Does not need to be integer
nx=6*nd+2+1;         % Number of cells in x direction
ny=6*nd+1;           % Number of cells in y direction
npml=40;             % Number of cells in PML
xdim=nx-2+2*npml;    % Total number of columns
ydim=ny-2+2*npml;    % Total number of rows
cn=1/sqrt(2);
dt=cn*dx/c/sqrt(2);  % Time step
cdtdx=c*dt/dx;       % c * dt / dx. Needed for ricker wavelet
nt=round(2*md/(fp*dt)+...
         sqrt(xdim^2+ydim^2)/cdtdx);%Number of time steps
%%% FFT INITIALIZATIONS
fmax=1/(2*dt);       % Maximum frequency in fft
NFFT=2^nextpow2(nt);
disp(['fp=',num2str(fp),' d=',num2str(d),' m'])
disp(['nx=',num2str(nx),' dx=',num2str(dx),' nt=',num2str(nt)])
%% OTHER INITIALIZATIONS
col_s=npml+nd;       % X location of source
row_s=npml+ny-1-nd/2;% Y location of source
col_r1=col_s+2*nd+2; % X location of recorder1
col_r2=col_r1+nd;    % X location of recorder2
col_r3=col_r2+nd;    % X location of recorder3
source=arrayfun(@(it)ricker_wavelet(cdtdx,np,it,md),1:nt);
r1=zeros(nt,1);      % Recorder1 Pre alocation
r2=r1;               % Recorder2 Pre alocation
r3=r2;               % Recorder3 Pre alocation
fft_source=fft(source,NFFT)/nt; 
f=fmax*linspace(0,1,NFFT/2+1)';
f=f(f<2.5*fp);
fft_source=abs(fft_source(1:length(f)))';
if ifshow
    figure(1)% to show the fft of source
    plot(f,2*fft_source,'bo')
    xlabel('Frequency(Hz)')
    title('fft of source and recorders')
end
% Defining the fields
ox=zeros(ydim,xdim);
oy=ox;
p=ox;
px=spalloc(ydim,xdim,2*npml*ydim+(nx-2)*2*npml);
py=px; % px and py are sparse. They fit around p.
% Frame is just made to better see the boundaries and source
% Frame is prealocated for more elements that actually needed.
frame=sparse(ydim,xdim,3*ydim+3*xdim);
frame(npml,:)=1;frame(ny+npml-1,:)=1;
frame(:,npml)=1;frame(:,nx+npml-1)=1;
frame(npml+ny-2*nd-1:npml+ny-1,npml+2*nd+1)=1;
%% INITIALIZING THE DAMPING MATRICES
% Starting out with the ox and oy...
KO_x=spalloc(ydim,xdim,(2*npml-1)*ydim);
KO_y=spalloc(ydim,xdim,npml*xdim);
dpml=npml-0.5;%maximum length of PML
dist=(0.5:1:dpml);%All Os have 0.5 offset.
base_damp=repmat(k_max*(dist'/dpml).^m,1,xdim);
KO_y(1:npml,:)=flipud(base_damp);%Top PML - OY
KO_y(end-npml+2:end,:)=base_damp(1:end-1,:); % Bottom PML - OY
base_damp=repmat(k_max*(dist/npml).^m,ydim,1);%needed for Ox
KO_x(:,end-npml+1:end)=base_damp;%right part of Ox
% The left part of Ox has 1 element less
KO_x(:,1:npml-1)=fliplr(base_damp(:,1:end-1));
% Now calculating the damping matrices for p.
KP_x=spalloc(ydim,xdim,2*npml*ydim);
KP_y=spalloc(ydim,xdim,2*npml*ydim);
dist=(0:npml-1);% Here there is no offset.
base_damp=repmat(k_max*(dist'/dpml).^m,1,xdim);
KP_y(end-npml+1:end,:)=base_damp;
KP_y(1:npml,:)=flipud(base_damp);
base_damp=repmat(k_max*(dist/dpml).^m,ydim,1);%needed for px
KP_x(:,1:npml)=fliplr(base_damp);
KP_x(:,end-npml+1:end)=base_damp;
% Construction of some needed constant matrices in update equations.
KO_x_o=(1-KO_x*dt/2)./(1+KO_x*dt/2);%Multiplied by previous ox
KO_x_p=dt/dx./(1+KO_x*dt/2);%Multiplied by p.
KO_y_o=(1-KO_y*dt/2)./(1+KO_y*dt/2);%Multiplied by previous oy
KO_y_p=dt/dy./(1+KO_y*dt/2);%Multiplied by p.
KP_x_o=c^2*dt/dx./(1+KP_x*dt/2);
KP_x_p=(1-KP_x*dt/2)./(1+KP_x*dt/2);
KP_y_o=c^2*dt/dy./(1+KP_y*dt/2);
KP_y_p=(1-KP_y*dt/2)./(1+KP_y*dt/2);
%% TIME STEPPING ITERATION
if ifshow;
    figure(2);
end
display='it=0/nt t=0.000000';
fprintf(display)
for it=1:nt
    t=(it-1)*dt;
    fprintf(repmat('\b',1,length(display)))
    display=sprintf('it=%d/%d t=%.6fs',it,nt,t);
    fprintf(display);
    p(row_s,col_s)=p(row_s,col_s)+source(it);
    r1(it)=p(row_s,col_r1);
    r2(it)=p(row_s,col_r2);
    r3(it)=p(row_s,col_r3);
    
    % Update Equations
    ox=fin_diff([p+px+py,zeros(ydim,1)],ox,'x',KO_x_p,KO_x_o);
    oy=fin_diff([zeros(1,xdim);p+px+py],oy,'y',KO_y_p,KO_y_o);
    p=fin_diff([ox(:,1),ox],px+p,'x',KP_x_o,KP_x_p);

    px(:,[1:npml,nx+npml-1:end])=p(:,[1:npml,nx+npml-1:end]);
    px([1:npml,ny+npml-1:end],:)=p([1:npml,ny+npml-1:end],:);
    p(:,[1:npml,nx+npml-1:end])=0;
    p([1:npml,ny+npml-1:end],:)=0;
    
    p=fin_diff([oy;zeros(1,xdim)],py+p,'y',KP_y_o,KP_y_p);
    
    py(:,[1:npml,nx+npml-1:end])=p(:,[1:npml,nx+npml-1:end]); %#ok<*SPRIX>
    py([1:npml,ny+npml-1:end],:)=p([1:npml,ny+npml-1:end],:); %#ok<SPRIX>
    p(:,[1:npml,nx+npml-1:end])=0;
    p([1:npml,ny+npml-1:end],:)=0;
    
    % Boundries
    if ~if_free
        p (npml+ny-2*nd-1:npml+ny-1,npml+2*nd+1)=0; % Thin wall
        p (npml+ny-1,:)=0;  % Ground
        px(npml+ny-1,:)=0;
        py(npml+ny-1,:)=0;
    end
  
    if ifshow && mod(it,20)==0;
        pcolor(flipud(p+px+py+frame));
        colormap('gray')
        shading interp
        caxis([-.01 .01])
        colorbar
        xlim([1 xdim])
        ylim([1 ydim])
        title(display)
        getframe(); 
    end
end
fprintf('\n')
%% POST PROCESSING
return_args=struct('f',0,'r1',0,'r2',0,'r3',0,'E',0);

% Calculate remaining energy
ox=ox(npml:end-npml+1,npml:end-npml+1);
oy=oy(npml:end-npml+1,npml:end-npml+1);
energy=sum(sum(ox.^2+oy.^2));
return_args.E=energy;

% Calculate fft
fft_r1=fft(r1,NFFT)/nt;
fft_r2=fft(r2,NFFT)/nt;
fft_r3=fft(r3,NFFT)/nt;
fft_r1=abs(fft_r1(1:length(f)));
fft_r2=abs(fft_r2(1:length(f)));
fft_r3=abs(fft_r3(1:length(f)));
% Plot FFT of recorders if requested
if ifshow
    figure(1);%Back to fft plot
    hold on
    h1=plot(f,2*fft_r1,'ro');
    h2=plot(f,2*fft_r2,'go');
    h3=plot(f,2*fft_r3,'co');
    legend([h1,h2,h3],'Recorder 1','Recorder 2','Recorder 3','Location',...
           'NorthEast')
end
% Fit curve using spline
ft=fittype('splineinterp');
opts=fitoptions(ft);
opts.Normalize = 'on';
[f,fft_r1]=prepareCurveData(f,fft_r1);
[f,fft_r2]=prepareCurveData(f,fft_r2);
[f,fft_r3]=prepareCurveData(f,fft_r3);
[f,fft_source]=prepareCurveData(f,fft_source);
fit_r1=fit(f,fft_r1,ft,opts);
fit_r2=fit(f,fft_r2,ft,opts);
fit_r3=fit(f,fft_r3,ft,opts);
fit_source=fit(f,fft_source,ft,opts);
if ifshow
    figure(1);%Back to fft plot
    hold on
    plot(f,2*fit_r1(f),'r');
    plot(f,2*fit_r2(f),'g');
    plot(f,2*fit_r3(f),'c');
end
f=linspace(0.5*fp,1.5*fp,101);
f=f(1:end-1);
return_args.f=f;
return_args.r1=fit_r1(f)./fit_source(f);
return_args.r2=fit_r2(f)./fit_source(f);
return_args.r3=fit_r3(f)./fit_source(f);
end

function value=ricker_wavelet(cdtdx,np,it,md)
    value=pi^2*(cdtdx*it/np-md)^2;
    value=(1-2*value)*exp(-value);
end

function b=fin_diff(a,b,direction,scale_a,scale_b)
    switch direction
      case 'x'
        cols=1:size(b,2);
        b=full(scale_b.*b-scale_a.*(a(:,cols+1)-a(:,cols)));
      case 'y'
        rows=1:size(b,1);
        b=full(scale_b.*b-scale_a.*(a(rows,:)-a(rows+1,:)));
    end
end

