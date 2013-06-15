function return_args=wedge_fine(ifshow,k_max,m,fp,N)
if nargin<1
    close all
    clc
    ifshow=1;
    fp=68*9;
    k_max=5637*fp/500;
    m=1.4;
    N=3; % Discretization in fine grid
end
%% INITIALIZATION OF GRID AND SIMULATING PARAMETERS
c=340;          % Speed of sound
lp=c/fp;        % Associated wavelength with peak frequency
md=1;           % dr=md/fp
d=1.0;          % d distance in exercise (in m)
dx=c/fp/20;     % Spatial discretisation step, at least 20 samples
                % per wavelength
nd=ceil(d/dx);  % number of grids in d.
if mod(nd,2)==1 % nd should be even, because nd/2 is needed.
    nd=nd+1;
end
if mod(N,2)==0 % N should be odd
    N=N+1;
end
N2=(N-1)/2;          % N/2
dx=d/nd;             % Recalculate dx
dy=dx;               % Same grid in x,y
np=lp/dx;            % Does not need to be integer
nx=7*nd+1;           % Number of cells in x direction
ny=6*nd+1;           % Number of cells in y direction
npml=40;             % Number of cells in PML
xdim=nx-2+2*npml;    % Total number of columns for course grid
xdim_f=(nd+1)*N;     % Total number of columns for fine grid
ydim=ny-2+2*npml;    % Total number of rows
ydim_f=(nd/2+1)*N;   % Total number of rows for fine grid
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
col_r1=col_s+3*nd;   % X location of recorder1
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
ox_f=zeros(ydim_f,xdim_f+1);
oy_f=zeros(ydim_f+1,xdim_f);
p_f=zeros(ydim_f,xdim_f);
px=spalloc(ydim,xdim,2*npml*ydim+(nx-2)*2*npml);
py=px; % px and py are sparse. They fit around p.
% Frame is just made to better see the boundaries and source
% Frame is prealocated for more elements that actually needed.
frame=sparse(ydim,xdim,3*ydim+3*xdim);
frame(npml,:)=1;frame(ny+npml-1,:)=1;
frame(:,npml)=1;frame(:,nx+npml-1)=1;
frame(npml+ny-1-2*nd:npml+ny-1,npml+2*nd:npml+3*nd)=1;
slope=tril(ones(nd/2+1)); % Will be used to make the frame
wedge=[fliplr(slope(:,2:end)),slope];
frame(npml+ny-1-2.5*nd:npml+ny-1-2*nd,npml+2*nd:npml+3*nd)=wedge;
if ifshow % Prepare plots
    figure(2);
    coarse_h=subplot(1,3,[1,2]);
    set(coarse_h,'NextPlot','ReplaceChildren')
    set(coarse_h,'LooseInset',get(coarse_h,'TightInset'))
    axis(coarse_h,'equal')
    fine_h=subplot(1,3,3);
    set(fine_h,'NextPlot','ReplaceChildren')
    set(fine_h,'LooseInset',get(fine_h,'TightInset'))
    axis(fine_h,'equal')
    getframe();
end
% Calculate Row and Cols of fine grid in coarse grid
Row_f=[npml+ny-1-2.5*nd,npml+ny-1-2*nd];    % fine row start and end
Col_f=[npml+2*nd,npml+3*nd];                % fine col start and end
slope=triu(ones(nd/2+1),1); % Will be used to make the frame
wedge=[fliplr(slope(:,2:end)),slope];
slope_f=diag(ones(1,N*nd/2+1));
frame_f=[fliplr(slope_f(:,2:end)),slope_f];
frame_f=[zeros(N*nd/2+1,N2),frame_f,zeros(N*nd/2+1,N2)];
frame_f=[zeros(N2,N*(nd+1));frame_f;zeros(N2,N*(nd+1))];
slope_f=triu(ones(N*nd/2+1),1); % Will be used on fine grid 
wedge_f=[fliplr(slope_f(:,2:end)),slope_f];
wedge_f=[ones(N*nd/2+1,N2),wedge_f,ones(N*nd/2+1,N2)];
wedge_f=[ones(N2,N*(nd+1));wedge_f;ones(N2,N*(nd+1))];
wedge_f(end-N2+1:end,N2+1:end-N2)=0;
mask=zeros(N);
mask(N2+1,N2+1)=1;
mask=repmat(mask,nd/2+1,nd+1);
mask=sub2ind(size(mask),find(mask)); % Will be used later to address
% coarse values in fine grid
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
KP_y=spalloc(ydim,xdim,npml*ydim);
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
KP_x_o=c^2*dt/dx;
KP_x_p=1-KP_x*dt;
KP_y_o=c^2*dt/dy;
KP_y_p=1-KP_y*dt;
%% STARTING OUT WITH TIME LOOPS
display='it=0/nt t=0.000000';
fprintf(display);
for it_c=1:nt
    %% SOURCE AND RECORDERS
    t=(it_c-1)*dt;
    fprintf(repmat('\b',1,length(display)))
    display=sprintf('it=%d/%d t=%.6fs\n',it_c,nt,t);
    fprintf(sprintf(display))
    p(row_s,col_s)=p(row_s,col_s)+source(it_c); % Apply source
    r1(it_c)=p(row_s,col_r1); % Record the values
    r2(it_c)=p(row_s,col_r2); % Record the values
    r3(it_c)=p(row_s,col_r3); % Record the values
    %% INTERPOLATE P_FINE ON EDGE BILINEAR
    coarse_data=p(Row_f(1)-1:Row_f(2)+1,Col_f(1)-1:Col_f(1));
    [X,Y]=meshgrid([0,1],0:nd/2+2);
    XI=0.5-1/(2*N);
    YI=0.5+1/(2*N):1/N:nd/2+1.5-1/(2*N);
    left_f=interp2(X,Y,coarse_data,XI,YI);
    coarse_data=p(Row_f(1)-1:Row_f(2)+1,Col_f(2):Col_f(2)+1);
    [X,Y]=meshgrid([0,1],0:nd/2+2);
    XI=0.5+1/(2*N);
    YI=0.5+1/(2*N):1/N:nd/2+1.5-1/(2*N);
    right_f=interp2(X,Y,coarse_data,XI,YI);
    coarse_data=p(Row_f(1)-1:Row_f(1),Col_f(1)-1:Col_f(2)+1);
    [X,Y]=meshgrid(0:nd+2,[0,1]);
    XI=0.5+1/(2*N):1/N:nd+1.5-1/(2*N);
    YI=0.5-1/(2*N);
    top_f=interp2(X,Y,coarse_data,XI,YI);
    coarse_data=p(Row_f(2):Row_f(2)+1,Col_f(1)-1:Col_f(2)+1);
    [X,Y]=meshgrid(0:nd+2,[0,1]);
    XI=0.5+1/(2*N):1/N:nd+1.5-1/(2*N);
    YI=0.5+1/(2*N);
    bottom_f=interp2(X,Y,coarse_data,XI,YI);
    %% UPDATE FINE GRID OX_F OY_F P_F AND CALCULATE CORRECTION FACTOR
    df_l=-c^2*dt/dx*ox(Row_f(1):Row_f(2),Col_f(1)-1);
    df_r= c^2*dt/dx*ox(Row_f(1):Row_f(2),Col_f(2));
    df_t= c^2*dt/dy*oy(Row_f(1),Col_f(1):Col_f(2));
    df_b=-c^2*dt/dy*oy(Row_f(2)+1,Col_f(1):Col_f(2));
    for it_f=1:N
        ox_f=fin_diff([left_f,p_f,right_f],ox_f,'x',dt/dx);
        oy_f=fin_diff([top_f;p_f;bottom_f],oy_f,'y',dt/dy);
        p_f=fin_diff(ox_f,p_f,'x',c^2*dt/dx);
        p_f=fin_diff(oy_f,p_f,'y',c^2*dt/dy);
        p_f=wedge_f.*p_f;  % APPLY WEDGE BOUNDRY ON FINE GRID
        for it=0:nd/2
            df_l(it+1)=df_l(it+1)+c^2*dt/dx/(N^2)*sum(ox_f(1+it*N:(it+1)*N,1));
            df_r(it+1)=df_r(it+1)-c^2*dt/dx/(N^2)*sum(ox_f(1+it*N:(it+1)*N,end));
            df_t(it+1)=df_t(it+1)-c^2*dt/dy/(N^2)*sum(oy_f(1,1+it*N:(it+1)*N));
            df_b(it+1)=df_b(it+1)+c^2*dt/dy/(N^2)*sum(oy_f(end,1+it*N:(it+1)*N));
        end
    end   
    %% UPDATE COARSE GRID OX OY P
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
    %%  REPLACE COARSE VALUES WITH OVERLAPPING FINE VALUES
    %p(Row_f(1):Row_f(2),Col_f(1):Col_f(2))=reshape(p_f(mask),nd/2+1,nd+1);
    for it_fr=0:nd/2
        for it_fc=0:nd
            p(Row_f(1)+it_fr,Col_f(1)+it_fc)=sum(sum(p_f(1+it_fr*N:...
                (it_fr+1)*N,1+it_fc*N:(it_fc+1)*N)))/(N^2);
        end
    end
    %% APPLY CORRECTION FACTOR FROM FINE GRID ON COARSE GRID
    p(Row_f(1):Row_f(2),Col_f(1)-1)=p(Row_f(1):Row_f(2),Col_f(1)-1)-df_l;
    p(Row_f(1):Row_f(2),Col_f(2)+1)=p(Row_f(1):Row_f(2),Col_f(2)+1)-df_r;
    p(Row_f(1)-1,Col_f(1):Col_f(2))=p(Row_f(1)-1,Col_f(1):Col_f(2))-df_t;
    p(Row_f(2)+1,Col_f(1):Col_f(2))=p(Row_f(2)+1,Col_f(1):Col_f(2))-df_b;
    % Apply Boundry Condition on Coarse Grid
    p(npml+ny-1,:)=0;
    px(npml+ny-1,:)=0;
    py(npml+ny-1,:)=0;
    p(npml+ny-2*nd-1:npml+ny-1,npml+2*nd:npml+3*nd)=0; % Thick Wall
    p(npml+ny-1-2.5*nd:npml+ny-1-2*nd,npml+2*nd:npml+3*nd)=wedge.*...
        p(npml+ny-1-2.5*nd:npml+ny-1-2*nd,npml+2*nd:npml+3*nd); % Wedge
    %% VISUALIZATION
    if ifshow && mod(it_c,10)==0
        pcolor(coarse_h,flipud(p+frame+px+py))
        shading(coarse_h,'interp')
        colormap(coarse_h,'gray')
        caxis(coarse_h,[-.01 .01]) %Assuming A=1
        xlim(coarse_h,[1 xdim])
        ylim(coarse_h,[1 ydim])
        title(coarse_h,sprintf('it=%d/%d, fp=%d, t=%.5f',it_c,nt,fp,t))
        
        pcolor(fine_h,flipud(p_f+frame_f))
        shading(fine_h,'interp')
        colormap(fine_h,'gray')
        caxis(fine_h,[-0.01 0.01])
        xlim(fine_h,[1 xdim_f])
        ylim(fine_h,[1 ydim_f])
        title(fine_h,sprintf('Fine Descritization, N=%d',N))
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
if nargin==4
    scale_b=1;
end
switch direction
    case 'x'
        cols=1:size(b,2);
        b=full(scale_b.*b-scale_a.*(a(:,cols+1)-a(:,cols)));
    case 'y'
        
        rows=1:size(b,1);
        b=full(scale_b.*b-scale_a.*(a(rows,:)-a(rows+1,:)));
end
end

