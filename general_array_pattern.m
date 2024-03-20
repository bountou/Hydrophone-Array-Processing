%%% Generate the beam patterns for linear, planar and circular arrays

clear
close all 
clc

%% Parameters 

c0     = 1480;       % speed of sound
f0     = 1480;       % frequency to plot
k0     = 2*pi*f0/c0; % wavenumber
lambda = c0/f0;      % wavelength


%% Arbitrary linear array 

% xm = zeros(1,M); 
xm = [-1.4 -1.2 -0.8 -0.5 -0.1 0.2 0.3 0.6 0.9 1.3 1.4 ]; % sensors along x-axis

M = length(xm);

ym = zeros(1,M);
% ym = [-1.4 -1.2 -0.8 -0.5 -0.1 0.2 0.4 0.6 0.9 1.3 1.4]; % sensors along y-axis

zm = zeros(1,M);
% zm = [-1.4 -1.2 -0.8 -0.5 -0.1 0.2 0.4 0.6 0.9 1.3 1.4]; % sensors along z-axis

d = min(diff(xm)); % minimum sensor spacing
f_al = c0/2/d; % aliasing frequency

pm = [xm ;ym ;zm]; % coordinates / sensors along x-axis


%% Uniform linear array 

% % d = lambda/2; % sensor spacing
% d = 0.5;
% f_al = c0/2/d; % aliasing frequency
% 
% 
% M  = 8; % number of sensors
% m  = (0:M-1);
% 
% % xm = zeros(1,M);
% xm = (m - (M-1)/2)*d; % sensors along x-axis
% 
% ym = zeros(1,M);
% % ym = (m - (M-1)/2)*d; % sensors along y-axis
% 
% zm = zeros(1,M);
% % zm = (m - (M-1)/2)*d; % sensors along z-axis
% 
% pm = [xm ;ym ;zm]; % coordinates / sensors along x-axis


%% Planar array on xy-plane

% % dx     = lambda/2; % sensor spacing along x-axis
% dx = 0.5;
% % dy     = lambda/2;
% dy = 0.4;
% 
% Mx = 8; % number of sensors along x-axis
% My = 6; % number of sensors along y-axis
% M = Mx*My;
% 
% 
% f_al = c0/2/min(dx,dy); % aliasing frequency
% 
% mx  = (0:Mx-1);
% x_linear= (mx - (Mx-1)/2)*dx;
% 
% my  = (0:My-1);
% y_linear = (my - (My-1)/2)*dy;
% 
% xm = repmat(x_linear,1,My);
% ym = repelem(y_linear,Mx);
% 
% zm = zeros(1,M);
% 
% 
% pm = [xm ;ym ;zm]; % coordinates / sensors along x-axis


%% Planar array on xz-plane

% % dx     = lambda/2; % sensor spacing along x-axis
% dx = 0.5;
% % dz     = lambda/2;
% dz = 0.4;
% 
% Mx = 5; % number of sensors along x-axis
% Mz = 4; % number of sensors along z-axis
% M = Mx*Mz;
% 
% 
% f_al = c0/2/min(dx,dz); % aliasing frequency
% 
% mx  = (0:Mx-1);
% x_linear= (mx - (Mx-1)/2)*dx;
% 
% ym = zeros(1,M);
% 
% mz  = (0:Mz-1);
% z_linear = (mz - (Mz-1)/2)*dz;
% 
% xm = repmat(x_linear,1,Mz);
% zm = repelem(z_linear,Mx);
% 
% 
% pm = [xm ;ym ;zm]; % coordinates / sensors along x-axis


%% Circular array on xy-plane

% R = 0.65; % radius
% M  = 8; % number of sensors
% m  = 0:M-1;
% mic_azi = (2*pi*m/M); % amizuth of sensor positions
% mic_azi_deg = 180/pi*mic_azi;
% [xm,ym] = pol2cart(mic_azi,R);
% 
% d = 2*R*sin(pi/M); % sensor spacing
% f_al = c0/2/d;     % aliasing frequency
% 
% zm = zeros(1,M);
% 
% 
% pm = [xm ;ym ;zm]; % coordinates / sensors along x-axis


%% Circular array on xz-plane

% R = 0.65; % radius
% M  = 8; % number of sensors
% m  = 0:M-1;
% mic_azi = (2*pi*m/M); % amizuth of sensor positions
% mic_azi_deg = 180/pi*mic_azi;
% [xm,zm] = pol2cart(mic_azi,R);
% 
% d = 2*R*sin(pi/M); % sensor spacing
% f_al = c0/2/d;     % aliasing frequency
% 
% ym = zeros(1,M);
% 
% 
% pm = [xm ;ym ;zm]; % coordinates / sensors along x-axis


%% Beamforming

% Incident wave directions
theta = -pi/2:pi/360:pi/2; % elevation (as in Matlab)
theta_deg = theta*180/pi;
phi   = -pi:pi/360:pi;     % azimuth
phi_deg = phi*180/pi;
[PHI, THETA] = meshgrid(phi,theta);


% Steering direction
theta0 = 0*pi/180;
phi0   = 90*pi/180;

u0 = [cos(theta0)*cos(phi0); cos(theta0)*sin(phi0); sin(theta0)];
vs = exp(1i*k0*pm'*u0);
  
% Spatial window
w = ones(M,1); % rectangular window
% w = hann(M); 

w_DS = w.*vs/(vs'*diag(w)*vs); % delay-and-sum weights

% Pattern (proportional to amplitude)
Y = zeros(size(PHI)); 
for q=1:length(theta)
    for l=1:length(phi)
        
       u = [cos(THETA(q,l))*cos(PHI(q,l)); cos(THETA(q,l))*sin(PHI(q,l)); sin(THETA(q,l))];
       vk = exp(1i*k0*pm'*u); % array manifold vector

       Y(q,l) = w_DS'*vk;
    end
end

%% Array performance metrics

%%%%%%%%%%%%%% WNG
% WNG_theor= 10*log10(M); % theoretical value for linear arrays

% WNG = 10*log10((vs'*vs)^2/(vs'*vs)); % for rectangular shading window

WNG = 10*log10((w_DS'*vs)^2/(w_DS'*w_DS)); % for any shading window



%% Plots  


% Array geometry
figure (1)
plot3(pm(1,:),pm(2,:),pm(3,:),'ko','MarkerSize',8)
grid on
xlabel('x')
ylabel('y')
zlabel('z')
set(gca, 'fontsize', 20, 'linewidth',2)
set(findall(gcf, 'Type', 'Line'),'LineWidth',2); 
view(120,30)


% 2D plots
th = 0*pi/180; % plotting on the horizontal plane

[~,ind]=min(abs(theta-th));
Y_xy = Y(ind,:);

figure (10)
polarplot(phi,abs(Y_xy),'k')
% set(gca,'ThetaZeroLocation', 'top','ThetaDir','clockwise')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
set(findall(gcf, 'Type', 'Line'),'LineWidth',2); 
set(gca, 'fontsize', 20, 'linewidth',2)
rlim([0 1])



figure(11)
% polardB(phi,20*log(abs(Y_xy)),[-40 0],5,'k',1)
polardB(phi,20*log(abs(Y_xy)),[-50 0],6,'k',1)
% set(gca,'ThetaZeroLocation', 'top','ThetaDir','clockwise')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(findall(gcf, 'Type', 'Line'),'LineWidth',2);
set(gca, 'fontsize', 20, 'linewidth',2)



% 3D plot
% [x,y,z] = sph2cart(PHI,THETA,abs(Y)); % linear scale
% dist = sqrt(x.^2 + y.^2 +z.^2);


[x,y,z] = sph2cart(PHI,THETA,20*log10(abs(Y)+1e-2)+40); % dB scale
dist = sqrt(x.^2 + y.^2 +z.^2);

figure (20)
s=surf(x,y,z,dist);
xlabel('x')
ylabel('y')
zlabel('z')
set(gca, 'fontsize', 16, 'linewidth',2)
s.EdgeColor = 'none';
% axis([-1 1 -1 1 -1 1]) 
axis equal   
view(120,30)
colormap(jet)
