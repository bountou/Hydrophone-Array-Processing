clear
close all
clc


%% Array geometry

% Uniform linear array

Q = 12; % number of sensors
q = (0:Q-1)';
d = 0.5; % inter-sensor distance
zn = ((q-(Q-1)/2)*d)'; % sensor positions on z-axis

xn = zeros(Q,1);
yn = zeros(Q,1);
pn = [xn yn zn']; % microphone positions


% source direction
doa_elev_deg = 60;
doa_elev = doa_elev_deg*pi/180; % elevation follows Matlab definition [-90, 90]
doa_azi_deg = 0;
doa_azi= doa_azi_deg*pi/180; 

[x,y,z] = sph2cart(doa_azi, doa_elev, 1); % cartesian coordinates
U_doa = [x y z]; % unit vector

% Impulse response parameters
Fs = 8000;
Lfilt = 2000;

% simulate array response using getArrayResponse()
fDirectivity = @(angle) 1; % response of omnidirectional sensor
[h_mic, H_mic] = getArrayResponse(U_doa, pn, [], fDirectivity, Lfilt, Fs); % sensor orientation irrelevant in this case


%% %%%%%%%%%%%%%%%%%%---- Static source -----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dur = 3; % duration in seconds
N = round(Fs*dur); % duration in samples

Ps = 1; % signal power 
sig = sqrt(Ps)*randn(N,1);


xs = zeros(N,Q); % array signals
for i=1:Q
    xs(:,i) = fftfilt(h_mic(:,i),sig);
end

%% %%%%%%%%%%%%%%%%%---- White noise ----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR = 0;

Pn = Ps/(10^(SNR/10)); % noise power at each sensor

sn = sqrt(Pn)*randn(N,Q);

x_tot = xs + sn;
