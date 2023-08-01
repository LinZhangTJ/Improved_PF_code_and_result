% Load variables. The specific description about these variables see readme file.
% Load observation vector
load yy.mat;
% Load time vector
load time.mat; 
% The estimations of least square (without parameter denosing) as the initiations of the forward Kalman filtering
%load least_square.mat % Including variables xd Qxd, equal to 
[xl,Qxl]=least_square(time,y);

%Initiatize
M=size(y,1);    % The numbers of observations
n=size(time,1); % Months

% The estimations of the forward Kalman filtering
t=240;        % Correlation time
sig=10.^(14); % Variance
xxf0=xl;      % Initiatized parameters
PPf0=Qxl;     % Initiatized cariance matrix
[xxf,PPf,xrf,Prf]=Kalman_forward(time,y,t,sig,xxf0,PPf0);

% The final estimations of deterministic components without requiring further backward Kalman filtering
xd=xxf(1:7*M,1);
Pd=PPf(1:7*M,1:7*M);

% The estimations of the backward Kalman filtering
% Initialize
xrb0=xxf;
Prb0=PPf;
[xrb,Prb]=Kalman_backward(time,y,t,sig,xrb0,Prb0);

% Take the average values of the forward and backward Kalman refer to Ji&Herring.(2013)
for i=1:n
    K=diag(Prf(:,i))*inv(diag(Prb(:,i))+diag(Prf(:,i)));
    % The final estimations of irregular components
    xr(:,i)=xrf(:,i)+K*(xrb(:,i)-xrf(:,i));
    Pr(:,i)=diag(Prf(:,i))-K*diag(Prf(:,i));
end

% Paramters denosing

% Ymn is the matrix of spherical harmonic representation from spatial domain (1°×1°) to frequency domain (degree 60)
for j=1:360
lat(1:180,j)=89.5:-1:-89.5;
end
for j=1:180
lon(j,1:360)=0.5:1:359.5;
end
lat1=reshape(lat',180*360,1);
lon1=reshape(lon',180*360,1);
lat_lon_new=[lat1 lon1];
maxdeg=60;
for i=1:size(lat_lon_new,1)
Ymn(:,i)=grid2cs_coefficient(lat_lon_new(i,1),lat_lon_new(i,2),maxdeg);%无量纲
end
% Removing degree 0 and 1
Ymn1=Ymn(4:end,:);

lamda=10.^(-18); %Balance factor
% Denoising deterministic parameters
xd_filter=denoise_d(xd,Pd,Ymn1,lamda);
% Denoising irregular parameters
xr_filter=denoise_r(xr,Pr,Ymn1,lamda);

% Parameter reconstruction 
t0=time(1);
for i=1:n
% Observation matrix
A0=[1*IM (tt(i)-t0)*IM ((tt(i)-t0).^2)*IM cos(2*pi*(tt(i)-t0))*IM sin(2*pi*(tt(i)-t0))*IM cos(4*pi*(tt(i)-t0))*IM sin(4*pi*(tt(i)-t0))*IM];
% Reconstrction
y_new=A0*xd_filter+xr_filter(:,i); % The order of mass SHCs: {S22,S21,C20,C21,C22;S33,S32,S31,C30,C31,C32,C33;...}
cs=vector2cs(y_ner,60);            % Re-oreder mass SHCs to matrix |C\S|, and the maximum degree is 60
% Add back degree-one 
cs(1,2)=Degree_one(i,1,2);
cs(2,1)=Degree_one(i,2,1);
cs(2,2)=Degree_one(i,2,2);
grid_filter(:,:,i)=cs2grid(cs);
end


 