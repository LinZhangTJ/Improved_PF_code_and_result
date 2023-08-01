function [xxf,PPf,xrf,Prf]=Kalman_forward(time,y,t,sig,xxf0,PPf0)

% Input data:
% time: the study months in year
% y: the observations involving SHCs (monthly SHC numbers * months)
% t: the correlation time of FOGM
% sig: the variance of FOGM
% xxf0: initial parameter (usually obtain from least-square)
% PPf0: initial covariance matrix of parameter (usually obtain from least-square)

% Output data:
% xxf: Itegrated parameters involving determistic and irreegular parameters
% PPf: The covariance matrix of xx
% xrf: Irregular parameters
% Prf: The covariance matrices of xr


% Converts geoid coefficients (gc) to mass coefficients (mc)
T=diag(v_gc2mc(60));
ym=T*y*100; % Multiply 100 representing m into cm

%Initialize
M=size(y,1);% The SHCs numbers
IM=eye(M);
Q=diag(ones(1,8*M),0)*0.000; % Process noise
I=diag(ones(1,8*M),0);
t0=time(1);
X01=[xxf0;zeros(M,1)];           
P01=blkdiag(diag(diag(PPf0)),eye(M)*10.^(14)); 
%Design Matrix
H=[1*IM (time(1)-t0)*IM ((time(1)-t0)^2)*IM cos(2*pi*(time(1)-t0))*IM sin(2*pi*(time(1)-t0))*IM cos(4*pi*(time(1)-t0))*IM sin(4*pi*(time(1)-t0))*IM 1*IM];

% Prediction
load EE1.mat;  %The covariance matrix of original observation vector y
EE=T*EE1(:,:,1)*T'*10000; %The covariance matrix of observation vector ym
K=P01*H'*inv(H*P01*H'+EE(:,:,1)); % Gain matrix
% Update
v=ym(:,1)-H*X01;
xxf=X01+K*v;
PPf=(I-K*H)*P01;
Prf(:,1)=diag(PPf(7*M+1:8*M,7*M+1:8*M));
xrf(:,1)=xxf(7*M+1:8*M,1);

% The rest n-1 months
for i=2:length(time)
   X00=xxf;
   P00=PPf;
   H=[1*IM (time(i)-t0)*IM ((time(i)-t0)^2)*IM (cos(2*pi*(time(i)-t0)))*IM (sin(2*pi*(time(i)-t0)))*IM (cos(4*pi*(time(i)-t0)))*IM (sin(4*pi*(time(i)-t0)))*IM 1*IM];
   dt=(time(i)-time(i-1))*365.25;
   B=blkdiag(eye(7*M),eye(M)*diag(exp(-dt/t)));
   Q(7*M+1:8*M,7*M+1:8*M)=eye(M)*sig*(1-exp(-2*dt/t)); 
   % Prediction
   X01=B*X00;
   P01=B*P00*B'+Q;
% Due to the memory limitation of MATLAB software, in these codes, we divided the covariance matrices of observations from the Tongji-Grace2018 model during 2002-2016 into five parts: 
% EE1 (2002.02-2004.12); EE2 (2005.01-2007.06); EE3 (2007.07-2009.12); EE4 (2010.01-2012.09); EE5 (2012.10-2016.12).
if i<=30
      if ~exist('EE1','var')
         load EE1.mat;
      end 
      % Covariance propagation
      EE=T*EE1(:,:,i)*T'*10000; 
  
   elseif i>30&&i<=60
      if exist('EE1','var') && ~isempty(EE1)
          clear EE1;
      end
      if ~exist('EE2','var')
          load EE2.mat;
      end
      % Covariance propagation
      EE=T*EE2(:,:,i-30)*T'*10000; 
   
   elseif i>60&&i<=90
      if exist('EE2','var') && ~isempty(EE2)
           clear EE2;
      end
      if ~exist('EE3','var')
           load EE3.mat;
      end
      % Covariance propagation
      EE=T*EE3(:,:,i-60)*T'*10000; 
      
   elseif i>90&&i<=120
        if exist('EE3','var') && ~isempty(EE3)
             clear EE3;
        end
        if ~exist('EE4','var')
             load EE4.mat;
        end   
     % Covariance propagation
     EE=T*EE4(:,:,i-90)*T'*10000;   
                 
   elseif i>120
       if exist('EE4','var') && ~isempty(EE4)
             clear EE4;
       end         
       if ~exist('EE5','var')
             load EE5.mat;
       end 
     % Covariance propagation
     EE=T*EE5(:,:,i-120)*T'*10000; 
end
   %Gain matrix
   K=P01*H'*inv(H*P01*H'+EE);
   
   %Update
   v=ym(:,i)-H*X01;
   xxf=X01+K*v;
   PPf=(I-K*H)*P01;
   Prf(:,i)=diag(PPf(7*M+1:8*M,7*M+1:8*M));
   xrf(:,i)=xxf(7*M+1:8*M,1);
end
end
