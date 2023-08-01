function [xxb,PPb,xrb,Prb]=Kalman_backward(time,y,t,sig,xb0,Pb0)

% Input data:
% time: the study months in year
% y: the observations involving SHCs; size of the SHC numbers (e.g.,3717:degree 2-60) * months 
% t: the correlation time of FOGM
% sig: the variance of FOGM
% xrb0: initial parameter (the estimations in the final month of the forward Kalman filtering)
% Prb0: initial covariance matrix of parameter (the estimations in the final month of the forward Kalman filtering)

% Output data:
% xx: Itegrated parameters involving determistic and irreegular parameters
% PP: The covariance matrix of xx
% xr: Irregular parameters
% Pr: The covariance matrices of xr

% T: Converts geoid coefficients (gc) to mass coefficients (mc)
T=diag(v_gc2mc(60));
ym=T*y*100; % Multiply 100 representing m into cm

% Initialize
t0=time(1);
M=size(y,1);% The SHCs numbers
n=size(y,2);% Months 
IM=eye(M);
Q=diag(ones(1,8*M),0)*0.000; % Process noise
I=diag(ones(1,8*M),0);
xxb=xb0;
PPb=Pb0;
Prb(:,n)=diag(PPb(7*M+1:8*M,7*M+1:8*M));
xrb(:,n)=xxb(7*M+1:8*M,1);

% The rest n-1 months
for i=n-1:-1:1
   X00=xxb;
   P00=PPb;
   H=[1*IM (time(i)-t0)*IM ((time(i)-t0)^2)*IM (cos(2*pi*(time(i)-t0)))*IM (sin(2*pi*(time(i)-t0)))*IM (cos(4*pi*(time(i)-t0)))*IM (sin(4*pi*(time(i)-t0)))*IM 1*IM];
   dt=(time(i+1)-time(i))*365.25;
   B=blkdiag(eye(7*M),eye(M)*diag(exp(dt/t)));% The inverse of B in the forward Kalman process
   Q(7*M+1:8*M,7*M+1:8*M)=eye(M)*sig*(1-exp(-2*dt/t)); 
   % Prediction
   X01=B*X00;
   P01=B*(P00+Q)*B';
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
   % Gain matrix
   K=P01*H'*inv(H*P01*H'+EE);
   % Update
   v=ym(:,i)-H*X01;
   xxb=X01+K*v;
   PPb=(I-K*H)*P01;
   Prb(:,i)=diag(PPb(7*M+1:8*M,7*M+1:8*M));
   xrb(:,i)=xxb(7*M+1:8*M,1);
end
end













