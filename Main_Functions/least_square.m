function [xl,Qxl]=least_square(time,y)
% The least-square estimations don't incorporate irregular parameters 

% Input data:
% time: the study months in year
% y: the observations involving SHCs; size of the SHC numbers (e.g.,3717:degree 2-60) * months 

% Output data:
% xl: Estimated parameter
% Qxl: Covariance matrix of parameter

%Initialization
M=size(y,1);% the SHC numbers (degree 2-60)
n=size(y,2);% months
IM=eye(M);
N0=zeros(7*M,7*M);
b0=zeros(7*M,1);
t0=time(1);
tt=time;
% T: Converts geoid coefficients (gc) to mass coefficients (mc)
T=diag(v_gc2mc(60));
ym=T*y*100; %multiply 100 representing m into cm

for i=1:n
%Observation matrix
A0=[1*IM (tt(i)-t0)*IM ((tt(i)-t0).^2)*IM cos(2*pi*(tt(i)-t0))*IM sin(2*pi*(tt(i)-t0))*IM cos(4*pi*(tt(i)-t0))*IM sin(4*pi*(tt(i)-t0))*IM];
%Due to the memory limitation of MATLAB software, in these codes, we divided the covariance matrices of observations from the Tongji-Grace2018 model during 2002-2016 into five parts: 
%EE1 (2002.02-2004.12); EE2 (2005.01-2007.06); EE3 (2007.07-2009.12); EE4 (2010.01-2012.09); EE5 (2012.10-2016.12).
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
     EE=T*EE4(:,:,i-90)*T';   
                 
   elseif i>120
       if exist('EE4','var') && ~isempty(EE4)
             clear EE4;
       end         
       if ~exist('EE5','var')
             load EE5.mat;
       end 
     % % Covariance propagation
     EE=T*EE5(:,:,i-120)*T'; 
end
N0=N0+A0'*inv(EE)*A0;
b0=b0+A0'*inv(EE)*ym(:,i);
end

xl(:,1)=((N0)\eye(size(N0)))*b0;
Qxl=(N0)\eye(size(N0));
end