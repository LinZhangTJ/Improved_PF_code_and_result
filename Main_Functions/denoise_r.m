function xr_filter=denoise_r(xr,Pr,Ymn1,lamda)

% Input data:
% xr: unfiltered irregular parameters
% Pr: the covariance matrix of xr
% Ymn1: the matrix of transformation from grid to SHC. Grid is 180*360 [89.5:-89.5 0.5:359.5]
% lamda: Initial balance factor for the covariance matrices of signals and noise

% Output data:
% xr_filter: filtered irregular parameter

%Initialization
M=size(xr,1); % The numbers of SHCs 
n=size(xr,2); % Months
gnum=64800; 
gm=zeros(gnum,1);
num=1;
ee(1)=1; 
% Scale Matrix
D=get_scaleM(60);
% Initial signal covariance
C=D*Ymn1;
Qs=C*C';

% Iteration begins
while ee(num)>=0.002
num=num+1;
mmk=num-1 % Iteration numbers
gm0=gm(:,1); 
for ii=1:n
Qn=diag(Pr(:,ii)); % Noise covariance matrix
Syn=Qs+lamda*Qn; 
invsyn=(Syn)\eye(size(Syn));
zz(:,ii,mmk)=Qs*invsyn*xr(:,ii);
csm=vector2cs(zz(:,ii,mmk),60);
grid(:,:,ii)=cs2grid(csm);
end
% Update signal variance
for i=1:180
    for j=1:360
        ss(:,:)=grid(i,j,:);
        rmss(i,j)=sqrt(sum(ss.^2)/n);
    end
end
gg0=reshape(rmss',180*360,1);
gm=gg0.^2;
Qsg=gm;
Qsg1=Qsg.^(1/2);
C=D*Ymn1;
Qs0=C.*Qsg1';
Qs=Qs0*Qs0';
% Iteration differences
ee(num)=(abs(mean(gm(:,1))-mean(gm0)))/(mean(gm(:,1)))
end

%Adjustment factor
for i=1:n
[x_Tik,alpha]=MMSE_Tikhonov_1(10.^(-14),eye(M),xr(:,i),inv(diag(Pr(:,ii))),inv(Qs),zz(:,i));
xx(:,i)=x_Tik;
end

xr_filter=xx;
end
