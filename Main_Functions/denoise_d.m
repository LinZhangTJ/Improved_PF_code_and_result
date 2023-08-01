function xd_filter=denoise_d(xd,Pd,Ymn1,lamda)

% Input data:
% xd: unfiltered determinitic parameters after Kalman filtering
% Pd: the covariance matrix of xd
% Ymn1: the matrix of transformation from grid to SHC. Grid is 180*360 [89.5:-89.5 0.5:359.5]
% lamda: balance factor for the covariance matrices of signals and noise

% Output data:
% xd_filter: denoising deterministic parameter

%Initialization
M=size(xd,1)/7;
gnum=64800; 
num=1;
ee(num)=1;
gm0=zeros(gnum,1);
gm=zeros(gnum,7);
%Scale Matrix
D=get_scaleM(60);
% Initial signal covariance
C=D*Ymn1;
Sy0=C*C';
Qs=kron(eye(7),Sy0);

% Iteration starts
while ee(num)>=0.002
num=num+1;
mmk=num-1; % Iteration numbers
gm0=gm(:,2);
Qn=Pd; % Noise covariance
Syn=Qs+lamda*Qn; 
invsyn=(Syn)\eye(size(Syn));
zz(:,mmk)=Qs*invsyn*xd;% Filtered paramters at each iteration

for j=1:7
    csm=vector2cs(zz((j-1)*M+1:j*M,mmk),60);
    grid=cs2grid(csm);
    gg0=reshape(grid',180*360,1);
    gg(:,j)=gg0.^2; % Unit: m
end
gm(:,:)=gg; 
% Update signal covariance
for jj=1:7
Qsg=gm(:,jj);
Qsg1=Qsg.^(1/2);
C=D*Ymn1;
S1=C.*Qsg1';
Qs0(:,:,jj)=S1*S1';
end
Qs=blkdiag(Qs0(:,:,1),Qs0(:,:,2),Qs0(:,:,3),Qs0(:,:,4),Qs0(:,:,5),Qs0(:,:,6),Qs0(:,:,7));
% Iteration differences
ee(num)=(abs(mean(gm(:,2))-mean(gm0)))/(mean(gm(:,2)));
end
xd_filter=zz;
end