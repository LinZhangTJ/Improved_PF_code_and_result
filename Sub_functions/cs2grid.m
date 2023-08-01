function grid=cs2grid(cs0)

% Input data:
% cs0: SHCs in the frequence domain

% Output data:
% grid: spatial grids 180*360 [89.5:-89.5 0.5:359.5]

c11cmn=[0.5 89.5 359.5 -89.5];
degree_max=size(cs0,1)-1;
index_lmcosi=(degree_max+2)*(degree_max+1)/2;
lmcosi=zeros(index_lmcosi,4);

for kk=1:size(cs0,3)
    cs(:,:)=cs0(:,:,kk);
for ii=1:degree_max+1
    for jj=1:ii
        index_tmp=ii*(ii-1)/2+jj;
        lmcosi(index_tmp,1)=ii-1; % l
        lmcosi(index_tmp,2)=jj-1; % m
        lmcosi(index_tmp,3)=cs(ii,jj); % C_lm
        if jj==1   
            lmcosi(index_tmp,4)=0; % S_lm
        else
            lmcosi(index_tmp,4)=cs(jj-1,ii); % S_lm
        end
    end
end
end

C=lmcosi(:,3);
S=lmcosi(:,4);
th1=90-c11cmn(2);
th2=90-c11cmn(4);
lam1=c11cmn(1);
lam2=c11cmn(3);
d=1;
N=degree_max;

CC=C(1:(N+1)*(N+2)/2);
SS=S(1:(N+1)*(N+2)/2);
cos_lam=cos(pi/180*(0:N)'*(lam1:d:lam2)); 
sin_lam=sin(pi/180*(0:N)'*(lam1:d:lam2));
cos_th=cos(pi/180*(th1:d:th2));

n_lam=round((lam2-lam1)/d+1);
n_th=round((th2-th1)/d+1);
M_coslam=zeros((N+1)*(N+2)/2,n_lam);
M_sinlam=M_coslam;
M_pnmth=zeros((N+1)*(N+2)/2,n_th);
for n=0:N
    M_coslam(n*(n+1)/2+1:(n+1)*(n+2)/2,:)=cos_lam(1:n+1,:);
    M_sinlam(n*(n+1)/2+1:(n+1)*(n+2)/2,:)=sin_lam(1:n+1,:);
    M_pnmth(n*(n+1)/2+1:(n+1)*(n+2)/2,:)=func_pnm(n,cos_th);
end

grid(:,:,kk)=M_pnmth'*(CC*ones(1,(lam2-lam1)/d+1).*M_coslam+SS *ones(1,(lam2-lam1)/d+1).*M_sinlam);

end
   



