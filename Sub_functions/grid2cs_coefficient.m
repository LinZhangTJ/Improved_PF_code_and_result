function coe=grid2cs_coefficient(lat,lon,maxdeg) 

% Input data:
% lat: Latiude
% lon: Lontitude

% Output data:
% coe: The spherical harmonic representation from spatial domain to frequency domain

l=zeros(maxdeg+1,2*maxdeg+1);
m=zeros(maxdeg+1,2*maxdeg+1);
q=zeros(maxdeg+1,2*maxdeg+1);
for i=2:maxdeg+1
    l(i,maxdeg+1-(i-1):maxdeg+1+(i-1))=i-1;
    m(i,maxdeg:-1:maxdeg+1-(i-1))=1:i-1;
    m(i,maxdeg+1:maxdeg+1+(i-1))=0:i-1;
    m(i,maxdeg+1)=1000;
    q(i,maxdeg:-1:maxdeg+1-(i-1))=1;
    q(i,maxdeg+1:maxdeg+1+(i-1))=2;
end
l0=l(2:maxdeg+1,:);
l1=l0';
ll=reshape(l1,maxdeg*(2*maxdeg+1),1);
ll(find(ll==0))=[];

m0=m(2:maxdeg+1,:);
m1=m0';
ml=reshape(m1,maxdeg*(2*maxdeg+1),1);
ml(find(ml==0))=[];
ml(find(ml==1000))=0;

q0=q(2:maxdeg+1,:);
q1=q0';
ql=reshape(q1,maxdeg*(2*maxdeg+1),1);
ql(find(ql==0))=[];

lmq(:,1)=ll;
lmq(:,2)=ml;
lmq(:,3)=ql;

ae =  6378136.3;
rho_w = 1000;
rho_ave = 5517;
dens = ae*rho_ave/(3.*rho_w);
load 'HanWahrLoveNumbers.mat'
k(2) = 0.021;
th=90-lat;
cos_th=cos(pi/180*th);
load 'F_lm.mat';

    for n=0:maxdeg
    M_pnmth(n*(n+1)/2+1:(n+1)*(n+2)/2,:)=func_pnm(n,cos_th);
    end
    pnm=[lm M_pnmth];
    si=func_weight_area(th,th,1);

for i=1:size(lmq(:,1),1)
 
    for jj=1:size(pnm,1)
        if pnm(jj,1)==lmq(i,1)&&pnm(jj,2)==lmq(i,2)
            P(i)=pnm(jj,3);
        end
    end
    if lmq(i,3)==1
        S(i)=sin(pi/180*lmq(i,2)*lon);
        C(i)=0;
    else
        C(i)=cos(pi/180*lmq(i,2)*lon);
        S(i)=0;
    end
    coe(i)=1/4/pi*(si*P(i)*S(i)+si*P(i)*C(i));
end
end
