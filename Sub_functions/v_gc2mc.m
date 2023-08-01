function D=v_gc2mc(N)

% Input data:
% N: The maximum degree of SHC

% Output data:
% D: Converts geoid coefficients (gc) to mass coefficients (mc)

l=zeros(N+1,2*N+1);
for i=3:N+1
    l(i,N+1-(i-1):N+1+(i-1))=i-1;
end
l0=l(3:N+1,:);
l1=l0';
ll=reshape(l1,(N-1)*(2*N+1),1);
ll(find(ll==0))=[];
lmq(:,1)=ll;

ae =  6378136.3;
rho_w = 1000;
rho_ave = 5517;
dens = ae*rho_ave/(3.*rho_w);
load 'HanWahrLoveNumbers.mat'
k(2) = 0.021;
load 'F_lm.mat';

for ll = 1:size(lmq,1)
    D0(ll)  = dens*(2*lmq(ll,1)+1 )/(1+k(lmq(ll,1)+1));
end

D=D0;

end



