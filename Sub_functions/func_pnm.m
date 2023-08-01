function y=func_pnm(n,x)

% Input data:
% n: the degree of SHCs
% x: co-latitude

% Output data
% y: nomalized Plm

P=legendre(n,x,'norm');
y=[sqrt(2)*P(1,:);2*P(2:n+1,:)];
end