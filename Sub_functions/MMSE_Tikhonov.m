function x_Tikhonov = MMSE_Tikhonov( sigma2,A, y, P, R, x_prior)

% Input data:
% sigma2: Unit Weight Variance
% A: Design Matrix
% y: Observation Vector
% P: Weight Matrix of Observation
% R: Regularization Matrix
% x_prior: Prior Parameter

% OUTPUTS:
% x_Tikhonov: estimated parameter

%% Initialization
[m, n] = size(A);
x_prior1 = x_prior;

%% Normalized equation
[Uc, Sc, Vc] = svd(R);
Vc = 1/2*(Uc+Vc);
R1= Vc*diag(sqrt(diag(Sc)))*Vc'; 
invR1=Vc*diag(sqrt(diag(1./Sc)))*Vc'; 
P1 = chol(P);     
y1 = P1 * y;
A1 = P1 * A * invR1;
x_prior = R1 * x_prior;

%% SVD decomposition
[U, S, V] = svd(A1);
S = S(1:min(m,n), :);
U = U(:, 1:min(m,n));
Lambda = diag(S);
    
%% regularization parameter estimation based on MMSE
alpha = 1;
left_alpha = 0;
right_alpha = 1e20;
times = 1000;
for iter = 1:times
first_deri = sum(Lambda.^2.*(alpha*(V'*x_prior).^2.-sigma2)./(Lambda.^2.+alpha).^3);        
if first_deri>=0
right_alpha = alpha;
alpha = (alpha+left_alpha)/2;            
else
left_alpha = alpha;
alpha = (alpha+right_alpha)/2;
end
end
Alpha = alpha * ones(n,1);

%% Tikhonov
x_Tikhonov = V * diag(1./diag(S'*S+ diag(Alpha))) * S * U' * y1;
x_Tikhonov  = invR1 * x_Tikhonov ;
end

