function [X,info] = TVNesterovRestore(B,M,eps_rel)
if nargin < 2
    error('Too few input parameters');
elseif nargin == 2
    eps_rel = 0.004;
end
if size(M) ~= size(B)
    error('The mask M must have the same dimensions as B');
end

Ic = int32( find(M(:) == 0) );
I  = int32( find(M(:) ~= 0) );

alpha = sum(B(Ic))/numel(Ic);
X = alpha*ones(size(B));
mdelta = norm(X-B,'fro');
if  mdelta < 0
    disp('The algorithm might experience numerical problem, because the');
    disp('solution X is almost alpha*ones(size(B)). Proceed with care');
end

R = max(B(Ic));
S = min(B(Ic));

gamma = 0.5*(R-S)*sqrt(length(I));
d = (R-S)/2 + S;

mn = numel(B);
epsilon = R*mn*eps_rel;
mu = epsilon/mn;
Lmu = 8/mu ;
N = int32( ceil(4*sqrt(2*mn)*sqrt(gamma^2)/epsilon) );

[X,~,~] = tv_restore(B,I,Ic,gamma,d,epsilon,Lmu,mu,N);
