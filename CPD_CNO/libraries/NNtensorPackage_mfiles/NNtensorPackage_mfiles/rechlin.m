function mu=rechlin(X, M, gradient, d, nbFactors, penalization, alpha, Miter, mu0)

% This function is part of the positive tensor package
% P. Comon
%
% This package was developped according to the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors
%
% -------------------------------
% mu=rechlin(X, M, gradient, d, nbFactors, penalization, alpha, Miter, mu0)
%
% X: 3rd order tensor.
% M: Loading matrices A, B and C stacked and vectorized.
% d: Search direction.
% nbFactors: Desired tensor rank.
% penalization: Desired penalization (1 for L1 penalization, 2 for L2
% penalization, 0 for no penalization).
% alpha : Value of the penalty term.
%
% Optional:
% Miter : Number of iterations to compute the backtracking
% mu0: Initial value of the stepsize mu
%
% Returns the stepsize mu

if nargin<8,Miter=20;end;
if nargin<9,mu0=1;end;

dimCube = size(X);
[A, B, C] = extractLoadings(M, dimCube, nbFactors);
[directionA, directionB, directionC] = extractLoadings(d, dimCube, nbFactors);

% XDep = deplier(X, [1 2 3]);

gamma=0.01;k=1;seuil=1;dimin=0;
mu=mu0;

Yold = computeError(X, A .* A, B .* B, C .* C, penalization, alpha);
%Yold=sum(sum((XDep - (A .* A) * kr(C .* C, B .* B).').^2));
gold=gradient;

while (k<Miter)&&(dimin<seuil)
%   modele = ((A + mu * directionA) .* (A + mu * directionA)) * ...
%     kr((C + mu * directionC) .* (C + mu * directionC), ...
%        (B + mu * directionB) .* (B + mu * directionB)).';
%    
%   Ynew=sum(sum((XDep - modele).^2));
  Ynew = computeError(X, ...
                      (A + mu * directionA) .* (A + mu * directionA), ...
                      (B + mu * directionB) .* (B + mu * directionB), ...
                      (C + mu * directionC) .* (C + mu * directionC), ...
                      penalization,...
                      alpha);
  dimin=Yold-Ynew;
  seuil=-gamma*mu*trace(gold'*d); %Virer la trace
  mu=mu*0.05;k=k+1;
  %fprintf('k=%g, mu=%g, Ynew=%g \n',k,mu,Ynew);
end;