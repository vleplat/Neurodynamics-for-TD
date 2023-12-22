function error = computeError(X, A, B, C, penalization, alpha)

% This function is part of the positive tensor package.
% June, 2011.
% J-P. Royer, P. Comon, N. Thirion-Moreau.
%
% This package was developped according to the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors.
%
% ----------------------------
% error = computeError(X, A, B, C, penalization, alpha) computes the
% squared error between the data X and the model.
% A penalization can be added to the cost function. A value of 0 means no
% penalization, 1 means a L-1 penalization, and 2 a L-2 penalization.
% Alpha is the penalty coefficient.

XUnf = unfold(X, [1 2 3]);
model = (A) * khatriRao(C, B).';

if(penalization == 0)
    error=sum(sum((XUnf - model).^2));
elseif(penalization == 1)
    error=sum(sum((XUnf - model).^2)) + alpha * sum(sum(A)) + alpha * sum(sum(B)) + alpha * sum(sum(C));
%    modele = (A .* A) * kr(C .* C, B .* B).';
%    erreur=sum(sum((XDep - modele).^2)) + alpha * norm(A .* A, 1) + alpha * norm(B .* B, 1) + alpha * norm(C .* C, 1);
else
    error=sum(sum((XUnf - model).^2)) + alpha * norm(A, 'fro')^2 + alpha * norm(B, 'fro')^2 + alpha * norm(C, 'fro')^2;
    %error('non fait');
end