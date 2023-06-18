function [gradA, gradB, gradC] = computeGradientsPositivity(X, A, B, C, penalization, alpha)

% This function is part of the positive tensor package.
% June, 2011.
% J-P. Royer, P. Comon, N. Thirion-Moreau.
%
% This package was developped according to the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors.
%
% --------------------------------
% [gradA, gradB, gradC] = computeGradientsPositivity(X, A, B, C,
% penalization, alpha)
%
% Computes the three gradient matrices.
%
% X is the 3rd order tensor.
% A, B, C are the 3 loading matrices.
% penalization is 1 or 2 if applied to the problem (L1 or L2 penalty), 
% 0 if not.
% alpha is a small scalar defining the penalty coefficient (if applied).

if(penalization == 0)
    alpha = 0;
end

mode1 = [1 2 3];
mode2 = [2 3 1];
mode3 = [3 1 2];

%1er dépliement
XUnf = unfold(X, mode1);
gradA = gradientCPPositivity(A, XUnf, C, B);

%2ème dépliement
XUnf = unfold(X, mode2);
gradB = gradientCPPositivity(B, XUnf, A, C);

%3ème dépliement
XUnf = unfold(X, mode3);
gradC = gradientCPPositivity(C, XUnf, B, A);

if(penalization == 1)
    gradA = gradA + 2 * alpha * A;
    gradB = gradB + 2 * alpha * B;
    gradC = gradC + 2 * alpha * C;
elseif(penalization == 2)
    gradA = gradA + 4 * alpha * A .* A .* A;
    gradB = gradB + 4 * alpha * B .* B .* B;
    gradC = gradC + 4 * alpha * C .* C .* C;
end