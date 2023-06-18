function [A B C] = initializeLoadings(X, nbFactors)

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
% [A B C] = initializeLoadings(X, nbFactors)
%
% Random initialization of the loadings. X is the 3rd order tensor, 
% nbFactors is the desired tensor rank.

A = randn(size(X, 1), nbFactors);
B = randn(size(X, 2), nbFactors);
C = randn(size(X, 3), nbFactors);