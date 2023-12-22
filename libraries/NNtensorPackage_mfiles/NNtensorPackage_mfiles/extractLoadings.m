function [A, B, C] = extractLoadings(M, dimCube, nbFactors)

% This function is part of the positive tensor package
% June, 2011
% J-P. Royer, P. Comon, N. Thirion-Moreau
%
% This package was developped according to the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors
%
% ------------------------------------
% [A, B, C] = extractLoadings(M, dimCube, nbFactors)
% 
% Extracts the 3 factor matrices A, B and C.
%
% M is a big vector containing the matrices A, B, and C stacked under 
% a vectorized form
%
% dimCube is the size of the 3rd order tensor we want to decompose.
% dimCube = size(X);
%
% nbFactors is the tensor rank. It corresponds to the number of columns of
% the matrices A, B, and C.

val = 1;
A = reshape(M(val : dimCube(1) * nbFactors, :), dimCube(1), nbFactors);

val = val + dimCube(1) * nbFactors;
B = reshape(M(val : (val + dimCube(2) * nbFactors - 1), :), dimCube(2), nbFactors);

val = val + dimCube(2) * nbFactors;
C = reshape(M(val : end, :), dimCube(3), nbFactors);