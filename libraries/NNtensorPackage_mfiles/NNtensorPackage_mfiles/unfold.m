function XUnf = unfold(X, dimOrder)

% This function is part of the positive tensor package
% June, 2011
% J-P. Royer, P. Comon, N. Thirion-Moreau
%
% This package was developped according to the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors
%
% ----------------------------------
% XUnf = unfold(X, order)
%
% Converts the 3rd order tensor X into a matrix, using dimOrder vector, 
% with values [i j k], 1 <= i, j, k <= 3 , et i != j != k.
%
% XUnf is a size(X, i) x ( size(X, j) * size(X, k) ) matrix.
% 
% Ex: 
%   X = rand(5, 6, 8);
%   XUnf = unfold(X, [3 1 2]); % XUnf is a 8 * 30 matrix 

X2 = permute(X, dimOrder);
XUnf = reshape(X2, size(X2, 1), size(X2, 2) * size(X2, 3)); 