function M = updateMatrix(M, mu, direction)

% This function is part of the positive tensor package
% June, 2011
% J-P. Royer, P. Comon, N. Thirion-Moreau
%
% This package was developped according to the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors
%
% -----------------------------------
% M = updateMatrix(M, mu, direction) 
% 
% Updates the matrix M using a scalar stepsize mu, and a 
% vectorized direction.

M = M + mu * direction;