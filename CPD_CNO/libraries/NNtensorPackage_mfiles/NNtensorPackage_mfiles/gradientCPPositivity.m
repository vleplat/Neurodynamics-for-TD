function deltaL = gradientCPPositivity(L, X, L2, L3)

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
% deltaL = gradientCPPositivity(L, X, C, B)
%
% Determines the positive gradient versus L (L can represent either matrix 
% A, B or C in the decomposition).
%
% X is the tensor unfolded into the chosen mode, depending on the gradient 
% matrix we want to estimate.
% 
% L2 and L3 are the two other matrices, corresponding to the two other 
% modes.

katriR = khatriRao(L2 .* L2, L3 .* L3);
deltaL = 4 * (L) .* ((- X + (L .* L) * (katriR).') * (katriR));