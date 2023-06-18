function options = createOptions(penalization, alpha, backtracking, maxIt)

% This function is part of the positive tensor package
% June, 2011
% J-P. Royer, P. Comon, N. Thirion-Moreau
%
% This package was developped according to the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors.
%
% --------------------------------
% options = createOptions(penalization, alpha, backtracking, maxIt) 
% 
% Creates and returns the structure containing the options: the desired 
% penalization (0, 1, 2 for no penalization, L1 penalization, L2 
% penalization), the alpha penalty term, the use of backtracking (0 or 1), 
% the number of iterations.

options = struct('penalization', penalization, 'alpha', alpha, 'backtracking', backtracking, 'maxIt', maxIt);