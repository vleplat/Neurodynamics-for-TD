function [A B C error] = bfgsP(X, nbFactors, options, AInit, BInit, CInit)

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
% [A B C error] = bfgsP(X, nbFactors, options, AInit, BInit, CInit)
%
% Computes the CP decomposition of a nonnegative third-order tensor. The 
% minimization is based on a quasi-Newton approach, more precisely the 
% Broyden–Fletcher–Goldfarb–Shanno(BFGS) algorithm. It computes an
% approximation of the hessian matrix. The resulting loadings are
% nonnegative.
%
% [A B C error] = bfgsP(X, nbFactors)
% Given a nonnegative third order tensor X and the desired number of 
% factors, it returns the three positive factor matrices and a vector 
% containing the reconstruction error along the iterations.
%
% [A B C error] = bfgsP(X, nbFactors, options) allows to customize some
% parameters in the structure options. Put [] to use default settings. 
%
% [A B C error] = bfgsP(X, nbFactors, options, AInit, BInit, CInit) 
% uses fixed loadings for the initial values. Note that the initial
% loadings can contain negative entries, and it is not a problem or a bad 
% choice for this nonnegative algorithm.
% However the final loadings WILL be positive.
%
% See also gradP, cgP, createOptions.

if(nargin < 3 || isempty(options))
    options = createOptions(0, 0, 0, 550);
end

options.algo = 'BFGS';

if(nargin < 4)
    [A B C error] = quasiNewtonP(X, nbFactors, options);
else
    [A B C error] = quasiNewtonP(X, nbFactors, options, AInit, BInit, CInit);
end