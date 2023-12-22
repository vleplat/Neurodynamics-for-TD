function [A B C error] = gradP(X, nbFactors, options, AInit, BInit, CInit)

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
% [A B C error] = gradP(X, nbFactors, options, AInit, BInit, CInit)
%
% Computes the CP decomposition of a nonnegative third-order tensor. The 
% minimization is based on a gradient algorithm, and the resulting loadings
% are nonnegative.
%
% [A B C error] = gradP(X, nbFactors) 
% Given a nonnegative third order tensor X and the desired number of 
% factors, it returns the three positive factor matrices and a vector 
% containing the reconstruction error along the iterations.
%
% [A B C error] = gradP(X, nbFactors, options) allows to customize some
% parameters in the structure options. Put [] to use default settings. 
%
% [A B C error] = gradP(X, nbFactors, options, AInit, BInit, CInit) uses
% fixed loadings for the initial values. Note that the initial loadings 
% can contain negative entries, and it is not a problem or a bad choice for 
% this nonnegative algorithm.
% However the final loadings WILL be positive.
%
% See also cgP, bfgsP, createOptions.

if(nargin < 3 || isempty(options))
    options = createOptions(0, 0, 0, 550);
end

newErr = 1000;
crit = 10^-16;
diffNorm = 1;
error = zeros(1, options.maxIt);
XDep = unfold(X, [1 2 3]);

if(nargin < 4)
    [A B C] = initializeLoadings(X, nbFactors);
    
    disp('Random initialization');
else    
    A = AInit;
    B = BInit;
    C = CInit;
    
    disp('Fixed initialisation');
end

displayOptions(options);

[A, B, C] = normalize(A, B, C); %Might be commented
M = [vec(A); vec(B); vec(C)];
itTot = 0;

while (diffNorm > crit && itTot < options.maxIt && newErr > 2 * eps)
    itTot = itTot + 1;

    [gradA, gradB, gradC] = computeGradientsPositivity(X, A, B, C, options.penalization, options.alpha);
    gradA(find(abs(gradA) < eps)) = eps; gradB(find(abs(gradB) < eps)) = eps; gradC(find(abs(gradC) < eps)) = eps;
    gradient = [vec(gradA); vec(gradB); vec(gradC)];
    direction = - gradient;

    if(mod(itTot, 10) == 0 || itTot == 1 || options.backtracking == 0)
        mu = computeOptimalStepsize(X, M, direction, nbFactors, options.penalization, options.alpha);
    else
        mu=rechlin(X, M, gradient, direction, nbFactors, options.penalization, options.alpha);
    end

    M = updateMatrix(M, mu, direction);
    [A, B, C] = extractLoadings(M, size(X), nbFactors);
    [A, B, C] = normalize(A, B, C); %Might be commented
    M = [vec(A); vec(B); vec(C)];

    oldErr = newErr;
    newErr = computeError(XDep, A .* A, B .* B, C .* C, options.penalization, options.alpha);
    error(itTot) = newErr;
    diffNorm = abs((newErr - oldErr) / newErr);

    if mod(itTot, 10) == 0
        fprintf('errorr : %12.10f  \tit : %d\n', error(itTot), itTot);
    end
end

error = error(1 : itTot);

A = A .* A;
B = B .* B;
C = C .* C;


function displayOptions(options)

if(options.penalization == 0)
    disp('No penalization');
else
    fprintf('Penalization L%d with alpha = %d\n', options.penalization, options.alpha);
end

if(options.backtracking == 1)
    disp('Backtracking alternated with ELS');
else
    disp('No backtracking');
end

fprintf('%d iterations max\n', options.maxIt);

disp(' ')