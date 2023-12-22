function [A B C error] = quasiNewtonP(X, nbFactors, options, AInit, BInit, CInit)

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
% [A B C error] = quasiNewtonP(X, nbFactors, options, AInit, BInit, CInit)
%
% Computes the CP decomposition of a nonnegative third-order tensor. The 
% minimization is based on quasi-Newton algorithms.
%
% [A B C error] = quasiNewtonP(X, nbFactors)
% Given a nonnegative third order tensor X and the desired number of 
% factors, it returns the three factor matrices and a vector containing 
% the reconstruction error along the iterations.
%
% [A B C error] = quasiNewtonP(X, nbFactors, options) allows to customize some
% parameters in the structure options. Put [] to use default settings.
%
% [A B C error] = quasiNewtonP(X, nbFactors, options, AInit, BInit, CInit) 
% uses fixed loadings for the initial values.
%
% See also bfgsP, createOptions.

errIt = 1000;
crit = 10^-16;
diff = 1;
error = zeros(1, options.maxIt);
XDep = unfold(X, [1 2 3]);
H = eye(sum(size(X)) * nbFactors);
isBFGS = 0;
myEpsilon = 1e-10;

if(strcmp(options.algo, 'BFGS'))
    isBFGS = 1;
end

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

A(find(A == 0)) = myEpsilon;
B(find(B == 0)) = myEpsilon;
C(find(C == 0)) = myEpsilon;

[A, B, C] = normalize(A, B, C); %%Might be commented
M = [vec(A); vec(B); vec(C)];
it = 0;
mu = 0;

%We compute the gradients and stack them into a large vector
[gradA, gradB, gradC] = computeGradientsPositivity(X, A, B, C, options.penalization, options.alpha);
gradient = [gradA(:); gradB(:); gradC(:)];

%Go
while (diff > crit && it < options.maxIt && errIt > 2 * eps)
    it=it+1;

    %The new descent direction
    if(strcmp(options.algo, 'LM'))
        direction = -inv(H + 0.001) * gradient;
    else
        direction = -H * gradient;
    end

    %Computation of mu
    if(mod(it, 10) == 0 || it == 1 || options.backtracking == 0)
        mu = computeOptimalStepsize(X, M, direction, nbFactors, options.penalization, options.alpha);
    else
        mu=rechlin(X, M, gradient, direction, nbFactors, options.penalization, options.alpha);
    end

    %Update step
    M = updateMatrix(M, mu, direction);
    
    %Normalization step
    [A, B, C] = extractLoadings(M, size(X), nbFactors);
    [A, B, C] = normalize(A, B, C); %%Might be commented
    M = [vec(A); vec(B); vec(C)];

    [gradA, gradB, gradC] = computeGradientsPositivity(X, A, B, C, options.penalization, options.alpha);% monEpsilon);
    oldGradient = gradient;
    gradient = [vec(gradA); vec(gradB); vec(gradC)];

    %if the gradient is nearly a matrix of zeros, we stop
    if (abs(mean(mean(gradient))) < eps)
        break;
    end

    deltaM = mu * direction;
    deltaGradient = gradient - oldGradient;

    %Hessian matrix
    if(isBFGS)
        %BFGS H-1
        rho = 1 / (deltaGradient.' * deltaM);
        H = H + rho * (1 + rho * deltaGradient.' * H * deltaGradient) * deltaM * deltaM.' ...
            - rho * deltaM * deltaGradient.' * H - rho * H * deltaGradient * deltaM.';
    elseif(strcmp(options.algo, 'DFP'))
        %DFP H-1
        H = H + (deltaM * deltaM.') / (deltaM.' * deltaGradient)...
            - (H * deltaGradient * deltaGradient.' * H) / (deltaGradient.' * H * deltaGradient);
    else
        %BFGS H ?
        H = H + (deltaGradient * deltaGradient.') / trace(deltaGradient.' * deltaM)...
            - ((H * deltaM) * (H * deltaM).') /  trace((H * deltaM).' * deltaM); 
    end

    %We evaluate the new reconstruction error and we display it
    oldError = errIt;
    errIt = computeError(XDep, A .* A, B .* B, C .* C, options.penalization, options.alpha);
    error(it) = errIt;

    diff = abs((errIt - oldError) / oldError);

    if mod(it, 10) == 0
        %fprintf('erreur : %d\tit : %d\n', erreur(it), it);
        fprintf('error : %12.10f  \tit : %d\n', error(it), it);
    end
end

error = error(1 : it);

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