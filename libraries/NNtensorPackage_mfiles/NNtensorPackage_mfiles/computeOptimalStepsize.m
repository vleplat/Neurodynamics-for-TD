function mu = computeOptimalStepsize(X, M, direction, nbFactors, penalization, alpha)

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
% mu = computeOptimalStepsize(X, M, direction, nbFactors, penalization,
% alpha)
%
% Computes the optimal stepsize, minimizing the polynomial with respect to
% mu.
%
% X is the 3rd order tensor.
% M are the loadings stacked and vectorized.
% nbFactors is the desired tensor rank.
% penalization is 1 or 2 if applied to the problem (L1 or L2 penalty), 
% 0 if not.
% alpha is a small scalar defining the penalty coefficient (if applied).

dimCube = size(X);

[A, B, C] = extractLoadings(M, dimCube, nbFactors);
[directionA, directionB, directionC] = extractLoadings(direction, dimCube, nbFactors);

XUnf = unfold(X, [1 2 3]);

b0 = A .* A;
b1 = 2 * A .* directionA;
b2 = directionA .* directionA;

a4 = khatriRao(directionC .* directionC, directionB .* directionB);

a3 = 2 * khatriRao(C .* directionC, directionB .* directionB) ...
    + 2 * khatriRao(directionC .* directionC, B .* directionB);

a2 = 4 * khatriRao(C .* directionC, B .* directionB) ...
    + khatriRao(directionC .* directionC, B .* B) ...
    + khatriRao(C .* C, directionB .* directionB);

a1 = + 2 * khatriRao(C .* directionC, B .* B) ...
    + 2 * khatriRao(C .* C, B .* directionB);

a0 = khatriRao(C .* C, B .* B);

K6 = - b2 * a4.';
K5 = - b1 * a4.' - b2 * a3.';
K4 = - b0 * a4.' - b1 * a3.';
K3 = - b0 * a3.' - b1 * a2.' - b2 * a1.';
K2 = - b0 * a2.' - b1 * a1.' - b2 * a0.';
K1 = - b0 * a1.' - b1 * a0.';
K0 = XUnf - b0 * a0.';

% tic
% c12 = calculerCoeff(K6, K6');
% c11 = calculerCoeff(2*K6, K5');
% c10 = calculerCoeff(2*K6, K4') + calculerCoeff(K5, K5');
% c9 = calculerCoeff(2*K6, K3') + calculerCoeff(2*K5, K4');
% c8 = calculerCoeff(2*K6, K2') + calculerCoeff(2*K5, K3') + calculerCoeff(K4, K4');
% c7 = calculerCoeff(2*K6, K1') + calculerCoeff(2*K5, K2') + calculerCoeff(2*K4, K3');
% c6 = calculerCoeff(2*K6, K0') + calculerCoeff(2*K5, K1') + calculerCoeff(2*K4, K2') + calculerCoeff(K3, K3');
% c5 = calculerCoeff(2*K5, K0') + calculerCoeff(2*K4, K1') + calculerCoeff(2*K3, K2');
% c4 = calculerCoeff(2*K4, K0') + calculerCoeff(2*K3, K1') + calculerCoeff(K2, K2');
% c3 = calculerCoeff(2*K3, K0') + calculerCoeff(2*K2, K1');
% c2 = calculerCoeff(2*K2, K0') + calculerCoeff(K1, K1');
% c1 = calculerCoeff(2*K1, K0');
% c0 = calculerCoeff(K0, K0');
% toc

c12 = trace(K6 * K6');
c11 = trace(2 * K6 * K5');
c10 = trace(2 * K6 * K4' + K5 * K5');
c9 = trace(2 * (K6 * K3' + K5 * K4'));
c8 = trace(2 * (K6 * K2' + K5 * K3') + K4 * K4');
c7 = trace(2 * (K6 * K1' + K5 * K2' + K4 * K3'));
c6 = trace(2 * (K6 * K0' + K5 * K1' + K4 * K2') + K3 * K3');
c5 = trace(2 * (K5 * K0' + K4 * K1' + K3 * K2'));
c4 = trace(2 * (K4 * K0' + K3 * K1') + K2 * K2');
c3 = trace(2 * (K3 * K0' + K2 * K1'));
c2 = trace(2 * K2 * K0' + K1 * K1');
c1 = trace(2 * K1 * K0');
c0 = trace(K0 * K0');

bB0 = B .* B;
bB1 = 2 * B .* directionB;
bB2 = directionB .* directionB;

bC0 = C .* C;
bC1 = 2 * C .* directionC;
bC2 = directionC .* directionC;

if(penalization == 1)
    c2 = c2 + alpha * (sum(sum(b2)) + sum(sum(bB2)) + sum(sum(bC2)));
    c1 = c1 + alpha * (sum(sum(b1)) + sum(sum(bB1)) + sum(sum(bC1)));
    c0 = c0 + alpha * (sum(sum(b0)) + sum(sum(bB0)) + sum(sum(bC0)));
elseif(penalization == 2)
    c4 = c4 + alpha * (trace(b2.' * b2) + trace(bB2' * bB2) + trace(bC2' * bC2));
    c3 = c3 + alpha * (trace(2 * b2' * b1) + trace(2 * bB2' * bB1) + trace(2 * bC2' * bC1));
    c2 = c2 + alpha * (trace(2 * b2' * b0 + b1' * b1) + trace(2 * bB2' * bB0 + bB1' * bB1) + trace(2 * bC2' * bC0 + bC1' * bC1));
    c1 = c1 + alpha * (trace(2 * b1' * b0) + trace(2 * bB1' * bB0) + trace(2 * bC1' * bC0)); 
    c0 = c0 + alpha * (trace(b0' * b0) + trace(bB0' * bB0) + trace(bC0' * bC0)); 
end

c11 = c11 / c12;
c10 = c10 / c12;
c9 = c9 / c12;
c8 = c8 / c12;
c7 = c7 / c12;
c6 = c6 / c12;
c5 = c5 / c12;
c4 = c4 / c12;
c3 = c3 / c12;
c2 = c2 / c12;
c1 = c1 / c12;
c0 = c0 / c12;
c12 = 1;

coeffs = [ 12 * c12, 11 * c11, 10 * c10, 9 * c9, 8 * c8, 7 * c7, 6 * c6, 5 * c5, 4 * c4, 3 * c3, 2 * c2, c1];

racines = real(roots(coeffs));

racines = racines(find(racines >= 0));
[mini ind] = min(polyval([c12 c11 c10 c9 c8 c7 c6 c5 c4 c3 c2 c1 c0], racines));
mu = racines(ind);

if (isempty(mu) == 1 || mu == 0)%size(mu, 1) == 0)
    mu = eps;
end

function coeff = calculerCoeff(varargin)

tailleCol = size(varargin{2}, 2);
res = zeros(1, tailleCol);

for ii = 1 : nargin / 2
    A = varargin{ii};
    B = varargin{ii + 1};
    for jj = 1 : tailleCol   
        res(jj) = res(jj) + A(jj, :) *  B(:, jj);
    end
end

clear A B;
coeff = sum(res);
        
        