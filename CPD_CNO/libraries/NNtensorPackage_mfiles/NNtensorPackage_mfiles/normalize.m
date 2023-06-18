function [A, B, C] = normalize(A, B, C)

% This function is part of the positive tensor package.
% June, 2011.
% J-P. Royer, P. Comon, N. Thirion-Moreau.
%
% This package was developped from the work presented here:
% http://dx.doi.org/10.1016/j.sigpro.2011.03.006
% and entitled: Computing the polyadic decomposition of nonnegative third
% order tensors.
%
% --------------------------------
% [A, B, C] = normalize(A, B, C)
%
% Normalizes the loadings columnwise and returns the result.

fac = size(A, 2);
for jj = 1 : fac
    A(:, jj) = A(:, jj) * norm(B(:, jj)) * norm(C(:, jj));
end
for jj = 1 : fac
    B(:, jj) = B(:, jj) / norm(B(:, jj));
end
for jj = 1 : fac
    C(:, jj) = C(:, jj) / norm(C(:, jj));
end

% for jj = 1 : fac
%     A(:, jj) = A(:, jj) .* A(:, jj) * norm(B(:, jj) .* B(:, jj)) * norm(C(:, jj) .* C(:, jj));
% end
% for jj = 1 : fac
%     B(:, jj) = (B(:, jj) .* B(:, jj)) / norm(B(:, jj) .* B(:, jj));
% end
% for jj = 1 : fac
%     C(:, jj) = (C(:, jj) .* C(:, jj)) / norm(C(:, jj) .* C(:, jj));
% end
% 
% A = sqrt(A);
% B = sqrt(B);
% C = sqrt(C);