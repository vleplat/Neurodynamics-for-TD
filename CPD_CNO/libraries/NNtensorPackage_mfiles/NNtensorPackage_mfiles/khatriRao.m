function K = khatriRao(A, B)

% K = khatriRao(A, B)
% 
% Computes the khatri-Rao product between two matrices with the same number
% of columns.

if(size(A, 2) ~= size(B, 2))
    error('A and B must have the same number of columns');
end

f = size(A, 2);
tA = size(A, 1);
tB = size(B, 1);
K = zeros(tA * tB, f);

for ii = 1 : f
    column = B(:, ii) * A(:, ii).';
    K(:, ii) = column(:);
end

% faster in some cases...
% for jj = 1 : f
%     for ii = 1 : tA
%     %K(:, f) = A(1 : end, f) .* B(:, f);
%         K((ii - 1) * tB + 1 : ii * tB, jj) = A(ii, jj) * B(:, jj);
%     end
% end