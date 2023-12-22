function  [A,B,C,X,Y] = tensor_gen(dim,R,sigma)
%-------------------------------------------------------------------------%
% [A,B,C,X,Y] = tensor_gen(dim,R,sigma)
%
% Random CP model generation. All factors are nonnegative. Factors A and B
% have columns with unitary L1 norm. The measured 3-way block is corrupted
% by additive i.i.d. Gaussian noise.
%
% The inputs are:
% 
% - dim         : tensor dimensions.
% - R           : rank.
% - sigma       : Gaussian noise standard deviation.
%
% The outputs are:
% 
% - A           : first factor.
% - B           : second factor.
% - C           : third factor.
% - X           : tensor without noise.
% - Y           : tensor with noise.
%
% List of updates                 -     15/05/2014  -     Rodrigo Cabral
%                                       Creation of the file 
%-------------------------------------------------------------------------%

%--------------------------Tensor generation------------------------------%
% Unnormalized nonnegative factors
A   =     abs(randn(dim(1),R));
B   =     abs(randn(dim(2),R));
C   =     abs(randn(dim(3),R));

% Normalized nonnegative factors
%A         =     A./ repmat(sum(A), dim(1), 1);
%B         =     B./ repmat(sum(B), dim(2), 1);

% First mode tensor without noise
X_m1=     A*kr(C,B)';
% Tensor without noise
X   =     reshape(X_m1,dim);
% Tensor measurement
Y   =     X+sigma*randn(dim);
%-------------------------------------------------------------------------%