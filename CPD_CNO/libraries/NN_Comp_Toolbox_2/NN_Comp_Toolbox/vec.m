function a = vec(G)
%-------------------------------------------------------------------------%
% [a]     =     vec(G)
%
% Vectorize a tensor G as in Kolda, 2008
%
% The inputs are:
%
% - G           : tensor to be vectorized
%
% The outputs are:
% 
% - a           : vectorized tensor
%
% List of updates                 -     06/2014     -     Jeremy Cohen
%                                       Creation of the file
%                                 -     07/01/2015  -     Jeremy Cohen
%                                       Reshaping
%-------------------------------------------------------------------------%


[K,L,M]   =     size(G);
a   =     reshape(permute(G,[3,2,1]),K*L*M,1,1);

end