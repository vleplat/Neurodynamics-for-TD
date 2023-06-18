function G = unvec(a,K,L,M)
%-------------------------------------------------------------------------%
% [G]     =     unvec(a,K,L,M)
%
% Transforms a vector a to a KxLxM tensor
%
% The inputs are:
%
% - a           : initial vector of size KLM
% - K,L,M       : dimensions of the built tensor
%
% The outputs are:
% 
% - G           : unvectorized tensor
%
% List of updates                 -     06/2014     -     Jeremy Cohen
%                                       Creation of the file
%                                 -     07/01/2015  -     Jeremy Cohen
%                                       Reshaping
%-------------------------------------------------------------------------%

G   =     permute(reshape(a,M,L,K),[3,2,1]);

end