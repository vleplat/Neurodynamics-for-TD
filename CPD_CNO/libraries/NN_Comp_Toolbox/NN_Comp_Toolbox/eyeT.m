function [ Y ] = eyeT ( n,m )
%-------------------------------------------------------------------------%
% [Y]     =     eyeT(n,m)
%
% Superdiagonal tensor with entries 1 or 0
%
% The inputs are:
%
% - n           : rank of Y
% - m           : size of Y
%
% The outputs are:
% 
% - Y           : Superdiagonal tensor 
%
% List of updates                 -     06/2014     -     Jeremy Cohen
%                                       Creation of the file
%                                 -     07/01/2015  -     Jeremy Cohen
%                                       Reshaping
%-------------------------------------------------------------------------%

Y    =     zeros(m,m,m);

for i=1:n

     Y(i,i,i)=1;

end

end

