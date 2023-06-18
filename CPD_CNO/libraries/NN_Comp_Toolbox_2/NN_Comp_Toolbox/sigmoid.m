function [ N ] = sigmoid( char,alpha,M )
%-------------------------------------------------------------------------%
% [ N ] = sigmoid( char,alpha,M )
%
% Computes coefficient-wise the designated function to matrix M
%  
%
% The inputs are:
% 
% - char        : penalisation function. Choices are 
%      + 'exp'  : exponential barrier. Recommanded.
%      + 'tanh' : hyperbolic tangent penalisation.
%      + 'atan' : arctangent penalisation.
%      + 'sigm' : sigmoid penalisation.
%      + 'log'  : logarithmic barrier.
%
% - alpha       : vector with stiffness constraints parameters for each mode
% - M           : input matrix
%
% The outputs are:
% 
% - N           : output matrix 'char'(M)
%
% List of updates                 -        08/2014  -    Jérémy Cohen
%                                       Creation of the file
%                                 -        10/2014  -    Jeremy Cohen
%                                       Reshaping, modified outputs
%-------------------------------------------------------------------------%


switch(char)
    
    case('tanh')
        
          N     =     (1-tanh(alpha*M))/2;
        
    case('atan')
        
          N     =     1/2-atan(alpha*M)/pi;
        
    case('sigm')
        
          N     =     1-1./(1+exp(-alpha*M));
        
    case('log')
        
          N     =     -1/alpha*log(M);
        
    case('exp')
        
          N     =     exp(-alpha*M);
        
end


end

