function [ Grad ] = sigmoid_diff( char,alpha,U,V,W,A,B,C )
%-------------------------------------------------------------------------%
% [ N ] = sigmoid( char,alpha,M )
%
% Compute the gradient of the L1 norm of the penality function with
% stiffness alpha, applied on M element-wise.
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
% - U,V,W       : decompression operators
% - A,B,C       : compressed factor
%
% The outputs are:
% 
% - Grad        : penalization gradient
%
% Note : The unknowns order is as follows : theta=[A.1, B.1, C.1, A.2, ..., C.R]
%
%
% List of updates                 -        08/2014  -    Jérémy Cohen
%                                       Creation of the file
%                                 -        10/2014  -    Jeremy Cohen
%                                       Added exponential barrier
%-------------------------------------------------------------------------%


% Computation of the penalized elements
M1=U*A;
M2=V*B;
M3=W*C;


% Computation of the gradient

switch(char)
    case('tanh')
        
        M1=U'*alpha(1)*-1/2*(1-tanh(alpha(1)*M1).^2);
        M2=V'*alpha(2)*-1/2*(1-tanh(alpha(2)*M2).^2);
        M3=W'*alpha(3)*-1/2*(1-tanh(alpha(3)*M3).^2);
        
        Grad=ParFold(M1,M2,M3);
        
    case('atan') 
        
        M1=U'*alpha(1)*-2/pi*(1./(1+(alpha(1)*M1).^2));
        M2=V'*alpha(2)*-2/pi*(1./(1+(alpha(2)*M2).^2));
        M3=W'*alpha(3)*-2/pi*(1./(1+(alpha(3)*M3).^2));
        
        Grad=ParFold(M1,M2,M3);       
        
        
    case('sigm')

        M1=U'*alpha(1)*-(exp(-alpha(1)*M1)./(1+exp(-alpha(1)*M1)).^2);
        M2=V'*alpha(2)*-(exp(-alpha(2)*M2)./(1+exp(-alpha(2)*M2)).^2);
        M3=W'*alpha(3)*-(exp(-alpha(3)*M3)./(1+exp(-alpha(3)*M3)).^2);
        
        Grad=ParFold(M1,M2,M3); 
        
    case('log')
        
        M1=U'*-1/alpha(1)*(1./M1);
        M2=V'*-1/alpha(2)*(1./M2);
        M3=W'*-1/alpha(3)*(1./M3);
        
        Grad=ParFold(M1,M2,M3);
        
    case('exp')
        
        M1=U'*-alpha(1)*exp(-alpha(1)*M1);
        M2=V'*-alpha(2)*exp(-alpha(2)*M2);
        M3=W'*-alpha(3)*exp(-alpha(3)*M3);
        
        Grad=ParFold(M1,M2,M3);
        
end
end

