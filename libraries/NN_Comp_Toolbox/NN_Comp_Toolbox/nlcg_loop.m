function [ theta,Dir,s,Gout,eps] = nlcg_loop( theta,char,alpha,gamma,Dir_prev ,s_prev, Eps, G, Gin, R, U,V,W, b)
%-------------------------------------------------------------------------%
% [ theta,Dir,s,Gout,eps] = nlcg_loop( theta,char,alpha,gamma,Dir_prev ,s_prev, Eps, G, Gin, R, U,V,W, b)
%
% Executes one iteration of the non-linear conjugate gradient
% using Polak-Ribiere direction. 
%
% The inputs are:
% 
% - theta       : parameter vector
% - char        : penalisation function. See 'sigmoid' for various choices,
% - alpha       : vector with stiffness constraints parameters for each mode
% - gamma       : penalisation magnitude
% - Dir_prev    : descent direction of the previous iteration
% - s_prev      : conjugated descent direction of the previous iteration
% - Eps         : table containing current step size for...
% ...the descent and the minimum and maximum stepsize.
% - G           : target compressed tensor
% - Gin         : previous estimate of the target tensor
% - R           : rank of the decomposition
% - U,V,W       : decompression operators
% - b           : increase factor for stepsize in Armijo

% The outputs are:
% 
% - theta       : parameter vector
% - Dir         : descent direction 
% - s           : conjugated descent direction
% - Gout        : current estimate of the target tensor
% - eps         : current stepsize
%
% List of updates                 -        08/2014  -    Jérémy Cohen
%                                       Creation of the file
%                                 -        10/2014  -    Jeremy Cohen
%                                       Reshaping, modified outputs
%                                 -     21/01/2015  -    Jeremy Cohen
%                                       Added zero-barrier of beta
%-------------------------------------------------------------------------%



% RECOVERING THE PARAMETERS

% Dimensions of the Initial and Compressed Spaces

[K,R1]    =     size(U);
[L,R2]    =     size(V);
[M,R3]    =     size(W);

% Stepsizes

eps       =     Eps(1);
mini      =     Eps(2);
maxi      =     Eps(3);

% Penalisation Stiffnesses

alpha1    =     alpha(1);
alpha2    =     alpha(2);
alpha3    =     alpha(3);

% Factor matrices A B C

[A,B,C]   =     ParUnfold(theta,R1,R2,R3,R);
tG        =     unvec(G,R1,R2,R3);


% GRADIENT COMPUTATION

% CP Gradient

KRaoCB    =     kr(C,B);
KRaoCA    =     kr(C,A);
KRaoBA    =     kr(B,A);

GradA     =     (-reshape(tG,R1,R2*R3) + A*KRaoCB') * KRaoCB; 
GradB     =     (-reshape(permute(tG,[2,1,3]),R2,R1*R3) + B*KRaoCA') * KRaoCA;
GradC     =     (-reshape(permute(tG,[3,1,2]),R3,R1*R2) + C*KRaoBA') * KRaoBA;

Grad=ParFold(GradA,GradB,GradC);

% Steepest descent direction

Dir       =     - Grad - gamma/(R*(K+L+M))*sigmoid_diff(char,alpha,U,V,W,A,B,C);

% Conjugate vector computation (Polak-Ribière)

beta      =     max(Dir'*(Dir-Dir_prev)/(Dir_prev'*Dir_prev),0); 

s   =     Dir+beta*s_prev;

%% Armijo method for stepsize selection

% Conjugate gradient descent

theta2      =     theta+eps*s;
    
% New approximation computation

[A0,B0,C0]  =     ParUnfold(theta2,R1,R2,R3,R);
Gout  =     vec(MultProd(eyeT(R,R),A0,B0,C0));
    
% Cost function computation
    
S     =     sum(sum((G-Gin).^2))+gamma/(R*(K+L+M))*(sum(sum(sigmoid(char,alpha1,U*A)))+sum(sum(sigmoid(char,alpha2,V*B)))+sum(sum(sigmoid(char,alpha3,W*C))));
S_2   =     sum(sum((G-Gout).^2))+gamma/(R*(K+L+M))*(sum(sum(sigmoid(char,alpha1,U*A0)))+sum(sum(sigmoid(char,alpha2,V*B0)))+sum(sum(sigmoid(char,alpha3,W*C0))));
    
% Randomized Armijo condition
  
m     =     0;
nu    =     1/2;
    

    while (S-S_2 < nu*(s*eps)'*Dir) && eps>=mini
       
    eps   =    eps*b;
    m     =    m+1;

    % repetition of the previous computations

    theta2     =    theta+eps*s;
    [A0,B0,C0] =    ParUnfold(theta2,R1,R2,R3,R);
    Gout  =    vec(MultProd(eyeT(R,R),A0,B0,C0)); 
    S_2   =    sum(sum((G-Gout).^2))+gamma/(R*(K+L+M))*(sum(sum(sigmoid(char,alpha1,U*A0)))+sum(sum(sigmoid(char,alpha2,V*B0)))+sum(sum(sigmoid(char,alpha3,W*C0))));
    
    end
 
% Increase the stepsize if Abijur immediatly verified   
    
    if m==0 && (S-S_2 >= nu*(s*eps)'*Dir)
        eps     =     min(eps*10,maxi);
    end 

 
    
theta    =     theta2;

end

