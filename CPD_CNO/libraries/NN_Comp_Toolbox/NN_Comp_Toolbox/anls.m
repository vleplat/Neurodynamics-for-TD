function [A_aux,B_aux,C_aux,T_als,err_als]=     anls(Y,iter,A0,B0,C0)
%-------------------------------------------------------------------------%
% [A_aux,B_aux,C_aux,err_als_vec]=     als_proj_1(Y,iter,A0,B0,C0)
%
% Alternating least squares algorithm for 3-way nonnegative tensor decompo-
% sition. The result of each iterate is projected on the nonnegative
% orthant.
%
% The inputs are:
% 
% - Y           : data block.
% - iter        : maximum number of iterates.
% - A_0,B_0,C_0 : initial parameters for the factors.
%
% The outputs are:
% 
% - A_aux       : first factor.
% - B_aux       : second factor.
% - C_aux       : third factor.
% - err_als_vec : vector of errors.
%
% List of updates                 -     19/06/2014  -     Rodrigo Cabral
%                                       Creation of the file 
%                                 -     02/10/2014  -     Jeremy Cohen
%                                       Reshaping, modified outputs
%-------------------------------------------------------------------------%


%-----------------------Algorithm initialization--------------------------%
% Size of the tensor
N         =     size(Y);
% Index for n-mode matricization
J         =     prod(N)./N;
% Data normalization
norm_Y    =     sum(sum(sum(Y.^2)));
Y         =     Y/sqrt(norm_Y);
% Three different unfoldings of the tensor
Y_mode    =     cell(3);
% Data tensor unfolding
Y_mode{1} =     reshape(permute(Y,[1,2,3]),N(1),J(1));  
Y_mode{2} =     reshape(permute(Y,[2,1,3]),N(2),J(2));
Y_mode{3} =     reshape(permute(Y,[3,1,2]),N(3),J(3));
% Storing auxiliary variables
A_aux     =     A0;
B_aux     =     B0;
C_aux     =     C0;

%-------------------------------------------------------------------------%

fprintf('\n\n')
fprintf('\n Projected ALS - ANLS \n\n')

%--------------------------Algorithm iterate------------------------------%
t         =     0;

while      t<iter  
         
    % Iteration counter
    t=t+1;

    % Least squares solution for each factor and projection
    
    % A factor update
    A_aux       =     Y_mode{1}*pinv((kr(C_aux,B_aux))');
    % Projection on the nonnegative orthant
    A_aux(A_aux<0)    =     0;
    
    % B factor update
    B_aux       =     Y_mode{2}*pinv((kr(C_aux,A_aux))');
    % Projection on the nonnegative orthant
    B_aux(B_aux<0)    =     0;
    
    % C factor update
    C_aux             =     Y_mode{3}*pinv((kr(B_aux,A_aux))');
    % Projection on the nonnegative orthant
    C_aux(C_aux<0)    =     0;


    if (mod(t,50) == 0)
    % Printing the iterate error and step size
    fprintf('\tit : %d\n',t)
    end

    %Normalisation by stacking on C_c
C_aux =     C_aux.*repmat(sqrt(sum(A_aux.^2).*sum(B_aux.^2)),N(3),1);
A_aux =     A_aux.*repmat(1./sqrt(sum(A_aux.^2)),N(1),1);
B_aux =     B_aux.*repmat(1./sqrt(sum(B_aux.^2)),N(2),1);
    
end


[~,R]=size(A_aux);
Lambda    =     zeros(R,R,R); % initial diagonal core for CP
for r=1:R
    Lambda(r,r,r)=sqrt(norm_Y)*norm(A_aux(:,r))*norm(B_aux(:,r))*norm(C_aux(:,r));
end

A_aux   =     A_aux.*repmat(1./sqrt(sum(A_aux.^2)),N(1),1);
B_aux   =     B_aux.*repmat(1./sqrt(sum(B_aux.^2)),N(2),1);
C_aux   =     C_aux.*repmat(1./sqrt(sum(C_aux.^2)),N(3),1);

T_als     =     MultProd(Lambda,A_aux,B_aux,C_aux);

% Squared norm of the error

err_als =     sqrt(sum(sum(sum((T_als/sqrt(norm_Y)-Y).^2))));


fprintf('Reconstruction error : %g  ',err_als)

% Dealing with NaN
A_aux(isnan(A_aux))   =     0;
B_aux(isnan(B_aux))   =     0;
C_aux(isnan(C_aux))   =     0;
%-------------------------------------------------------------------------%    