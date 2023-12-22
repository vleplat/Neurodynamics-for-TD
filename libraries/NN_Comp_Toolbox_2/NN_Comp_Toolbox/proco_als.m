function [A,B,C,Lambda,T_als,err_als,ert]     =     proco_als(Y,iter,dim_c,A_0,B_0,C_0)
%-------------------------------------------------------------------------%
% [A_c,B_c,C_c,err_als_vec]=     proco_als(Y,iter,dim_c,A_0,B_0,C_0)
%
% Compressed and projected ALS algorithm for large-scale 3-way nonnegative
% tensor decomposition.
%
% The inputs are:
% 
% - Y           : data block.
% - iter        : maximum number of iterates.
% - dim_c       : vector with the dimensions of the compressed model.
% - A_0,B_0,C_0 : initial parameters for the factors.
%
% The outputs are:
% 
% - A           : first factor.
% - B           : second factor.
% - C           : third factor.
% - Lambda      : amplitude of the R components
% - T_als       : reconstructed tensor
% - err_als     : final error
%
% List of updates                 -     28/08/2014  -     Rodrigo Cabral
%                                       Creation of the file
%                                 -     02/10/2014  -     Jeremy Cohen
%                                       Reshaping, modified outputs
%                                 -     07/01/2015  -     Jeremy Cohen
%                                       Added Lambda as output
%-------------------------------------------------------------------------%


%--------------------------Tensor compression-----------------------------%
% Size of the tensor
dim_u     =     size(Y);
R    =    size(A_0,2);
% Index for n-mode matricization
J_u       =     prod(dim_u)./dim_u;
% Data normalization
norm_Y    =     sum(sum(sum(Y.^2)));
Y         =     Y/sqrt(norm_Y);


% Three different unfoldings of the tensor
Y_mode    =     cell(3);
% Data tensor unfolding
Y_mode{1} =     reshape(permute(Y,[1,2,3]),dim_u(1),J_u(1));  
Y_mode{2} =     reshape(permute(Y,[2,1,3]),dim_u(2),J_u(2));
Y_mode{3} =     reshape(permute(Y,[3,1,2]),dim_u(3),J_u(3));

% Singular Value Decompositions (single precision is faster)
% [U,~,~]=svd(single(Y_mode{1}),'econ');
% [V,~,~]=svd(single(Y_mode{2}),'econ');
% [W,~,~]=svd(single(Y_mode{3}),'econ');
% % Low rank approximation
% U=U(:,1:dim_c(1));
% V=V(:,1:dim_c(2));
% W=W(:,1:dim_c(3));

[U,~,~]   =     rsvd(Y_mode{1}, dim_c(1), 1);
[V,~,~]   =     rsvd(Y_mode{2}, dim_c(2), 1);
[W,~,~]   =     rsvd(Y_mode{3}, dim_c(3), 1);

% Index for n-mode matricization of compressed data
J_c  =    prod(dim_c)./dim_c;
% Truncated core tensors (compressed data)
G    =    MultProd(Y,U',V',W');

% Three different unfoldings of the tensor
G_mode    =     cell(3);
% Data tensor unfolding
G_mode{1} =     reshape(permute(G,[1,2,3]),dim_c(1),J_c(1));  
G_mode{2} =     reshape(permute(G,[2,1,3]),dim_c(2),J_c(2));
G_mode{3} =     reshape(permute(G,[3,1,2]),dim_c(3),J_c(3));


%-------------------------------------------------------------------------%

%-----------------------Algorithm initialization--------------------------%

% Storing initial condition on the compressed space
A_c  =    (U')*A_0;
B_c  =    (V')*B_0;
C_c  =    (W')*C_0;

%-------------------------------------------------------------------------%

fprintf('\n\n')
fprintf('\n Compressed Projected ALS - CP ALS \n\n')

%--------------------------Algorithm iterate------------------------------%
t         =     0;

ert=[];
while     t<iter
    
    % Iteration counter
    t=t+1;
    
    % Least squares solution for each factor
    
    % A factor update
    A_c =     G_mode{1}*pinv((kr(C_c,B_c))');
    % Projection to have nonnegative factors on the uncompressed space
    A_u   =     U*A_c;          % Uncompressed estimated factor 
    A_u(A_u<0)  =     0;        % Projection on the nonnegative orthant     
    A_c         =     U'*A_u;   % Return compressed estimated factor
    
    % B factor update
    B_c =     G_mode{2}*pinv((kr(C_c,A_c))');
    % Projection to have nonnegative factors on the uncompressed space
    B_u   =     V*B_c;
    B_u(B_u<0)  =     0;
    B_c         =     V'*B_u;
    
    % C factor update
    C_c =     G_mode{3}*pinv((kr(B_c,A_c))');
    % Projection to have nonnegative factors on the uncompressed space
    C_u   =     W*C_c;
    C_u(C_u<0)  =     0;
    C_c         =     W'*C_u;

    if (mod(t,50) == 0 )%&& t<iter)
    fprintf('\tit : %d\n',t)    
    end
    
    %Normalisation by stacking on C_c
C_c =     C_c.*repmat(sqrt(sum(A_c.^2).*sum(B_c.^2)),dim_c(3),1);
A_c =     A_c.*repmat(1./sqrt(sum(A_c.^2)),dim_c(1),1);
B_c =     B_c.*repmat(1./sqrt(sum(B_c.^2)),dim_c(2),1);
    

% Uncompressed factors
A   =     U*A_c;
B   =     V*B_c;
C   =     W*C_c;

A(A<0)    =     0;
B(B<0)    =     0;
C(C<0)    =     0;

[~,R]=size(A);
Lambda    =     zeros(R,R,R); % initial diagonal core for CP
for r=1:R
   Lambda(r,r,r)=sqrt(norm_Y)*norm(A(:,r))*norm(B(:,r))*norm(C(:,r));
end

A   =     A.*repmat(1./sqrt(sum(A.^2)),dim_u(1),1);
B   =     B.*repmat(1./sqrt(sum(B.^2)),dim_u(2),1);
C   =     C.*repmat(1./sqrt(sum(C.^2)),dim_u(3),1);


T_als     =     MultProd(Lambda,A,B,C);

err_als =     sqrt(sum(sum(sum((T_als/sqrt(norm_Y)-Y).^2))));

ert=[ert,err_als];

end

fprintf('Reconstruction error : %g  ',err_als) 

A(isnan(A))     =     0;
B(isnan(B))     =     0;
C(isnan(C))     =     0;
%-------------------------------------------------------------------------%    