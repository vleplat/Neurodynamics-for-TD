function [ A_ccg, B_ccg, C_ccg, Lambda, T_ccg, err_ccg,ert] = ccg( T,itermax,dim_c,A0,B0,C0,alpha,gamma,Re,char )
%-------------------------------------------------------------------------%
% [ A_ccg, B_ccg, C_ccg, Lambda, T_ccg, err_ccg] = ccg( T,itermax,dim_c,A0,B0,C0,alpha,gamma,Re,char )
%
% Conjugate compressed penalized non-linear gradient for Tensor
% factorisation.
%
% The inputs are:
% 
% - T           : Tensor data.
% - itermax     : maximum number of iterations.
% - dim_c       : vector with the dimensions of the compressed model.
% - A_0,B_0,C_0 : initial parameters for the factors.
% - alpha       : vector with stiffness constraints parameters for each mode
% - gamma       : penalisation magnitude
% - Re          : rank of the decomposition
% - char        : penalisation function. See 'sigmoid' for various choices,
% 'exp' is recommanded with the present version.
%
% The outputs are:
% 
% - A_ccg       : first factor.
% - B_ccg       : second factor.
% - C_ccg       : third factor.
% - Lambda      : amplitude of the R components
% - T_ccg       : reconstructed tensor
% - err_ccg     : final error.
%
% List of updates                 -     28/08/2014  -     Rodrigo Cabral
%                                       Creation of the file
%                                 -     02/10/2014  -     Jeremy Cohen
%                                       Reshaping, modified outputs
%-------------------------------------------------------------------------%



% Conjugate compressed penalized non-linear gradient for Tensor
% factorisation.
% ------ Inputs ----- %
% T : Tensor data
% itermax : Maximum number of iterations
% dim_c : Compressed dimensions
% A0,B0,C0 : Initial guesses for the ccg, to be positive in the
% uncompressed space.
% alpha : Vector containing stiffnesses alpha1 alpha2 alpha3.
% gamma : Penalisation magnitude.

dim_u     =    size(T);
J_u       =    prod(dim_u)./dim_u;

normT     =    sum(sum(sum(T.^2)));
T         =    T/sqrt(normT);

fprintf('\n\n')
fprintf('\n Compressed Projected Gradient - CCG \n\n')

%------------------------------------------------------------------------%

% COMPRESSION

% Three different unfoldings of the tensor
T_mode    =     cell(3);
% Data tensor unfolding
T_mode{1} =     reshape(permute(T,[1,2,3]),dim_u(1),J_u(1));  
T_mode{2} =     reshape(permute(T,[2,1,3]),dim_u(2),J_u(2));
T_mode{3} =     reshape(permute(T,[3,1,2]),dim_u(3),J_u(3));

% Fast HOSVD
[U,~,~]   =     rsvd(T_mode{1}, dim_c(1), 1);
[V,~,~]   =     rsvd(T_mode{2}, dim_c(2), 1);
[W,~,~]   =     rsvd(T_mode{3}, dim_c(3), 1);

G   =     MultProd(T,U',V',W');
%[G,U,V,W,HO_err]     =     HOSVD(T,dim_u,dim_c,1);

%Preprocessing for CCG
G   =     vec(G);

%------------------------------------------------------------------------%

% CCG

% Parameters of CP decomposition algorithm
    
    % Loops

tol    =     -40; % minimal stepsize (in log10) to trigger algorithm convergence criterion
err    =     1; % arbitrary initialisation of the error for 1st loop
i      =     0; % loop index initialisation

eps    =     10^(-10); % initial stepsize
mini   =     0;%10^-16; % minimal stepsize
maxi   =     400; % maximal stepsize

	% Penalisation : linear evolution of parameters (set for 'exp')
 
        % Stiffness
        
alpha01      =     alpha(1); 
alpha02      =     alpha(2);
alpha03      =     alpha(3);

init_dur     =     floor(0.1*itermax); % Duration of the constant initial value for alpha
tran_size    =     floor(0.9*itermax); % Duration of the increase of alpha
mult_factor  =     500;                % Final relative value of alpha (TUNE WITH RESPECT TO PENALIZATION FUNCTION)

shape_alpha  =[ones(1,init_dur),...
    1:((mult_factor-1)/(tran_size)...
    ):mult_factor,mult_factor*ones(1,itermax-init_dur-tran_size-1)]; % Temporal evolution of alpha

alpha_vec    =     [ alpha01 * shape_alpha ; alpha02 * shape_alpha ; alpha03 * shape_alpha ];

        % Magnitude

init_dur2    =     floor(1*itermax); % Duration of the constant initial value for gamma
tran_size2   =     floor(0.9*itermax); % Duration of the increase of gamma; IMPORTANT : Set to non-increasing & non-decreasing by default
mult_factor2 =     5; % Final relative value of gamma

shape_gamma  =     [ones(1,init_dur2),...
                   1:((mult_factor2-1)/(tran_size2)):mult_factor2,...
                   mult_factor2*ones(1,itermax-init_dur2-tran_size2-1)]; % Temporal evolution of gamma

gamma_vec    =     gamma*shape_gamma;

% Initialisation

    % Loop variables

s_prev       =     0; % initialisation of the conjugate gradient
Dir_prev     =     ones(Re*(dim_c(1)+dim_c(2)+dim_c(3)),1); % initialisation for the descent direction


    % First guess
    
Ac  =     (U')*A0;
Bc  =     (V')*B0;
Cc  =     (W')*C0;

    % Normalisation L1 of first guessed factors
 
Ac   =     Ac.*repmat(1./sqrt(sum(Ac.^2)),dim_c(1),1);
Bc   =     Bc.*repmat(1./sqrt(sum(Bc.^2)),dim_c(2),1);
Cc   =     Cc.*repmat(1./sqrt(sum(Cc.^2)),dim_c(3),1);


    % Unknowns folding
    
theta     =     ParFold(Ac,Bc,Cc); % Related parameters vector
Gin  =    vec(MultProd(eyeT(Re,Re),Ac,Bc,Cc)); % Vectorized initial core

Lambda    =     zeros(Re,Re,Re); % initial diagonal core for CP

%------------------------------------------------------------------------%

% Iterations

ert=[];
%&& eps>10^(tol)

while i<itermax 
        
    % Main loop
  
    [theta,Dir_prev,s_prev,Gin,eps]     =     nlcg_loop(theta,char,alpha_vec(:,i+1),gamma_vec(i+1),...
                                                        Dir_prev ,s_prev, [eps,mini,maxi], G, Gin, Re, U,V,W,1/2);
    i     =     i+1;
    
if mod(i,1)==0
fprintf('\tit  %d\n',i)
end

    % Determination of the CP factors and core


[A_hat,B_hat,C_hat]   =     ParUnfold(theta,dim_c(1),dim_c(2),dim_c(3),Re);


% Projection for error computation and final output

A_ccg     =     U*A_hat;
A_ccg(A_ccg<0)  =     0;

B_ccg     =     V*B_hat;
B_ccg(B_ccg<0)  =     0;

C_ccg     =     W*C_hat;
C_ccg(C_ccg<0)  =     0;

    % Final Normalisation
for r=1:Re
    Lambda(r,r,r)=sqrt(normT)*norm(A_ccg(:,r))*norm(B_ccg(:,r))*norm(C_ccg(:,r));
end

A_ccg     =     A_ccg.*repmat(1./sqrt(sum(A_ccg.^2)),dim_u(1),1);
B_ccg     =     B_ccg.*repmat(1./sqrt(sum(B_ccg.^2)),dim_u(2),1);
C_ccg     =     C_ccg.*repmat(1./sqrt(sum(C_ccg.^2)),dim_u(3),1);

T_ccg     =     MultProd(Lambda,A_ccg,B_ccg,C_ccg);

% Error and computation time

err_ccg =     sqrt(sum(sum(sum((T_ccg/sqrt(normT)-T).^2))));

ert=[ert,err_ccg];

end

fprintf('Reconstruction error : %g  ',err_ccg)



end

