
%------------------------------Run_Decomp.m-------------------------------%
% In this file, different versions of the compressed algorithms are tested.
% The compressed versions can be compared to the uncompressed projected ALS 
% algorithm. 
% 
% A randomized approximate version of the svd is computed for numerical co-
% -mplexity reasons.
%
% The squared L2 norm of the measurement reconstruction error is evaluated.
%-------------------------------------------------------------------------%

%----------------------------Cleaning variables---------------------------%
clc;
clear all;
close all;
%-------------------------------------------------------------------------%

%-----------------------------Model parameters----------------------------%

% --- YOUR DATA --- %

%load('snow_subsample3'); % the nonnegative data tensor has to be called Y
%Y=my_data
%Y=miss_data(Y) % To handle missing data by linear interpolation along the
%third mode.

% --- SIMULATED DATA --- %

% --- PARAMETERS --- %

% CP number of components
R         =     5;
% Uncompressed dimensions
dim_u     =     [4800 200 70];
% Compressed dimensions
dim_c     =     [30 10 5];
% Number of iterates for the ALS algorithms
itermax_ccg     =     500;
itermax_als     =     500;

sigma=10^(-4);
% CP model and noisy measurements (additive Gaussian noise)
[A,B,C,~,Y]     = tensor_gen(dim_u,R,sigma);



%-------------------------------------------------------------------------%

%-----------------------------Initialization------------------------------%
Re =     R; % estimated rank

A0 =     abs(randn(dim_u(1),Re));
B0 =     abs(randn(dim_u(2),Re));
C0 =     abs(randn(dim_u(3),Re));
C0 =     C0.*repmat(1./sqrt(sum(C0.^2)),dim_u(3),1);
A0 =     A0.*repmat(1./sqrt(sum(A0.^2)),dim_u(1),1);
B0 =     B0.*repmat(1./sqrt(sum(B0.^2)),dim_u(2),1);
%-------------------------------------------------------------------------%

%--------------------------------  CCG  ----------------------------------%
% Non-linear Compressed Penalized conjugate gradient

% Start time counter
t_ccg_start     =     tic;

% CCG
alpha     =     [5,5,5];

gamma     =     0.01;

[Accg, Bccg, Cccg, ~, Tccg, err_ccg] =  ccg(Y,itermax_ccg,dim_c,A0,B0,C0,alpha,gamma,Re,'exp');


% Stop time counter
t_ccg     =     toc(t_ccg_start);

%---------------------------------ProCo ALS----------------------------------%
% Standard CP-ALS

% Start time counter
t_cp_als_start              =     tic;

% Compressed projected ALS
[Aals, Bals, Cals, ~, Tals, err_cpals]    =     proco_als(Y,itermax_als,dim_c,A0,B0,C0);

% Stop time counter
t_cp_als  =     toc(t_cp_als_start);

%-------------------------------------------------------------------------%


% %----------------------------------ANLS---------------------------------%
% % (OPTIONAL)
% 
% % Start time counter
% t_anls_start                =     tic;
% 
% % Projected ALS
% [Aanls,Banls,Canls,Tanls,err_anls]=     anls(Y,itermax_als,A0,B0,C0);
% 
% % Stop time counter
% t_anls          =     toc(t_anls_start);

%-------------------------------------------------------------------------%


