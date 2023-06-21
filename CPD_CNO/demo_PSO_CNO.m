clear all;clc;close all;
warning("off");
%% ----------------------------------------------------
% Data loading
% -----------------------------------------------------
dataSetSel = 7; 
sizeSynthetic(1) = 9;
sizeSynthetic(2) = 13;
[X,R] = data_Loader(dataSetSel,sizeSynthetic);
% 1: ORL data set
% 2: Yale data set
% 3: COIL20 data set
% 4: Cuprite HSI
% 5: SanDiedo HSI
% 6: mnist data set
% 7: synthetic data set
X=X/max(X(:));
tX = tensor(X);


%% ----------------------------------------------------
% Setting parameters for Solvers
% -----------------------------------------------------
%%% Main parameters for CNO-PSO
options_gen.NN = 6;
options_gen.R = R;
options_gen.maxIter = 10;
options_gen.verbose = 1;
options_gen.selAlgo = 'continuous' ; 
% 'discrete'   : Algorithm 2
% 'continuous' : ALS + ODE45 
options_gen.initType = 1;

%%% Parameters for solver based on ode45
options_CS = odeset;
options_CS.NonNegative = 1;
options_CS.maxKrun = 10;
options_CS.R = R;
options_CS.algo_Sel = 'als2'; % 'als', 'als2', 'hals2', 'hals'
options_CS.maxTime = 50;
options_CS.epsilon(1)=10^(-4);
options_CS.epsilon(2)=10^(-4);
options_CS.epsilon(3)=10^(-4);

%%% Parameters for the Discrete Solver
options_DS.maxIter = 300;
options_DS.verbose = 0;
options_DS.initType = 3;
options_DS.R = R; 
options_DS.beta =.2;
options_DS.alpha =.2;
options_DS.delta = 0;
options_DS.AlgoSel = 3;
options_DS.preCondiSel = 1;

%% ----------------------------------------------------
% Call of Solvers with PSO
% -----------------------------------------------------
[pbest,pbest_val,DI,y] = CNO_PSO(X,options_gen,options_CS,options_DS);
