clear all;clc;close all;
warning("off");
rng(2024)
%% ----------------------------------------------------
% Data loading
% -----------------------------------------------------
dataSetSel = 8; 
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
% 8: Urban data set
X=X/max(X(:));
tX = tensor(X);


%% ----------------------------------------------------
% Setting parameters for Solvers
% -----------------------------------------------------
%%% Main parameters for CNO-PSO
options_gen.NN = 6;
options_gen.R = R;
options_gen.maxIter = 20;
options_gen.verbose = 1;
options_gen.selAlgo = 'continuous' ; 
% 'discrete'   : Algorithm 2
% 'continuous' : ALS + ODE45 
options_gen.initType = 1;

%%% Parameters for solver based on ode45
options_CS = odeset;
options_CS.NonNegative = 1;
options_CS.maxKrun = 1;
options_CS.R = R;
options_CS.algo_Sel = 'als2'; % 'als', 'als2', 'hals2', 'hals'
options_CS.tSpanEnd = 0.004;
options_CS.NbComputationPoints = 100;
alpha = 0.5;
options_CS.epsilon(1)=alpha*10^(-4);
options_CS.epsilon(2)=alpha*10^(-4);
options_CS.epsilon(3)=alpha*10^(-4);

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

%% ----------------------------------------------------
% Post-processing
% -----------------------------------------------------
sizeX = size(X);
idx_s = 0;
for i=1:length(sizeX)
    B{i} = pbest(idx_s+1:idx_s+sizeX(i)*R,:);
    B{i} = reshape(B{i},[sizeX(i) R]);
    idx_s = sizeX(i)*R;
end

%%% Normalize 
sizeB3 = size(B{3});
for i=1:sizeB3(2)
    B{3}(:,i) = B{3}(:,i)/max(B{3}(:,i));
end
%%% Plot spectral signatures of endmembers
figure;plot(B{3})
grid on
xlabel("Wavelength Id","Interpreter","latex")
ylabel("Intensity","Interpreter","latex")

