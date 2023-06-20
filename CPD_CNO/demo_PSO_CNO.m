clear all;clc;close all;
warning("off");
%% ----------------------------------------------------
% Data loading
% -----------------------------------------------------
dataSetSel = 7; 
switch dataSetSel
    case 1
        disp('ORL data set selected...')
        load('ORL_32x32.mat');
        AA=permute(reshape(fea,[400,32,32]),[2,3,1]);
        K=400;
        %%% selection of a subset of data
        X = double(AA(:,:,1:K));
        R=10;

    case 2
        disp('Yale data set selected...')
        % load('Yale_32x32.mat');
        load('Yale_64x64.mat');
        si = 64; %64
        AA=permute(reshape(fea,[165,si,si]),[2,3,1]);
        K=165;
        %%% selection of a subset of data
        X = double(AA(:,:,1:K));
        R=11;

    case 3
        disp('COIL20 data set selected...')
        load('COIL20.mat');
        AA=permute(reshape(fea,[1440,32,32]),[2,3,1]);
        K=1440;
        %%% selection of a subset of data
        X = double(AA(:,:,1:K));
        R=11; 

    case 4
        disp('Cuprite HSI selected...')
        load('V.mat');
        AA = V;
        clear V;
        K=180;
        %%% selection of a subset of data
        X = double(AA(:,:,1:K));
        R=12; 

    case 5
        disp('SanDiedo HSI selected...')
        [Xsub] = SanDiego_preProc(1);
        AA=permute(reshape(Xsub,[5,400,400]),[2,3,1]);
        K=5;
        %%% selection of a subset of data
        X = double(AA(:,:,1:K));
        R=8; 
    case 6
        disp('mnist data set selected...')
        load('mnist_all.mat');
        AA_cell{1}=test0; AA_cell{2}=test1; AA_cell{3}=test2;
        AA_cell{4}=test3; AA_cell{5}=test4; AA_cell{6}=test5;
        AA_cell{7}=test6; AA_cell{8}=test7; AA_cell{9}=test8;
        AA_cell{10}=test9; 
        for j = 1:10
            for i=1:140
                AA(:,:,(j-1)*140+i) = (reshape(AA_cell{j}(i,:),[28,28]))';
            end
        end
        K=1400;
        %%% selection of a subset of data
        X = double(AA(:,:,1:K));
        R=10;
    
    case 7
        disp('synthetic data set selected...')
        n=30;  %30
        R=45; %45     
        % Generating a nonnegative tensor
        A{1}=rand(n,R);
        A{2}=rand(n,R);
        A{3}=rand(n,R);
        X=double(full(ktensor(ones(R,1),A{1},A{2},A{3})));

    otherwise
        disp('wrong selection for the data set...')
        return
end

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
options_gen.selAlgo = 'discrete' ; 
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
