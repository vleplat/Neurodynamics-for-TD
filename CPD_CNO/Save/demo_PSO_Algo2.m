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
% Init RNN
% -----------------------------------------------------
%%% Number of RNNs
NN=6;
sizeX = size(X);
%%% Initial points for each RNN
initType = 1;
for i=1:NN
    if initType == 1
        B{i,1}=rand(sizeX(1),R);
        B{i,2}=rand(sizeX(2),R);
        B{i,3}=rand(sizeX(3),R);
    else 
        opts = cp_init();
        opts.init = {'orth' 'fiber' 'rand' };
        Ui = cp_init(tensor(X),R,opts);
        B{i,1} = Ui{1};
        B{i,2} = Ui{2};
        B{i,3} = Ui{3};
    end
    y0{i,1}=[B{i,1}(:);B{i,2}(:);B{i,3}(:)];
end


%% ----------------------------------------------------
% Solve set of ODE by CNO along with PSO
% -----------------------------------------------------
%%% Main parameters
maxIter = 10;
verbose = 1;
selAlgo = 'discrete' ; % 'discrete','continuous'

%%% Init of misc variables
ee = zeros(NN,maxIter);
DI = zeros(1,maxIter);

%%% Parameters for solver based on ode45
epsilon.eps_1=10^(-4);
epsilon.eps_2=10^(-4);
epsilon.eps_3=10^(-4);
maxTime = 50;
tspan = linspace(0,0.01,maxTime); % [0 0.02]
options_CS = odeset;
options_CS.NonNegative = 1;
options_CS.maxKrun = 10;
options_CS.R = R;
options_CS.algo_Sel = 'als2'; % 'als', 'als2', 'hals2', 'hals'

%%% Parameters for the Discrete Solver
options_DS.maxIter = 300;
options_DS.verbose = 0;
options_DS.initType = 3;
options_DS.R = R; %11 for Yale, 12 for Cuprite
options_DS.beta =.2;
options_DS.alpha =.2;
options_DS.delta = 0;
options_DS.AlgoSel = 3;
options_DS.preCondiSel = 1;

%%% Init for velocity vectors for PSO method
fprintf('--------------------------------------------------------------------------------------')
fprintf('\n')
fprintf('------------------------------------ PSO Progress ------------------------------------')
fprintf('\n')
for i=1:NN
    V{i}= rand(sum(sizeX*R),1);
    text{i}=sprintf('RNN-%d           ', i);
    fprintf(text{i})
end
fprintf('\n')
fprintf('--------------------------------------------------------------------------------------')
fprintf('\n')

%%% Main Loop
for j=1:maxIter

    %%% Loop on each RNN
    for i=1:NN
        switch selAlgo
            case 'continuous'
                %%% Call of ODE45-based solver
                [err_ode,err_ode2,cpu_time,B] = ALS_ODE(X,epsilon,y0{i,1},tspan,options_CS);
                y{i,j} = [B{1}(:);B{2}(:);B{3}(:)];
            case 'discrete'
                %%% Call of Discrete Solver
                options_DS.U0{1} = reshape(y0{i,1}(1:sizeX(1)*R),[sizeX(1) R]);
                options_DS.U0{2} = reshape(y0{i,1}(sizeX(1)*R+1:sizeX(1)*R+sizeX(2)*R),[sizeX(2) R]);
                options_DS.U0{3} = reshape(y0{i,1}(sizeX(1)*R + sizeX(2)*R + 1:end),[sizeX(3) R]);
                [err_ode,cpu_time,A_1,A_2,A_3,lamTab] = Algorithm2(X,options_DS);
                y{i,j} = [A_1{end}(:);A_2{end}(:);A_3{end}(:)];
            otherwise
                disp('wrong selection for solver...')
            return
        end
    
        ee(i,j) = err_ode(end);
        
    end

    %%% Print errors for each particle
    if verbose == 1
        for i=1:NN
            text{i}=sprintf('%f        ', ee(i,j));
            fprintf(text{i});
        end
        fprintf('Iteration %d',j)
        fprintf('\n')
    end

    %%% PSO
    %%%% Computation of p_best
    if j==1
        tt_h=min(min(ee(:,1:j)));
        [tt_1,tt_2]=find(ee==tt_h);
        pbest =y{tt_1(1),tt_2(1)};
        pbest_val = tt_h;
    else
        tt_h=min(min(ee(:,j)));
        if tt_h < pbest_val
            [tt_1,tt_2]=find(ee==tt_h);
            pbest =y{tt_1(1),tt_2(1)};
            pbest_val = tt_h;

        end
    end
    
    %%%% Update of the velocity and particle positions
    w=0;
    for i=1:NN
        %%%% Computation of p_n : the best location of the n-th particle
        [ii,jj]=min(ee(i,1:j));
        p_n = y{i,jj};
        [y0{i,1},V{i}]=PSO_Code(V{i},y{i,j},p_n,pbest);
        w=w+norm(p_n-pbest);
    end
    DI(j)=w/NN;
    %%% Print the current DI
    fprintf('D(%d): %f',j,DI(j))
    fprintf('\n')
    fprintf('--------------------------------------------------------------------------------------------------------------')
    fprintf('\n')

end

