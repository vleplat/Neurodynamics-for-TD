function [pbest,pbest_val,DI,y] = CNO_PSO(X,options_gen,options_CS,options_DS)
%CNO_PSO Solve systems of ODE with PSO
%% ----------------------------------------------------
% Load parameters
% -----------------------------------------------------
%%% Number of RNNs
NN=options_gen.NN;
sizeX = size(X);

%%% Ranks
R = options_gen.R;

%%% Main parameters
maxIter = options_gen.maxIter;
verbose = options_gen.verbose;
selAlgo = options_gen.selAlgo;
initType = options_gen.initType;
maxTime = options_CS.maxTime;
epsilon.eps_1=options_CS.epsilon(1);
epsilon.eps_2=options_CS.epsilon(2);
epsilon.eps_3=options_CS.epsilon(3);

%% ----------------------------------------------------
% Init RNN
% -----------------------------------------------------
%%% Initial points for each RNN
for i=1:NN
    if initType == 1
        B{i,1}=rand(sizeX(1),R);
        B{i,2}=rand(sizeX(2),R);
        B{i,3}=rand(sizeX(3),R);
    else 
        opts = cp_init();
        opts.init = {'orth' 'fiber' 'rand' };
        Ui = cp_init(tensor(X),RR,opts);
        B{i,1} = Ui{1};
        B{i,2} = Ui{2};
        B{i,3} = Ui{3};
    end
    y0{i,1}=[B{i,1}(:);B{i,2}(:);B{i,3}(:)];
end

%% ----------------------------------------------------
% Solve set of ODE by CNO along with PSO
% -----------------------------------------------------
%%% Init of misc variables
ee = zeros(NN,maxIter);
DI = zeros(1,maxIter);

%%% Parameters for solver based on ode45
tspan = linspace(0,0.01,maxTime); % [0 0.02]

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
                [err_ode,~,~,B] = ALS_ODE(X,epsilon,y0{i,1},tspan,options_CS);
                y{i,j} = [B{1}(:);B{2}(:);B{3}(:)];
            case 'discrete'
                %%% Call of Discrete Solver
                options_DS.U0{1} = reshape(y0{i,1}(1:sizeX(1)*R),[sizeX(1) R]);
                options_DS.U0{2} = reshape(y0{i,1}(sizeX(1)*R+1:sizeX(1)*R+sizeX(2)*R),[sizeX(2) R]);
                options_DS.U0{3} = reshape(y0{i,1}(sizeX(1)*R + sizeX(2)*R + 1:end),[sizeX(3) R]);
                [err_ode,~,A_1,A_2,A_3,~] = Algorithm2(X,options_DS);
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
        [~,jj]=min(ee(i,1:j));
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

end