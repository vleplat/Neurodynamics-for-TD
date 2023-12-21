clc;clear all; close all
rng(2024)
% ------------------------------------------------------------------------
% Loading data set
%--------------------------------------------------------------------------
%warning off

ex='cuprite';
% ex='urban';
% ex='jasper';
% ex='SanDiego';
% ex='Samson';
switch  ex
    case 'cuprite'
        disp('Cuprite data set selected...')
        load('V.mat');
        X = V; clear V;
        % X=X(30:80,30:80,:); % X=X(:,:,1:180);
        X=X(:,:,1:180);
        X=X/max(X(:));
        X=max(X,eps);
        R=12;
     case 'urban'
        disp('Urban data set selected...')
        load('Urban.mat');
        X = reshape(A,[307 307 162]); clear A;
        X=X(120:120+119,120:120+119,:); %X=X(120:120+119,120:120+119,:);
%         X=X/max(X(:));
%         X=max(X,eps);
        R=6;

    case 'jasper'
        disp('Jasper Ridge data set selected...')
        load('jasperRidge2_R198.mat');
        X = reshape(Y',[100 100 198]); clear Y;
        % X=X/max(X(:));
        % X=max(X,eps);
        R=4;

    case 'SanDiego'
        disp('SanDiego data set selected...')
        load('SanDiego.mat');
        X = reshape(A,[400 400 158]); clear A;
        X=X(120:120+119,120:120+119,:);
        X=X/max(X(:));
        X=max(X,eps);
        R=4;

    case 'Samson'
        disp('Samson data set selected...')
        load('samson_1.mat');
        X = reshape(V',[95 95 156]); clear V;
        % X=X/max(X(:));
        %X=max(X,eps);
        R=3;

end
Szx=size(X);
N = ndims(X);
%% ------------------------------------------------------------------------
% Generate common initial points
%--------------------------------------------------------------------------
opts = ncp_hals;
opts.init = 'rand';
opts.maxiters = 10;
opts.tol = 1e-10;
[Yx,out] = ncp_hals(tensor(X),R,opts);
B0 = Yx.U;

%% ------------------------------------------------------------------------
% Call of ALS+ODE45 solver (proposed by Prof Phan)
%--------------------------------------------------------------------------
% Time constant for three-scale neurodynamics
alpha = 0.5;
epsilon.eps_1=alpha*1e-4;
epsilon.eps_2=alpha*1e-4;
epsilon.eps_3=alpha*1e-4;
% 
% epsilon.eps_1=1;
% epsilon.eps_2=1;
% epsilon.eps_3=1;


% Initial point for ODE
theta0 = fac2vec(B0);
tspan = linspace(0,0.004,100); % [0 0.02] % tspan = linspace(0,0.01,50);
% tspan = [0 0.01];
% Misc options 
options = odeset;
options.NonNegative = 1;
options.maxKrun = 1;
options.R = R;
options.algo_Sel = 'als2'; % 'als', 'als2', 'hals2', 'hals'

% Call of solver
[err_ode,err_ode2,cpu_time,B_ode] = ALS_ODE(X,epsilon,theta0,tspan,options);
time_0 = cpu_time;
err_ode(end)

%%
% Pn = normalize(ktensor(B_ode));
% theta0 = fac2vec(Pn.U);
% [err_ode,err_ode2,cpu_time,B_ode] = ALS_ODE(X,epsilon,theta0,tspan,options);
% time_0 = cpu_time;

%% ------------------------------------------------------------------------
% Call of ANLS solver (Put reference)
%--------------------------------------------------------------------------

% A_1=B0{1};
% A_2=B0{2};
% A_3=B0{3};
% tic
% [A_aux,B_aux,C_aux,T_anls,err_anls,ert_anls]= anls(X,550,A_1,A_2,A_3);
% % norm(T_als(:)-X(:))/norm(X(:))
% time_1=toc;


%% ------------------------------------------------------------------------
% Call of proco_als solver (Put reference)
%--------------------------------------------------------------------------
% A_1=B0{1};
% A_2=B0{2};
% A_3=B0{3};
% tic
% [A,B,C,Lambda,T_proco_als,err_proco_als,ert_proco_als] = proco_als(X,550,Szx,A_1,A_2,A_3);
% time_2=toc;

%% ------------------------------------------------------------------------
% Call of ccg solver (Put reference)
%--------------------------------------------------------------------------
% tic
% [ A_ccg, B_ccg, C_ccg, Lambda, T_ccg, err_ccg,ert_ccg] = ccg(X,550,Szx,A_1,A_2,A_3,[.01,.01,.01],0.1,R,'exp');
% time_3=toc;

%% ------------------------------------------------------------------------
% Call of cgP solver (Put reference)
%--------------------------------------------------------------------------
% A_1=B0{1};
% A_2=B0{2};
% A_3=B0{3};
% tic
% [A_cgP B_cgP C_cgP ert_cgP] = cgP(X, R, [],A_1,A_2,A_3);
% time_4=toc;

%% ------------------------------------------------------------------------
% Call of bfgsP solver (Put reference)
%--------------------------------------------------------------------------
% A_1=B0{1};
% A_2=B0{2};
% A_3=B0{3};
% tic
% [A_bfgsP B_bfgsP C_bfgsP error_bfgsP] = bfgsP(X, R, [], A_1,A_2,A_3);
% time_5=toc;

%% ------------------------------------------------------------------------
% Call of gradP solver (Put reference)
% %--------------------------------------------------------------------------
% A_1=B0{1};
% A_2=B0{2};
% A_3=B0{3};
% tic
% [A_gradP B_gradP C_gradP error_gradP] = gradP(X, R, [], A_1,A_2,A_3);
% time_6=toc;

%% ------------------------------------------------------------------------
% Call of HALS solver (Put reference)
%--------------------------------------------------------------------------
opts = ncp_hals;
opts.init = B0;
opts.maxiters = 100;
opts.tol = 1e-10;
tic
[Yx_ncp_hals,out_ncp_hals] = ncp_hals(tensor(X),R,opts);
time_7=toc;
err_hals = norm(X - full(ktensor(Yx_ncp_hals)))/norm(X(:));
yt_ncp_hals=1-out_ncp_hals.Fit(:,2);

%%
%% ----------------------------------------------------
% Setting parameters for Solvers
% -----------------------------------------------------
%%% Main parameters for CNO-PSO
% options_gen.NN = 6;
% options_gen.R = R;
% options_gen.maxIter = 10;
% options_gen.verbose = 1;
% %continuous  discrete
% options_gen.selAlgo = 'continuous' ; 
% % 'discrete'   : Algorithm 2
% % 'continuous' : ALS + ODE45 
% options_gen.initType = 1;
% 
% %%% Parameters for solver based on ode45
% options_CS = odeset;
% options_CS.NonNegative = 1;
% options_CS.maxKrun = 10;
% options_CS.R = R;
% options_CS.algo_Sel = 'als2'; % 'als', 'als2', 'hals2', 'hals'
% options_CS.maxTime = 50;
% options_CS.epsilon(1)=10^(-4);
% options_CS.epsilon(2)=10^(-4);
% options_CS.epsilon(3)=10^(-4);
% 
% %%% Parameters for the Discrete Solver
% options_DS.maxIter = 300;
% options_DS.verbose = 0;
% options_DS.initType = 3;
% options_DS.R = R; 
% options_DS.beta =.2;
% options_DS.alpha =.2;
% options_DS.delta = 0;
% options_DS.AlgoSel = 3;
% options_DS.preCondiSel = 1;
% 
% % Call of Solvers with PSO
% % -----------------------------------------------------
% [pbest,pbest_val,DI,y] = CNO_PSO(X,options_gen,options_CS,options_DS);


%% ------------------------------------------------------------------------
% Call of ncp_mls solver (Put reference)
%--------------------------------------------------------------------------
opts = ncp_mls;
opts.init = B0;
opts.maxiters = 100;
opts.tol = 1e-10;
tic
[Yx_ncp_mls,out_ncp_mls] = ncp_mls(tensor(X),R,opts);
time_8=toc;
err_hals = norm(X - full(ktensor(Yx_ncp_mls)))/norm(X(:));
yt_ncp_mls=1-out_ncp_mls.Fit(:,2);


%% ------------------------------------------------------------------------
% Post-processing
%--------------------------------------------------------------------------
fontSize = 14;
t=1:10:100; %t=1:10:550

figure(1)
clf
% semilogy([err_ode(1:end-10) err_ode2(1:end-10)],'-','LineWidth',3)
semilogy([err_ode(1:end-10) ],'-.','LineWidth',3)
% semilogy([err_ode(1:end-10)],'-','LineWidth',3)
% semilogy([err_ode],'--','LineWidth',3)
hold on;text{1} = "ALS-ODE45";
%semilogy(t,ert_anls(t),'->','LineWidth',3,'color','blue')
%hold on;text{2} = "ANLS";
%semilogy(t,ert_proco_als(t),'-^','LineWidth',3)
%hold on;text{3} = "Proco-ALS";
%semilogy(t,ert_ccg(t),'-<','LineWidth',3,'color','red')
%hold on;text{4} = "CCG";
%semilogy(t,ert_cgP(t),'-o','LineWidth',3,'color','green')
%hold on;text{5} = "CGP";
%semilogy(t,error_bfgsP(t),'-s','LineWidth',3)
%hold on;text{6} = "BFGSP";
%semilogy(t,error_gradP(t),'-diamond','LineWidth',3)
%hold on;text{7} = "GradP";
semilogy(t,yt_ncp_mls(t),'-hexagram','color','black','LineWidth',3)
hold on;text{9} = "MUR";
gg=size(yt_ncp_hals);
% gg =100;
t=1:10:gg;
semilogy(t,yt_ncp_hals(t),'-x','LineWidth',3)
hold on;text{8} = "HALS";
xlabel('iteration counter $k$',"Interpreter","latex")
ylabel('$\| \mathcal{X} - \mathcal{I} \times_1 U^{(1)} \times_2 U^{(2)} \times_3 U^{(3)} \|_F^2$',"Interpreter","latex")
grid on;
legend('ODE','MUR','HALS',"Interpreter","latex")


% legend(text,'interpreter','latex')
% xlabel('Number of iterations','fontsize',fontSize,'interpreter','latex')
% ylabel('$\frac{\|X - \tilde{X}\|_F}{\|X\|_F}$','fontsize',fontSize,'interpreter','latex')
% title('Benchmark results','fontsize',fontSize,'interpreter','latex')
% grid on
% 
% figure(2)
% Xg = categorical({'ODE','ANLS','Proco-ALS','CCG','GradP','CGP','BFGSP','HALS','MUR'});
% Xg = reordercats(Xg,{'ODE','ANLS','Proco-ALS','CCG','GradP','CGP','BFGSP','HALS','MUR'});
% Yg = [time_0 time_1 time_2 time_3 time_4 time_5 time_6 time_7 time_7];
% bar(Xg,Yg)
% ylabel('Running time (Second)')

figure;
for r=1:R
    plot(B_ode{3}(:,r)/max(B_ode{3}(:,r))); hold on;
end
grid on
xlabel("Wavelength Id","Interpreter","latex")
ylabel("Intensity","Interpreter","latex")
%%

% figure(2)
% Pn = normalize(ktensor(B_ode));
% 
% subplot(3,3,1)
% imagesc(Pn.U{1})
% subplot(3,3,2)
% imagesc(Pn.U{2})
% subplot(3,3,3)
% imagesc(Pn.U{3})
% 
% subplot(3,3,4)
% imagesc(Yx_ncp_hals.U{1})
% subplot(3,3,5)
% imagesc(Yx_ncp_hals.U{2})
% subplot(3,3,6)
% imagesc(Yx_ncp_hals.U{3})
% 
% subplot(3,3,7)
% imagesc(Yx_ncp_mls.U{1})
% subplot(3,3,8)
% imagesc(Yx_ncp_mls.U{2})
% subplot(3,3,9)
% imagesc(Yx_ncp_mls.U{3})
