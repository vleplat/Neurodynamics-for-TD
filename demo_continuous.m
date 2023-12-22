clc;clear all; close all
%% ------------------------------------------------------------------------
% Loading data set
%--------------------------------------------------------------------------
%warning off
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
% Call of ALS+ODE45 solver
%--------------------------------------------------------------------------
% Time constant for three-scale neurodynamics
epsilon.eps_1=1e-4;
epsilon.eps_2=1e-4;
epsilon.eps_3=1e-4;

% Initial point for ODE
theta0 = fac2vec(B0);
tspan = linspace(0,0.01,50); % [0 0.02]

% Misc options 
options = odeset;
options.NonNegative = 1;
options.maxKrun = 10;
options.R = R;
options.algo_Sel = 'als2'; % 'als', 'als2', 'hals2', 'hals'

% Call of solver
[err_ode,err_ode2,cpu_time,B_ode] = ALS_ODE(X,epsilon,theta0,tspan,options);
time_0 = cpu_time;

%% ------------------------------------------------------------------------
% Call of ANLS solver (Put reference)
%--------------------------------------------------------------------------

A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A_aux,B_aux,C_aux,T_anls,err_anls,ert_anls]= anls(X,550,A_1,A_2,A_3);
% norm(T_als(:)-X(:))/norm(X(:))
time_1=toc;


%% ------------------------------------------------------------------------
% Call of proco_als solver (Put reference)
%--------------------------------------------------------------------------
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A,B,C,Lambda,T_proco_als,err_proco_als,ert_proco_als] = proco_als(X,550,Szx,A_1,A_2,A_3);
time_2=toc;

%% ------------------------------------------------------------------------
% Call of ccg solver (Put reference)
%--------------------------------------------------------------------------
tic
[ A_ccg, B_ccg, C_ccg, Lambda, T_ccg, err_ccg,ert_ccg] = ccg(X,550,Szx,A_1,A_2,A_3,[.01,.01,.01],0.1,R,'exp');
time_3=toc;

%% ------------------------------------------------------------------------
% Call of cgP solver (Put reference)
%--------------------------------------------------------------------------
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A_cgP B_cgP C_cgP ert_cgP] = cgP(X, R, [],A_1,A_2,A_3);
time_4=toc;

%% ------------------------------------------------------------------------
% Call of bfgsP solver (Put reference)
%--------------------------------------------------------------------------
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A_bfgsP B_bfgsP C_bfgsP error_bfgsP] = bfgsP(X, R, [], A_1,A_2,A_3);
time_5=toc;

%% ------------------------------------------------------------------------
% Call of gradP solver (Put reference)
%--------------------------------------------------------------------------
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A_gradP B_gradP C_gradP error_gradP] = gradP(X, R, [], A_1,A_2,A_3);
time_6=toc;

%% ------------------------------------------------------------------------
% Call of HALS solver (Put reference)
%--------------------------------------------------------------------------
opts = ncp_hals;
opts.init = B0;
opts.maxiters = 550;
opts.tol = 1e-10;
tic
[Yx_ncp_hals,out_ncp_hals] = ncp_hals(tensor(X),R,opts);
time_7=toc;
err_hals = norm(X - full(ktensor(Yx_ncp_hals)))/norm(X(:));
yt_ncp_hals=1-out_ncp_hals.Fit(:,2);

%% ------------------------------------------------------------------------
% Call of ncp_mls solver (Put reference)
%--------------------------------------------------------------------------
opts = ncp_mls;
opts.init = B0;
opts.maxiters = 550;
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
t=1:10:550;

figure(1)
clf
semilogy([err_ode(1:end-10) err_ode2(1:end-10)],'--','LineWidth',3)
hold on;text{1} = "ALS-ODE45";
semilogy(t,ert_anls(t),'->','LineWidth',3,'color','blue')
hold on;text{2} = "ANLS";
semilogy(t,ert_proco_als(t),'-^','LineWidth',3)
hold on;text{3} = "Proco-ALS";
semilogy(t,ert_ccg(t),'-<','LineWidth',3,'color','red')
hold on;text{4} = "CCG";
semilogy(t,ert_cgP(t),'-o','LineWidth',3,'color','green')
hold on;text{5} = "CGP";
semilogy(t,error_bfgsP(t),'-s','LineWidth',3)
hold on;text{6} = "BFGSP";
semilogy(t,error_gradP(t),'-diamond','LineWidth',3)
hold on;text{7} = "GradP";
semilogy(t,yt_ncp_hals(t),'-x','LineWidth',3)
hold on;text{8} = "HALS";
semilogy(t,yt_ncp_mls(t),'-hexagram','color','black','LineWidth',3)
hold on;text{9} = "MUR";

legend(text,'interpreter','latex')
xlabel('Number of iterations','fontsize',fontSize,'interpreter','latex')
ylabel('$\frac{\|X - \tilde{X}\|_F}{\|X\|_F}$','fontsize',fontSize,'interpreter','latex')
title('Benchmark results','fontsize',fontSize,'interpreter','latex')
grid on

figure(2)
Xg = categorical({'ODE','ANLS','Proco-ALS','CCG','GradP','CGP','BFGSP','HALS','MUR'});
Xg = reordercats(Xg,{'ODE','ANLS','Proco-ALS','CCG','GradP','CGP','BFGSP','HALS','MUR'});
Yg = [time_0 time_1 time_2 time_3 time_4 time_5 time_6 time_7 time_7];
bar(Xg,Yg)
ylabel('Running time (Second)')

