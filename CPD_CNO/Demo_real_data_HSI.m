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
        X=X(:,:,1:180);
        X=X/max(X(:));
        X=max(X,eps);
        R=12;
     case 'urban'
        disp('Urban data set selected...')
        load('Urban.mat');
        X = reshape(A,[307 307 162]); clear A;
        X=X(120:120+119,120:120+119,:);
        R=6;

    case 'jasper'
        disp('Jasper Ridge data set selected...')
        load('jasperRidge2_R198.mat');
        X = reshape(Y',[100 100 198]); clear Y;
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
% Call of ALS+ODE45 solver
%--------------------------------------------------------------------------
% Time constant for three-scale neurodynamics
alpha = 0.5;
epsilon.eps_1=alpha*1e-4;
epsilon.eps_2=alpha*1e-4;
epsilon.eps_3=alpha*1e-4;

% Initial point for ODE
theta0 = fac2vec(B0);
tspan = linspace(0,0.004,100);
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
fontSize = 20;
t=1:10:100;

% % Error plots
figure(1)
set(0, 'DefaultAxesFontSize', fontSize);
clf
semilogy([err_ode(1:end-10) ],'-.','LineWidth',3)
hold on;text{1} = "ODE";
semilogy(t,yt_ncp_mls(t),'-hexagram','color','black','LineWidth',3)
hold on;text{2} = "MUR";
gg=size(yt_ncp_hals);
t=1:10:gg;
semilogy(t,yt_ncp_hals(t),'-x','LineWidth',3)
hold on;text{3} = "HALS";
% xlabel('iteration counter $k$',"Interpreter","latex",'FontSize',fontSize)
%ylabel('$\frac{\| \mathcal{X} - \mathcal{I} \times_1 U^{(1)} \times_2 U^{(2)} \times_3 U^{(3)} \|_F^2}{\| \mathcal{X} \|_F^2}$',"Interpreter","latex",'FontSize',fontSize)
xlabel('Number of iterations',"Interpreter","latex",'FontSize',fontSize)
ylabel('Relative Error',"Interpreter","latex",'FontSize',fontSize)
grid on;
legend(text,'FontSize',fontSize, 'Interpreter','latex')


% % Spectral signatures
figure;
for r=1:R
    plot(B_ode{3}(:,r)/max(B_ode{3}(:,r))); hold on;
end
grid on
xlabel("Wavelength Id","Interpreter","latex")
ylabel("Intensity","Interpreter","latex")

