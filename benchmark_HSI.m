clear all;clc;close all;

%% ------------------------------------------------------------------------
% Benchmark for HSI tests
%%-------------------------------------------------------------------------
% Tip: in the script demo_real_data_HSI, select the disired data set
nb_trials = 10;
ODE_err = [];
MUR_err = [];
HALS_err = [];
seed = [];

for it=1:nb_trials

    % setup a new seed for pseudo-random generator
    seed = [seed rng(it)];

    % call of demo_real_data_HSI
    run demo_real_data_HSI.m

    % load variables
    load variables.mat;

    % Save results
    ODE_err =  [ODE_err err_ode(1:end-10)'];
    MUR_err = [MUR_err yt_ncp_mls];
    HALS_err = [HALS_err yt_ncp_hals];

    close all;

end

%% ------------------------------------------------------------------------
% Post-processing
%%-------------------------------------------------------------------------
% Computation of the median
median_MUR = median(MUR_err,2);
median_HALS = median(HALS_err,2);
median_ODE = median(ODE_err,2);

% Graph generation
fontSize = 20;
t=1:10:100;

% % Error plots
figure(1)
set(0, 'DefaultAxesFontSize', fontSize);
clf
semilogy(median_ODE,'-.','LineWidth',3)
hold on;text{1} = "ODE";
semilogy(t,median_MUR(t),'-hexagram','color','black','LineWidth',3)
hold on;text{2} = "MUR";
gg=size(median_HALS);
t=1:10:gg;
semilogy(t,median_HALS(t),'-x','LineWidth',3)
hold on;text{3} = "HALS";
% xlabel('iteration counter $k$',"Interpreter","latex",'FontSize',fontSize)
%ylabel('$\frac{\| \mathcal{X} - \mathcal{I} \times_1 U^{(1)} \times_2 U^{(2)} \times_3 U^{(3)} \|_F^2}{\| \mathcal{X} \|_F^2}$',"Interpreter","latex",'FontSize',fontSize)
xlabel('Number of iterations',"Interpreter","latex",'FontSize',fontSize)
ylabel('Relative Error',"Interpreter","latex",'FontSize',fontSize)
grid on;
title(['Median curves over ' num2str(nb_trials) ' runs'],'FontSize',fontSize, 'Interpreter','latex')
legend(text,'FontSize',fontSize, 'Interpreter','latex')