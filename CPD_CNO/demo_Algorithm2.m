clc;clear all
%% ---------------------------------------------------------
% Data init.
% ----------------------------------------------------------
dataSetSel = 4; % Select data set 4
switch dataSetSel
    case 1
        disp('ORL data set selected...')
        load('ORL_32x32.mat');
        AA=permute(reshape(fea,[400,32,32]),[2,3,1]);
        K=400;

    case 2
        disp('Yale data set selected...')
        % load('Yale_32x32.mat');
        load('Yale_64x64.mat');
        si = 64; %64
        AA=permute(reshape(fea,[165,si,si]),[2,3,1]);
        K=165;

    case 3
        disp('COIL20 data set selected...')
        load('COIL20.mat');
        AA=permute(reshape(fea,[1440,32,32]),[2,3,1]);
        K=1440;

    case 4
        disp('Cuprite HSI selected...')
        load('V.mat');
        AA = V;
        K=180;

    case 5
        disp('SanDiedo HSI selected...')
        [Xsub] = SanDiego_preProc(1);
        AA=permute(reshape(Xsub,[5,400,400]),[2,3,1]);
        K=5;

    otherwise
        disp('wrong selection for the data set...')
        return
end

%%% selection of a subset of data
X=double(AA(:,:,1:K));
X=X/max(X(:));
tX = tensor(X);

%%% Parameters for the Discrete Solver
options.maxIter = 500;
options.verbose = 1;
options.initType = 3;
options.R =12; %11 for Yale, 12 for Cuprite
options.beta =.2;
options.alpha =.2;
options.delta = 0;

%%% Init. for the factors
% Random
sizeX = size(X);
A_10 = rand(sizeX(1),options.R);
A_20 = rand(sizeX(2),options.R);
A_30 = rand(K,options.R);
options.U0{1} = A_10;
options.U0{2} = A_20;
options.U0{3} = A_30;

%% ---------------------------------------------------------
% Call of Algorithm 2 variants
% ----------------------------------------------------------
% Algorithm 2 - Hessian preconditionning, semi-implicit
options.AlgoSel = 3;
options.preCondiSel = 1;
[err,cpu_time,A_1,A_2,A_3,lamTab] = Algorithm2(X,options);

% Algorithm 2 - Identity preconditionning, fully explicit
options.AlgoSel = 1;
options.preCondiSel = 2;
[err_t2,cpu_time_t2,A_1_t2,A_2_t2,A_3_t2,lamTab_t2] = Algorithm2(X,options);

% Algorithm 2 - Cubic-regularized method
options.AlgoSel = 2;
[err_t3,cpu_time_t3,A_1_t3,A_2_t3,A_3_t3,lamTab_t3] = Algorithm2(X,options);

% HALS
opts = ncp_hals;
opts.init = {'rand' 'rand' 'rand' };
opts.maxiters = options.maxIter;
[P,output] = ncp_hals(tX,options.R,opts);

%% ---------------------------------------------------------
% Post-processing
% ----------------------------------------------------------
R = options.R;
iter = options.maxIter;
t=1:iter;
close all

%%% Convergence plots
fontSize = 14;
figure(1);clf; semilogy(output.Fit(:,1),1-output.Fit(:,2),'-<','LineWidth',0.5,'color','red');hold on;text{1} = "HALS";
semilogy(err(1:end),'->','LineWidth',0.5,'color','blue'); hold on;text{2} = "Algo2 - SemiImplicit - $P=H$";
semilogy(err_t2(1:end),'-o','LineWidth',0.5,'color','green'); hold on;text{3} = "Algo2 - Full Explicit - $P=I$";
semilogy(err_t3(1:end),'-diamond','LineWidth',0.5); hold on;text{4} = "Algo2 - Cubic Reg.";

legend(text,'interpreter','latex')
xlabel('Number of iterations','fontsize',fontSize,'interpreter','latex')
ylabel('$\frac{\|X - \tilde{X}\|_F}{\|X\|_F}$','fontsize',fontSize,'interpreter','latex')
title('Benchmark results','fontsize',fontSize,'interpreter','latex')
grid on



%%% Transient responses - Semi-implicit

for i=1:iter
    Ahh1(i,:)=A_1{i}(:)';
    Ahh2(i,:)=A_2{i}(:)';
    Ahh3(i,:)=A_3{i}(:)';
end

figure;
subplot(2,2,1);
plot(t,Ahh1(:,1:end));
title('TR - Factor 1 - SE')
ylabel('A^{(1)}')
xlabel('Time (t)')

subplot(2,2,2);
plot(t,Ahh2(:,1:end));
title('TR - Factor 2 - SE')
ylabel('A^{(2)}')
xlabel('Time (t)')

subplot(2,2,[3,4]);
plot(t,Ahh3(:,1:end));
title('TR - Factor 3 - SE')
xlabel('Time (t)')
ylabel('A^{(3)}')
%figure(2)
% plot(t,Ahh1(:,1:end),'LineWidth',1.5)
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(1)}')

%figure(3)
% plot(t,Ahh2(:,1:end),'LineWidth',1.5)
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(2)}')

% figure(4)
% plot(t,Ahh3(:,1:end),'LineWidth',1.5)
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(3)}')



%%% Transient responses - Fully explicit

for i=1:iter
    Ahh1(i,:)=A_1_t2{i}(:)';
    Ahh2(i,:)=A_2_t2{i}(:)';
    Ahh3(i,:)=A_3_t2{i}(:)';
end
figure;
subplot(2,2,1);
plot(t,Ahh1(:,1:end));
title('TR - Factor 1 - FE')
ylabel('A^{(1)}')
xlabel('Time (t)')

subplot(2,2,2);
plot(t,Ahh2(:,1:end));
title('TR - Factor 2 - FE')
ylabel('A^{(2)}')
xlabel('Time (t)')

subplot(2,2,[3,4]);
plot(t,Ahh3(:,1:end));
title('TR - Factor 3 - FE')
xlabel('Time (t)')
ylabel('A^{(3)}')
% figure(5)
% plot(t,Ahh1(:,1:end),'LineWidth',1.5)
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(1)}')
% 
% figure(6)
% plot(t,Ahh2(:,1:end),'LineWidth',1.5)
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(2)}')
% 
% figure(7)
% plot(t,Ahh3(:,1:end),'LineWidth',1.5)
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(3)}')




%%% cuprite - Display estimated spectral signatures
if dataSetSel == 4
    figure(5)
    for r=1:R
        subplot(3,4,r); 
        plot(A_3{i+1}(:,r))
    end
end


