clc;clear all; close all
%warning off
%ex = 'ORL';
%ex = 'Yale';
ex = 'rand';
switch  ex
    case 'rand'
        n=9;  
        R=10;      
        % Generating a nonnegative tensor
        A{1}=rand(n,R);
        A{2}=rand(n,R);
        A{3}=rand(n,R);
      
        X=double(full(ktensor(ones(R,1),A{1},A{2},A{3})));
        % X=nasob23(3,3,3);R = 23;

    case  'ORL'
        % for mnist
        disp('ORL data set selected...')
        load('ORL_32x32.mat');
        AA=permute(reshape(fea,[400,32,32]),[2,3,1]);
        K=400;
        X=double(AA(:,:,1:K));
        %X=X/max(X(:));
        R=35;

     case 'Yale'
        disp('Yale data set selected...')
        % load('Yale_32x32.mat');
        load('Yale_64x64.mat');
        si = 64; %64
        AA=permute(reshape(fea,[165,si,si]),[2,3,1]);
        K=30;
        X=double(AA(:,:,1:K));
        %X=X/max(X(:));
        R=20;
end
Szx=size(X);
N = ndims(X);
%% Generate initial point
opts = ncp_hals;
opts.init = 'rand';
opts.maxiters = 10;
opts.tol = 1e-10;
[Yx,out] = ncp_hals(tensor(X),R,opts);
B0 = Yx.U;
%% note that h in ODE45 is 6e-7
epsilon.eps_1=1e-4;
epsilon.eps_2=1e-4;
epsilon.eps_3=1e-4;

theta = fac2vec(B0);
tspan = linspace(0,0.01,50); % [0 0.02]
err_ode = [];
normX = norm(X(:));
options = odeset;
options.NonNegative = 1;

%%
tic
for krun = 1:1
    [t_ode,yode] = ode45(@(t,xx) ODE(t,xx,X,epsilon,R), tspan, theta,options);

    clear err_ode2;
    for k = 1:size(yode,1)
        %Bnew = vec2fac(yode(k,:),[size(X) 1]);
        %err_ode(k) = norm(X - full(ktensor(Bnew{4}(:),Bnew(1:3))))/norm(X(:));
        Bnew = vec2fac(yode(k,:),size(X));
        Pb = ktensor(Bnew(1:3));
        err_ode2(k) = sqrt(normX^2+ norm(Pb)^2-2*innerprod(tensor(X),Pb))/normX;
    end
    theta = yode(end,:)';
    
    if ~isempty(err_ode)
        if abs(err_ode2(end)-err_ode(end))<=1e-4*err_ode(end)
            break
        end
    end

    err_ode = [err_ode err_ode2];
     
%     %% Alg1 : EPC
%     ssc = 0;
%     if ssc
%         Y = tensor(X);
%         Up = vec2fac(yode(end,:),size(X));
%         opts = cp_anc;
%         opts.printitn = 1;
%         opts.maxiters = 30;
%         opts.tol = 1e-5;% 1e-5
%         delta0 = norm(tensor(Y) - full(ktensor(Up))) ;
%         expandfactor = min(1.1,0.9*norm(Y)/delta0);
%         delta = delta0*expandfactor;
%         Yx = Up;
%         ss_anc = [1 cp_sensitivity(Yx)];
%         wss_anc = [1 cp_wsensitivity(Up)];
%         err_anc = [1 delta0];
%         for kc = 1:30
%             opts.init = Yx;
%             [Yx,output] = cp_anc(tensor(Y),R,delta,opts);
%             ss_anc = [ss_anc; ss_anc(end,1)+floor(numel(output.cost)/N)+1 cp_sensitivity(Yx)];
%             U = Yx.U;
%             U = cellfun(@(x) x*diag(Yx.lambda.^(1/N)),U,'uni',0);
%             wss_anc = [wss_anc; wss_anc(end,1)+floor(numel(output.cost)/N)+1 cp_wsensitivity(U)];
%             err_anc = [err_anc; err_anc(end,1)+floor(numel(output.cost)/N)+1 norm(Y-full(Yx))];
%         end
%         Uepc = U;
%         theta = fac2vec(U);
% 
%     end
end
time_0=toc
%% HALS
% opts = ncp_hals;
% opts.init = B0;
% opts.maxiters = numel(err_ode);
% opts.tol = 1e-10;
% [Yx,out] = ncp_hals(tensor(X),R,opts);
% err_hals = norm(X - full(ktensor(Yx)))/norm(X(:));
% 
% [err_ode(end) err_hals]

%%
figure(1)
clf
semilogy([err_ode(1:end-10) err_ode2(1:end-10)],'--','LineWidth',3)
hold on
%loglog(1-out.Fit(:,2))
%legend({'ODE' 'HALS'})
%xlabel('Iterations')
%ylabel('Relative Error')
grid on

%%
t=1:10:550
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A_aux,B_aux,C_aux,T_als,err_als,ert]= anls(X,550,A_1,A_2,A_3);
norm(T_als(:)-X(:))/norm(X(:))
time_1=toc
semilogy(t,ert(t),'->','LineWidth',3,'color','blue')
hold on
%%
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A,B,C,Lambda,T_als,err_als,ert] = proco_als(X,550,Szx,A_1,A_2,A_3);
time_2=toc
semilogy(t,ert(t),'-<','LineWidth',3,'color','red')
hold on
%%
tic
[ A_ccg, B_ccg, C_ccg, Lambda, T_ccg, err_ccg,ert] = ccg(X,550,Szx,A_1,A_2,A_3,[.01,.01,.01],0.1,R,'exp');
time_3=toc
semilogy(t,ert(t),'-^','LineWidth',3)
%%
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A B C ert] = cgP(X, R, [],A_1,A_2,A_3);
time_4=toc
semilogy(t,ert(t),'-o','LineWidth',3,'color','green')
hold on
%%
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A B C error] = bfgsP(X, R, [], A_1,A_2,A_3);
time_5=toc
semilogy(t,error(t),'-s','LineWidth',3)
hold on
%%
A_1=B0{1};
A_2=B0{2};
A_3=B0{3};
tic
[A B C error] = gradP(X, R, [], A_1,A_2,A_3);
time_6=toc
semilogy(t,error(t),'-diamond','LineWidth',3)
hold on
%%
opts = ncp_hals;
opts.init = B0;
opts.maxiters = 550;
opts.tol = 1e-10;
tic
[Yx,out] = ncp_hals(tensor(X),R,opts);
time_7=toc
err_hals = norm(X - full(ktensor(Yx)))/norm(X(:));
yt=1-out.Fit(:,2);
semilogy(t,yt(t),'-x','LineWidth',3)
%%
opts = ncp_mls;
opts.init = B0;
opts.maxiters = 550;
opts.tol = 1e-10;
tic
[Yx,out] = ncp_mls(tensor(X),R,opts);
time_8=toc
err_hals = norm(X - full(ktensor(Yx)))/norm(X(:));
yt=1-out.Fit(:,2);
semilogy(t,yt(t),'-hexagram','color','black','LineWidth',3)


legend('ODE','ANLS','Proco-ALS','CCG','GradP','CGP','BFGSP','HALS','MUR')
xlabel('Number of iterations')
ylabel('Relative errors')
%title('')

figure(2)
X = categorical({'ODE','ANLS','Proco-ALS','CCG','GradP','CGP','BFGSP','HALS','MUR'});
X = reordercats(X,{'ODE','ANLS','Proco-ALS','CCG','GradP','CGP','BFGSP','HALS','MUR'});
Y = [time_0 time_1 time_2 time_3 time_4 time_5 time_6 time_7 time_7];
bar(X,Y)
ylabel('Running time (Second)')
%%
function dydt=ODE(t,xx,X,epsilon,R)
%%Tensor decomposition based on the colloborative neurodynamic optimization
%% size of dynamic system
% [n,R]=size(A_1);
szX = size(X);
dydt=zeros(sum(szX)*R,1);
B = vec2fac(xx,szX);

%% epsilon

% ALS
alg = 'als2';
%alg = 'hals2';
switch alg
    case 'als'

        epsilon_1=epsilon.eps_1;
        epsilon_2=epsilon.eps_2;
        epsilon_3=epsilon.eps_3;
        x3tx3 = B{3}'*B{3};
        x2tx2 = B{2}*B{2};
        x1tx1 = B{1}'*B{1};
        H1 = ((x3tx3).*((x2tx2)));% Hessian for B1
        H2 = ((x3tx3).*((x1tx1)));% Hessian for B2
        H3 = ((x2tx2).*((x1tx1))); % Hessian for B3

        tX = tensor(X);
        tX1 = mttkrp(tX,B,1);
        tX2 = mttkrp(tX,B,2);
        tX3 = mttkrp(tX,B,3);

        x1_als = max(tX1/H1,0);% x1als = max(0,Y1*khatr(B3,B2)*inv(H1))
        x2_als = max(tX2/H2,0);
        x3_als = max(tX3/H3,0);

        dydt = [(1/epsilon_1)*(x1_als(:)-B{1}(:))
                (1/epsilon_2)*(x2_als(:)-B{2}(:))
                (1/epsilon_3)*(x3_als(:)-B{3}(:))];

    case 'als2'
        tX = tensor(X);
        Q = zeros(R,R,numel(B));
        for k = 1:numel(B)
            Q(:,:,k) = B{k}'*B{k};
        end

        % ALS update
        Bnew = B;
        nu = 1e-6;
        for k = 1:numel(B)
            Qk = prod(Q(:,:,[1:k-1 k+1:end]),3);
            Bnew{k} = mttkrp(tX,Bnew,k)/(Qk+nu*eye(R));
            Bnew{k} = max(1e-8,Bnew{k});
            Q(:,:,k) = Bnew{k}'*Bnew{k};
        end

%         Yx = normalize(ktensor(Bnew));
%         Bnew = Yx.U;
%         for k = 1:numel(Bnew)
%             Bnew{k} = Bnew{k}*diag(Yx.lambda.^(1/numel(Bnew)));
%         end

        dB = B;
        for k = 1:numel(B)
            dB{k} = (Bnew{k}-B{k})/epsilon.(sprintf('eps_%d',k));
        end
        dydt = fac2vec(dB);

    case 'hals2'

        tX = tensor(X);
        B = vec2fac(xx,size(X));
        Q = zeros(R,R,numel(B));
        for k = 1:numel(B)
            Q(:,:,k) = B{k}'*B{k};
        end

        % ALS update
        Bnew = B;
        for k = 1:numel(B)
            Qk = prod(Q(:,:,[1:k-1 k+1:end]),3);
            for r = 1:size(B{k},2)
                br = cellfun(@(x) x(:,r),Bnew([1:k-1 k+1:end]),'uni',0);
                Bnew{k}(:,r) = (ttv(tX,br,-k) - Bnew{k}(:,[1:r-1 r+1:end])*Qk([1:r-1 r+1:end],r))/Qk(r,r);
                Bnew{k}(:,r) = max(1e-8,Bnew{k}(:,r));
            end
            Q(:,:,k) = Bnew{k}'*Bnew{k};
        end

        dB = B;
        for k = 1:numel(B)
            dB{k} = (Bnew{k}-B{k})/epsilon.(sprintf('eps_%d',k));
        end
        dydt = fac2vec(dB);


    case 'hals'

        %B = {xlda(:)' x1 x2 x3};
        opts = ncp_hals;
        opts.init = B;
        opts.maxiters = 1;
        opts.tol = 1e-10;
        [Yx,out] = ncp_hals(tensor(X,[1 size(X)]),R,opts);
        %[Yx,out] = ncp_fLM(tensor(X),R,opts);

        Bnew = Yx.U;
        Bnew{1} = Bnew{1}*diag(Yx.lambda);
        dB = B;
        for k = 1:numel(Bnew)
            %Bnew{k} = Bnew{k}*diag(Yx.lambda.^(1/3));
            dB{k} = (Bnew{k}-B{k})/epsilon.eps_1;
        end
        dB = dB([2:end 1]);
        dydt = fac2vec(dB);

%     case 'opt' 
%         opts = struct();
%         opts.init = B;
%         opts.maxiters = 1;
%         opts.tol = 1e-10;
%         params = ncg('defaults');
%         params.MaxIters = 2;
%         [Yx,out] = cp_opt(tensor(X),R,'init',B,'alg_options',params);
%         %[Yx,out] = ncp_fLM(tensor(X),R,opts);
% 
%         Bnew = Yx.U;
%         dB = B;
%         for k = 1:3
%             Bnew{k} = Bnew{k}*diag(Yx.lambda.^(1/3));
%             %dB{k} = (Bnew{k}-B{k})/epsilon.(sprintf('eps_%d',k));
%             dB{k} = (Bnew{k}-B{k})/epsilon.eps_1;
%         end
%         dydt = fac2vec(dB);

end




end


function X = get_mnist(digits)

addpath('C:\Users\huyph\Documents\MATLAB\MNIST\')

% clear all
% digits = [1 7 ];

No_digits_ = 1000; % number of images per digit
imageSize =[28 28];
orientationsPerScale = [8 8 8 8]; % assume number of orientations at scales are the same
numberBlocks=imageSize(1);% number of blocks after downsampling

%%
clear X
newsize = [14 14];
true_labels  = [];
No_digits_2 = 100;
for kd = 1:numel(digits)
    load(sprintf('mnist_gabor_no%d_%d.mat',digits(kd),No_digits_));
    F = F(:,:,:,:,1:No_digits_2);
    szF = size(F);

    F = reshape(F,imageSize(1),imageSize(2),[]);
    F = imresize(F,newsize);
    newszF = szF;
    newszF(1:2) = newsize;
    F = reshape(F,newszF);
    X(:,kd) = F(:);

    true_labels = [true_labels  ; kd*ones(No_digits_,1)];
end
No_digits = size(F,5);
SzF = size(F);

No_digits_ = 1000; % number of images per digit
imageSize = newsize;
orientationsPerScale = [8 8 8 8]; % assume number of orientations at scales are the same
numberBlocks=imageSize(1);% number of blocks after downsampling

X = reshape(X,[SzF(1)*SzF(2) SzF(3) SzF(4) SzF(5)*numel(digits)]);
% centralize X
N = ndims(X);
%Xm = mean(X,N);
%X = bsxfun(@minus,X,Xm);

y = true_labels;
X0 = X;

%% take three views
% X = squeeze(X0(:,1:2:5,1,:));
end