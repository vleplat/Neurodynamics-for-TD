clc;clear all
%% ---------------------------------------------------------
% Data init.
% ----------------------------------------------------------
dataSetSel = 6; % Select data set 6
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
        % AA=permute(reshape(fea,[165,32,32]),[2,3,1]);
        
    
    case 5
        disp('cifar10 data set selected...')
        load('cifar.mat');
        AA=permute(reshape(data,[50000,32,32*3]),[2,3,1]);
        K=1000;

    case 6
        disp('Cuprite HSI selected...')
        load('V.mat');
        % AA=permute(reshape(data,[50000,32,32*3]),[2,3,1]);
        AA = V;
%         disp('Synthetic HSI selected...')
%         load ('SyntheticHSI_R12.mat')
%         AA=permute(reshape(syntheticImage,[180,100,100]),[2,3,1]);

        K=180;

    case 7
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

%%% facto rank and LS parameters
R=12; %11 for Yale, 12 for Cuprite
beta=.2;
alpha=.2;

%%% Init. for the factors
% Random
sizeX = size(X);
A_10 = rand(sizeX(1),R);
A_20 = rand(sizeX(2),R);
A_30 = rand(K,R);
A_1{1} = A_10;
A_2{1} = A_20;
A_3{1} = A_30;

% Using mixed strategies (from Anh Huy Phan)
% opts = cp_init();
% opts.init = {'rand' 'fiber' 'fiber' }; %rand fiber fiber
% Ui = cp_init(tensor(X),R,opts);
% A_10 = Ui{1};
% A_20 = Ui{2};
% A_30 = Ui{3};
% A_1{1} = A_10;
% A_2{1} = A_20;
% A_3{1} = A_30;

%%% Errors arrays
err=[];
X_app=(double(full(ktensor(ones(R,1),A_1{1},A_2{1},A_3{1}))));
Error=X-X_app;
relaerror=norm(Error(:))/norm(X(:));
err=[err,relaerror];
%% ---------------------------------------------------------
% Main Optimization Loop
% ----------------------------------------------------------
delta = 0;
iter = 500;
t=1:iter;

for i=1:iter
    %%% Restart values for parameters lambda_i
    lambda_1=1;
    lambda_2=1;
    lambda_3=1;
   

    %%% Update of factor one
    x3tx3 = A_3{i}'*A_3{i};
    x2tx2 = A_2{i}'*A_2{i};
    H1 = ((x3tx3).*((x2tx2)));
    H=H1+delta*eye(R);
    gradient_1=A_1{i}*H1-mttkrp(tX,{A_1{i},A_2{i},A_3{i}},1);
    
    
    A_1{i+1}=(A_1{i}+lambda_1*(+max(A_1{i}-gradient_1/H,0)))/(1+lambda_1);

    % Line-search
    Error_11=X-double(full(ktensor(ones(R,1),A_1{i},A_2{i},A_3{i})));
    Error_12=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
    
    if (norm(Error_12(:))^2-norm(Error_11(:))^2)>=alpha*lambda_1*gradient_1(:)'*(A_1{i+1}(:)-A_1{i}(:))
    while (norm(Error_12(:))^2-norm(Error_11(:))^2)>alpha*lambda_1*gradient_1(:)'*(A_1{i+1}(:)-A_1{i}(:))
        lambda_1=lambda_1*beta;
        A_1{i+1}=(A_1{i}+lambda_1*(+max(A_1{i}-gradient_1/H,0)))/(1+lambda_1);
        Error_11=X-double(full(ktensor(ones(R,1),A_1{i},A_2{i},A_3{i})));
        Error_12=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
        if lambda_1<1e-9
            break
        end
    end
    end
    
    %%% Update of factor two
    x1tx1 = A_1{i+1}'*A_1{i+1};
    H2 = ((x3tx3).*((x1tx1)));
    H  = H2+delta*eye(R);
    gradient_2 = A_2{i}*H2-mttkrp(tX,{A_1{i+1},A_2{i},A_3{i}},2);
    A_2{i+1}=(A_2{i}+lambda_2*(+max(A_2{i}-gradient_2/H,0)))/(1+lambda_2);

    % Line-search
    Error_21=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
    Error_22=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i})));
    if (norm(Error_22(:))^2-norm(Error_21(:))^2)>=alpha*lambda_2*gradient_2(:)'*(A_2{i+1}(:)-A_2{i}(:))
    while (norm(Error_22(:))^2-norm(Error_21(:))^2)>alpha*lambda_2*gradient_2(:)'*(A_2{i+1}(:)-A_2{i}(:))
        lambda_2=lambda_2*beta;
        A_2{i+1}=(A_2{i}+lambda_2*(+max(A_2{i}-gradient_2/H,0)))/(1+lambda_2);
        Error_21=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
        Error_22=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i})));
        if lambda_2<1e-9
            break
        end
    end
    end
  
    %%% Update of factor three
    x2tx2 = A_2{i+1}'*A_2{i+1};
    H3 = ((x2tx2).*((x1tx1)));
    H = H3+delta*eye(R);
    gradient_3 = A_3{i}*H3-mttkrp(tX,{A_1{i+1},A_2{i+1},A_3{i}},3);
    A_3{i+1}=(A_3{i}+lambda_3*(max(A_3{i}-gradient_3/H,0)))/(1+lambda_3);

    % Line-search
    Error_31=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i})));
    Error_32=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i+1})));
    if (norm(Error_32(:))^2-norm(Error_31(:))^2)>=alpha*lambda_3*gradient_3(:)'*(A_3{i+1}(:)-A_3{i}(:))
    while (norm(Error_32(:))^2-norm(Error_31(:))^2)>alpha*lambda_3*gradient_3(:)'*(A_3{i+1}(:)-A_3{i}(:))
        lambda_3=lambda_3*beta;
        A_3{i+1}=(A_3{i}+lambda_3*(max(A_3{i}-gradient_3/H,0)))/(1+lambda_3);
        Error_31=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i})));
        Error_32=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i+1})));
        if lambda_3<1e-9
            break
        end
    end
    end

    %%% Errors computations
    Err=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i+1})));
    ET=norm(Err(:))/norm(X(:));
    err=[err,ET];

    %%% Monitoring
    lamTab = [lambda_1 lambda_2 lambda_3];
    if i > 1 && mod(i,50)==0
        fprintf('Iter: %d   | Rel.Error  %.4d  | min(lam): %.4d  |max(lam): %d \n',i,err(end), min(lamTab),max(lamTab))
    end
end
%% ---------------------------------------------------------
% Post-processing
% ----------------------------------------------------------

%%% Clustering
switch dataSetSel
    case 1
        %ORL
        species=cell(K,1);
        species(1:10)={num2str(1)};
        j=10;
        for i=1:39
            species(1+j:j+10)={num2str(i+1)};
            j=(i+1)*10;
        end

    case 2
        %YALE
        species=cell(K,1);
        species(1:15)={num2str(1)};
        j=15;
        for i=1:10
            species(1+j:j+15)={num2str(i+1)};
            j=(i+1)*15;
        end

    case 3
        %COIL20
        species=cell(K,1);
        species(1:72)={num2str(1)};
        j=72;
        for i=1:19
            species(1+j:j+72)={num2str(i+1)};
            j=(i+1)*72;
        end
    
    case 4
        %MNIST
        species=cell(K,1);
        species(1:140)={num2str(1)};
        j=140;
        for i=1:9
            species(1+j:j+140)={num2str(i+1)};
            j=(i+1)*140;
        end
    
    case 6
        %HSI-Cuprite
        species=cell(K,1);
        species(1:15)={num2str(1)};
        j=15;
        for i=1:11
            species(1+j:j+15)={num2str(i+1)};
            j=(i+1)*15;
        end    
    
    otherwise
        disp('wrong selection for the data set...')
        return
end


% species=cell(K,1);
% species(1:11)={num2str(1)};
% j=11;
% for i=1:14
%     species(1+j:j+11)={num2str(i+1)};
%     j=(i+1)*11;
% end

%%% Display Clustering
clr = hsv(R);
Y = tsne(A_3{end});
figure;
gscatter(Y(:,1),Y(:,2),species,clr)



%%% Congergence plots
figure;
semilogy(err(3:end),'LineWidth',1.5)
legend('CNO-CPD')

figure(2)
for i=1:iter
    Ahh1(i,:)=A_1{i}(:)';
end
plot(t,Ahh1(:,1:end),'LineWidth',1.5)
title('Transient behaviors of CNO')
xlabel('Time (t)')
ylabel('A^{(1)}')

figure(3)
for i=1:iter
    Ahh2(i,:)=A_2{i}(:)';
end
plot(t,Ahh2(:,1:end),'LineWidth',1.5)
title('Transient behaviors of CNO')
xlabel('Time (t)')
ylabel('A^{(2)}')

figure(4)
for i=1:iter
    Ahh3(i,:)=A_3{i}(:)';
end
plot(t,Ahh3(:,1:end),'LineWidth',1.5)
title('Transient behaviors of CNO')
xlabel('Time (t)')
ylabel('A^{(3)}')

%% cuprite 
figure(5)
for r=1:R

    subplot(3,4,r); 
    plot(A_3{i+1}(:,r))
end

%% Comparison vs HALS
opts = ncp_hals;
opts.init = {'rand' 'rand' 'rand' };
opts.maxiters = 500;
[P,output] = ncp_hals(tX,R,opts);
figure(1);clf; semilogy(output.Fit(:,1),1-output.Fit(:,2));hold on;
semilogy(err(1:end));legend('cpd HALS','CNO')
xlabel('Iterations'); ylabel('Relative Error')

%%% Display Clustering
clr = hsv(R);
Yhals = tsne(P{end});
figure;
gscatter(Yhals(:,1),Yhals(:,2),species,clr)
