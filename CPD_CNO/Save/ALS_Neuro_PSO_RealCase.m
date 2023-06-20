clc;clear all
%% ----------------------------------------------------
% Data loading
% -----------------------------------------------------
% load('ORL_32x32.mat');
% load('COIL20.mat');   %1440 R=10
% dim = [1440,32,32];
% AA=permute(reshape(fea,dim),[2,3,1]);
% K=1440;
load('V.mat');
AA = V;
clear V;
K=180;

X=double(AA(:,:,1:K));
X=X/max(X(:));
R=12;
dim(1) = K;
tim = tic 
%% ----------------------------------------------------
% CNO initi.
% -----------------------------------------------------
NN=6;
RR=R;
tspan = [0 10*1e-3];
sizeX = size(X);
dim = [sizeX(1),sizeX(1),sizeX(2),sizeX(2),sizeX(3),sizeX(3)];
for i=1:NN
    V{i}= rand(1,sum(sizeX*R));
end

for i=1:NN
%     B{i,1}=rand(sizeX(1),RR);
%     B{i,2}=rand(sizeX(2),RR);
%     B{i,3}=rand(sizeX(3),RR);
    opts = cp_init();
    opts.init = {'orth' 'fiber' 'rand' };
    Ui = cp_init(tensor(X),RR,opts);
    B{i,1} = Ui{1};
    B{i,2} = Ui{2};
    B{i,3} = Ui{3};
    y0{i,1}=[B{i,1}(:);B{i,2}(:);B{i,3}(:)];
end

B_prev = B;

%%
nu=2*ones(NN,1); tau=1;
N = 3;
for i=1:NN
    ell2 = zeros(N,RR);
    for n = 1:N
        ell2(n,:) = sum(abs(B{i,n}).^2);
    end
    mm = zeros(N,R);
    for n = 1:N
        mm(n,:) = prod(ell2([1:n-1 n+1:N],:),1);
    end
    mu(i)=min(1,tau*max(mm(:)));
    
end
mu0 = mu;
warning('off');
%% ----------------------------------------------------
% Solve set of ODE by CNO
% -----------------------------------------------------
% Main optimization loop
for j=1:50
    % Loop on NN
    for i=1:NN
        epsilon.eps_1=10^(-4);
        epsilon.eps_2=10^(-4);
        epsilon.eps_3=10^(-4);
        %%-Call of ODE45 solver
        [t{i,j},y{i,j}] = ode45(@(t,xx) ODE(t,xx,X,epsilon,dim(i),RR,mu(i)), tspan, y0{i,j});

        z{i}=y{i,j};
        % zz{k}=y{i}(end,:);
        %%Reinitialization
        % y0=y(end,:);
        [~,nn]=size(y{i,j});
        %%Updating the factors
        idx1 = sizeX(1)*R;
        B{i,1}=reshape(z{i}(end,1:idx1),[sizeX(1),RR]);
        idx2 = sizeX(1)*R+sizeX(2)*R;
        B{i,2}=reshape(z{i}(end,idx1+1:idx2),[sizeX(2),RR]);
        B{i,3}=reshape(z{i}(end,idx2+1:nn),[sizeX(3),RR]);
        Error=X-double(full(ktensor(ones(RR,1),B{i,1},B{i,2},B{i,3})));
        % Display erros
        ee(i,j)=norm(Error(:))/norm(X(:))
        Grad{i} = gradientCPD(B(i,:),X,RR);


%         [~,nn]=size(y{i,j});
%         % figure(2)
%         % subplot(3,1,1)
%         figure(1)
%         plot(t{i,j},y{i,j}(:,1:idx1),'LineWidth',1.5)
%         % drawnow
%         % title('Transient behaviors of CNO')
%         xlabel('Time (t)','fontweight','bold','FontSize',15)
%         ylabel('A^{(1)}','fontweight','bold','FontSize',15)
%         % figure(3)
%         % subplot(3,1,2)
%         figure(2)
%         plot(t{i,j},y{i,j}(:,idx1+1:idx2),'LineWidth',1.5)
%         % drawnow
%         % title('Transient behaviors of CNO')
%         xlabel('Time (t)','fontweight','bold','FontSize',15)
%         ylabel('A^{(2)}','fontweight','bold','FontSize',15)
%         % figure(4)
%         % subplot(3,1,3)
%         figure(3)
%         plot(t{i,j},y{i,j}(:,idx2+1:nn),'LineWidth',1.5)
%         % drawnow
%         % title('Transient behaviors of CNO')
%         xlabel('Time (t)','fontweight','bold','FontSize',15)
%         ylabel('A^{(3)}','fontweight','bold','FontSize',15)
%         keyboard
%         close all



    end
    % Call of PSO method
    tt_h=min(min(ee));
    [tt_1,tt_2]=find(ee==tt_h);
    GH=y{tt_1(1),tt_2(1)}(end,:);
    w=0;
    for i=1:NN
        [ii,jj]=min(ee(i,:));
        [y0{i,j+1},V{i}]=PSO_Code(V{i},z{i}(end,:),y{i,jj}(end,:),GH);
        w=w+norm(y{i,jj}(end,:)-GH);
    end
    % Display PSO advancement
    DI=w/NN

end
toc(tim)
%%
% [nn,~]=size(y);
% subplot(1,3,1)
% plot(t(1:nn/3,:),y(1:nn/3,:))
% drawnow
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A_1')
% subplot(1,3,2)
% plot(t(nn/3+1:(2*nn)/3,:),y(nn/3+1:(2*nn)/3,:))
% drawnow
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A_2')
% subplot(1,3,3)
% plot(t((2*nn)/3+1:nn,:),y((2*nn)/3+1:nn,:))
% drawnow
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A_3')

%end
% [~,nn]=size(y);
% figure(2)
% subplot(2,1,1)
% plot(t,y(:,1:nn),'LineWidth',1.5)
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(1)}')
% figure(3)
% subplot(2,1,1)
% plot(t,y(:,nn/3+1:(2*nn)/3),'LineWidth',1.5)
%  title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(2)}')
% figure(4)
% subplot(2,1,1)
% plot(t,y(:,(2*nn)/3+1:nn),'LineWidth',1.5)
%  title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A^{(3)}')
% sgtitle('Transient behaviors of CNO')
% [t,y] = ode45(@ODE,[0 20],ones(30,1));
% [nn,~]=size(y);
% figure(1)
% plot(t(1:nn/3,:),y(1:nn/3,:),'LineWidth',1.5)
% xlim([0 0.5])
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A_1')
% figure(2)
% plot(t(nn/3+1:(2*nn)/3,:),y(nn/3+1:(2*nn)/3,:),'LineWidth',1.5)
% xlim([0 0.5])
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A_2')
% figure(3)
% plot(t((2*nn)/3+1:nn,:),y((2*nn)/3+1:nn,:),'LineWidth',1.5)
% xlim([0 0.5])
% title('Transient behaviors of CNO')
% xlabel('Time (t)')
% ylabel('A_3')
%%
function [dydt,G]=ODE(t,xx,X,epsilon,n,R,mu)
%%Tensor decomposition based on the colloborative neurodynamic optimization
    %% size of dynamic system
    sizeX = size(X);
    dydt=zeros(sum(sizeX*R),1);
    nn=sum(sizeX*R);
    idx1 = sizeX(1)*R;
    x1=reshape(xx(1:idx1),[sizeX(1),R]);
    idx2 = sizeX(1)*R+sizeX(2)*R;
    x2=reshape(xx(idx1+1:idx2),[sizeX(2),R]);
    x3=reshape(xx(idx2+1:nn),[sizeX(3),R]);
    %% epsilon
    epsilon_1=epsilon.eps_1;
    epsilon_2=epsilon.eps_2;
    epsilon_3=epsilon.eps_3;
    %% Gradient computation
    x3tx3 = x3'*x3;
    x2tx2 = x2'*x2;
    x1tx1 = x1'*x1;
    H1 = ((x3tx3).*((x2tx2)));
    tX = tensor(X);
    gradient_1=x1*H1-mttkrp(tX,{x1,x2,x3},1);
    H2 = ((x3tx3).*((x1tx1)));
    gradient_2=x2*H2-mttkrp(tX,{x1,x2,x3},2);
    H3 = ((x2tx2).*((x1tx1)));
    gradient_3=x3*H3-mttkrp(tX,{x1,x2,x3},3);
    
    %% Neurodynamic model
    %---------
    % Factor 1
    % --------
    % Preconditioning ALS approach
    delta = mu;
    H=H1+delta*eye(R);
    v = max(x1 - gradient_1/H,0);
    % v = SimplexProj(x1 - gradient_1/H);
    
    % Cubic regularized approach
    % [alphak,step] = step_CubReg_Newton(H,gradient_1);
    % v = max(x1 - alphak*step,0);
    
    % computing dydt for ODE solver
    dydt(1:idx1)=(1/epsilon_1)*(-x1(:)+v(:));
    
    %---------
    % Factor 2
    % --------
    % Preconditioning ALS approach
    H=H2+delta*eye(R);
    v = max(x2 - gradient_2/H,0);
    % v = SimplexProj(x2 - gradient_2/H);
    
    % Cubic regularized approach
    % [alphak,step] = step_CubReg_Newton(H,gradient_2);
    % v = max(x2 - alphak*step,0);
    
    % computing dydt for ODE solver
    dydt(idx1+1:idx2)=(1/epsilon_2)*(-x2(:)+ v(:));
    
    %---------
    % Factor 3
    % --------
    % Preconditioning ALS approach
    H=H3+delta*eye(R);
    v = max(x3 - gradient_3/H,0);
    % v = SimplexProj(x3 - gradient_3/H);
    
    % Cubic regularized approach
    % [alphak,step] = step_CubReg_Newton(H,gradient_3);
    % v = max(x3 - alphak*step,0);
    
    % computing dydt for ODE solver
    dydt(idx2+1:nn)=(1/epsilon_3)*(-x3(:)+v(:));

end

function Grad = gradientCPD(xx,X,R)
    %% size of dynamic system
    % [n,R]=size(A_1);
    sizeX = size(X);
    x1 = xx{1};
    x2 = xx{2};
    x3 = xx{3};
    
    %% Gradient computation
    x3tx3 = x3'*x3;
    x2tx2 = x2'*x2;
    x1tx1 = x1'*x1;
    H1 = ((x3tx3).*((x2tx2)));
    tX = tensor(X);
    gradient_1=x1*H1-mttkrp(tX,{x1,x2,x3},1);
    H2 = ((x3tx3).*((x1tx1)));
    gradient_2=x2*H2-mttkrp(tX,{x1,x2,x3},2);
    H3 = ((x2tx2).*((x1tx1)));
    gradient_3=x3*H3-mttkrp(tX,{x1,x2,x3},3);
    Grad = [vec(gradient_1);vec(gradient_2);vec(gradient_3)];
end
%% --------------
% Data visu.
%----------------
%ORL
% species=cell(K,1);
% species(1:10)={num2str(1)};
% j=10;
% for i=1:39
%     species(1+j:j+10)={num2str(i+1)};
%     j=(i+1)*10;
% end

% COIL
% species=cell(K,1);
% species(1:72)={num2str(1)};
% j=72;
% for i=1:19
%     species(1+j:j+72)={num2str(i+1)};
%     j=(i+1)*72;
% end
% 
% Y = tsne(B{end,3});
% gscatter(Y(:,1),Y(:,2),species)
