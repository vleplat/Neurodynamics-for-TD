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
% X = X*5;
% X = max(X,eps);
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

for j=1:50
%%

    for i=1:NN
        epsilon.eps_1=10^(-4);
        epsilon.eps_2=10^(-4);
        epsilon.eps_3=10^(-4);
        %for i=1:1

        %%
% tspan = [0 2*1e-3];
% option.MaxStep=1e-3;
% option.InitialStep=1e-3;
      %  [t{i,j},y{i,j}] = ode45(@(t,xx) ODE(t,xx,X,epsilon,n,RR), tspan, y0{i,j},option);
        [t{i,j},y{i,j}] = ode45(@(t,xx) ODE(t,xx,X,epsilon,dim(i),RR,mu(i)), tspan, y0{i,j});
        %%
        
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
        %%Orthonormalize the factors
        % A_1(:,i)=A_1(:,i)/norm(A_1(:,i));
        % A_2(:,i)=A_2(:,i)/norm(A_2(:,i));
        % A_3(:,i)=A_3(:,i)/norm(A_3(:,i));
        Error=X-double(full(ktensor(ones(RR,1),B{i,1},B{i,2},B{i,3})));
        ee(i,j)=norm(Error(:))/norm(X(:))
        % ee=[e{i,j},ee];
        % norm(Error(:))/norm(X(:))
        % norm(y{i}(end,:)-X(:))
        % k=k+1;

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
    tt_h=min(min(ee));
    [tt_1,tt_2]=find(ee==tt_h);
    GH=y{tt_1(1),tt_2(1)}(end,:);
    w=0;
    for i=1:NN
        [ii,jj]=min(ee(i,:));
        [y0{i,j+1},V{i}]=PSO_Code(V{i},z{i}(end,:),y{i,jj}(end,:),GH);
        w=w+norm(y{i,jj}(end,:)-GH);
    end
    DI=w/NN

    %% update rho
    % Rho to compute the damp parameter
    % mu = mu0/sqrt(j);
    
    %% update damp parameter
%     if j>1
%         for i=1:NN
%              Rr = Grad{i};
%              err= ee(i,j-1)*norm(X(:));
%              err2 = ee(i,j)*norm(X(:));
%              d = cellfun(@(x,y) x(:)-y(:),B(i,:),B_prev(i,:),'uni',0);
%              d = cell2mat(d(:));
%              U = B(i,:);
%              rho=real((err-err2)/(d(:)'*(Rr+mu(i)*d(:))));   
%               
%             if err2>err               %%% step is not accepted
%                 mu(i)=mu(i)*nu(i); nu(i)=2*nu(i);                                      % Eq. (5.7)
%                 B(i,:) = B_prev(i,:);
%             else
%                 nu(i)=2;
%                 mu(i)=mu(i)*max([1/3 1-(2*rho-1)^3]);                         % Eq. (5.7)
%                  
%             end
%             mu(i)=min(200,mu(i));
%         end
        
%   end

    B_prev = B;
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
% [n,R]=size(A_1);
sizeX = size(X);
% dydt=zeros(3*n*R,1);
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
% gradient_1=x1*((x3'*x3).*((x2'*x2)))-unfold(X,1)*(khatrirao_fast(x3,x2));
% gradient_2=x2*((x3'*x3).*((x1'*x1)))-unfold(X,2)*(khatrirao_fast(x3,x1));
% gradient_3=x3*((x2'*x2).*((x1'*x1)))-unfold(X,3)*(khatrirao_fast(x2,x1));
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
% v = SimplexProj(x1 - gradient_1);
delta = mu;
[k,~] = size(x1);
% dIp = sparse(delta*eye(k*R));
%H = sparse(kron(H1',eye(k))+delta*eye(k*R));% + dIp;
H=H1+delta*eye(R);
% [alphak,step] = step_CubReg_Newton(H,gradient_1);
v = max(x1 - gradient_1/H,0);
% v = max(x1 - alphak*step,0);

%v = max(x1(:) - H\gradient_1(:),0);
dydt(1:idx1)=(1/epsilon_1)*(-x1(:)+v(:));
% v = SimplexProj(x2 - gradient_2);
[k,~] = size(x2);
% dIp = sparse(delta*eye(k*R));
%H = sparse(kron(H2',eye(k))+delta*eye(k*R));% + dIp;
%v = max(x2(:) - H\gradient_2(:),0);
H=H2+delta*eye(R);
% [alphak,step] = step_CubReg_Newton(H,gradient_2);
v = max(x2 - gradient_2/H,0);
% v = max(x2 - alphak*step,0);

dydt(idx1+1:idx2)=(1/epsilon_2)*(-x2(:)+ v(:));
[k,~] = size(x3);
% dIp = sparse(delta*eye(k*R));
%H = sparse(kron(H3',eye(k))+delta*eye(k*R));% + dIp;
H=H3+delta*eye(R);
% [alphak,step] = step_CubReg_Newton(H,gradient_3);
v = max(x3 - gradient_3/H,0);
% v = max(x3 - alphak*step,0);

%v = max(x3(:) - H\gradient_3(:),0);
% v = SimplexProj(x3-gradient_3);
dydt(idx2+1:nn)=(1/epsilon_3)*(-x3(:)+v(:));

end

function Grad = gradientCPD(xx,X,R)
    %% size of dynamic system
% [n,R]=size(A_1);
sizeX = size(X);

nn=sum(sizeX*R);
idx1 = sizeX(1)*R;
x1 = xx{1};
x2 = xx{2};
x3 = xx{3};

%% Gradient computation
% gradient_1=x1*((x3'*x3).*((x2'*x2)))-unfold(X,1)*(khatrirao_fast(x3,x2));
% gradient_2=x2*((x3'*x3).*((x1'*x1)))-unfold(X,2)*(khatrirao_fast(x3,x1));
% gradient_3=x3*((x2'*x2).*((x1'*x1)))-unfold(X,3)*(khatrirao_fast(x2,x1));
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
