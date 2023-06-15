%% ---------------------------------------------------------
% Data init.
% ----------------------------------------------------------
clc;clear all
n=100;  
R=10;    
tim = tic ;

%%% number of RNN
NN=6;
RR=R;
epsilon.eps_1=10^(-4);
epsilon.eps_2=10^(-4);
epsilon.eps_3=10^(-4);
maxTimestep = 15;

%%% Init for ODE45
tspan = [0 2*1e-3];
for i=1:NN
    V{i}= rand(1,3*n*RR);
end

%%% Generating a nonnegative tensor
A{1}=rand(n,R);
A{2}=rand(n,R);
A{3}=rand(n,R);
X=double(full(ktensor(ones(R,1),A{1},A{2},A{3})));

%%% Initi of the factors (random)
for i=1:NN
    B{i,1}=rand(n,RR);
    B{i,2}=rand(n,RR);
    B{i,3}=rand(n,RR);
    y0{i,1}=[B{i,1}(:);B{i,2}(:);B{i,3}(:)];
end

%%% error matrix
ee=[];

%% ---------------------------------------------------------
% Main Optimization Loop
% ----------------------------------------------------------

% time loop
for j=1:maxTimestep
   
    % loop over RNNs
    for i=1:NN 
        
        %%Call of ODE solver
        % tspan = [0 2*1e-3];
        % option.MaxStep=1e-3;
        % option.InitialStep=1e-3;
        % [t{i,j},y{i,j}] = ode45(@(t,xx) ODE(t,xx,X,epsilon,n,RR), tspan, y0{i,j},option);
        [t{i,j},y{i,j}] = ode45(@(t,xx) ODE(t,xx,X,epsilon,n,RR), tspan, y0{i,j});
        z{i}=y{i,j};
        [~,nn]=size(y{i,j});

        %%Updating the factors
        B{i,1}=reshape(z{i}(end,1:nn/3),[n,RR]);
        B{i,2}=reshape(z{i}(end,nn/3+1:(2*nn)/3),[n,RR]);
        B{i,3}=reshape(z{i}(end,(2*nn)/3+1:nn),[n,RR]);
        
        %%Errors computation
        Error=X-double(full(ktensor(ones(RR,1),B{i,1},B{i,2},B{i,3})));
        ee(i,j)=norm(Error(:))/norm(X(:))


%         [~,nn]=size(y{i,j});
%         % figure(2)
%         % subplot(3,1,1)
%         figure(1)
%         plot(t{i,j},y{i,j}(:,1:nn),'LineWidth',1.5)
%         % drawnow
%         % title('Transient behaviors of CNO')
%         xlabel('Time (t)','fontweight','bold','FontSize',15)
%         ylabel('A^{(1)}','fontweight','bold','FontSize',15)
%         % figure(3)
%         % subplot(3,1,2)
%         figure(2)
%         plot(t{i,j},y{i,j}(:,nn/3+1:(2*nn)/3),'LineWidth',1.5)
%         % drawnow
%         % title('Transient behaviors of CNO')
%         xlabel('Time (t)','fontweight','bold','FontSize',15)
%         ylabel('A^{(2)}','fontweight','bold','FontSize',15)
%         % figure(4)
%         % subplot(3,1,3)
%         figure(3)
%         plot(t{i,j},y{i,j}(:,(2*nn)/3+1:nn),'LineWidth',1.5)
%         % drawnow
%         % title('Transient behaviors of CNO')
%         xlabel('Time (t)','fontweight','bold','FontSize',15)
%         ylabel('A^{(3)}','fontweight','bold','FontSize',15)



    end
    
    %%Call of PSO - communication steps between the RNNs
    tt_h=min(min(ee));
    [tt_1,tt_2]=find(ee==tt_h);
    GH=y{tt_1(1),tt_2(1)}(end,:);
    w=0;
    for i=1:NN
        [ii,jj]=min(ee(i,:));
        [y0{i,j+1},V{i}]=PSO_Code(V{i},z{i}(end,:),y{i,jj}(end,:),GH);
        w=w+norm(y{i,jj}(end,:)-GH);
    end
    %%Display divergence
    DI=w/NN

end
toc(tim)
%% ---------------------------------------------------------
% Post-processing
% ----------------------------------------------------------
[~,nn]=size(y{i,j});

%%%Dynamic of the systems over time
subplot(3,1,1)
plot(t{i,j},y{i,j}(:,1:nn),'LineWidth',1.5)
title('Transient behaviors of CNO - factor one')
xlabel('Time (t)','fontweight','bold','FontSize',15)
ylabel('A^{(1)}','fontweight','bold','FontSize',15)

subplot(3,1,2)
plot(t{i,j},y{i,j}(:,nn/3+1:(2*nn)/3),'LineWidth',1.5)
title('Transient behaviors of CNO - factor two')
xlabel('Time (t)','fontweight','bold','FontSize',15)
ylabel('A^{(2)}','fontweight','bold','FontSize',15)

subplot(3,1,3)
plot(t{i,j},y{i,j}(:,(2*nn)/3+1:nn),'LineWidth',1.5)
title('Transient behaviors of CNO - factor three')
xlabel('Time (t)','fontweight','bold','FontSize',15)
ylabel('A^{(3)}','fontweight','bold','FontSize',15)
%% ---------------------------------------------------------
% Library of functions
% ----------------------------------------------------------
function dydt=ODE(t,xx,X,epsilon,n,R)
%%Tensor decomposition based on the colloborative neurodynamic optimization
%%size of dynamic system
    dydt=zeros(3*n*R,1);
    nn=3*n*R;
    x1=reshape(xx(1:nn/3),[n,R]);
    x2=reshape(xx(nn/3+1:(2*nn)/3),[n,R]);
    x3=reshape(xx((2*nn)/3+1:nn),[n,R]);
    
    %%epsilon
    epsilon_1=epsilon.eps_1;
    epsilon_2=epsilon.eps_2;
    epsilon_3=epsilon.eps_3;
    
    %%Gradient computation
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
    
    %%Neurodynamic model
    delta = 0;
    % v = SimplexProj(x1 - gradient_1);
    H = H1+delta*eye(R);
    v = max(x1 - gradient_1/H,0);
    dydt(1:n*R)=(1/epsilon_1)*(-x1(:)+v(:));

    % v = SimplexProj(x2 - gradient_2);
    H = H2+delta*eye(R);
    v = max(x2 - gradient_2/H,0);
    dydt(n*R+1:2*n*R)=(1/epsilon_2)*(-x2(:)+ v(:));

    H = H3+delta*eye(R);
    v = max(x3 - gradient_3/H,0);
    % v = SimplexProj(x3-gradient_3);
    dydt(2*n*R+1:3*n*R)=(1/epsilon_3)*(-x3(:)+v(:));
end