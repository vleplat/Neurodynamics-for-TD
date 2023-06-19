function dydt=ODE_sol(t,xx,X,epsilon,R,algo_Sel)
%% Function used to solve the system of ODE's

%size of dynamic system
szX = size(X);
dydt=zeros(sum(szX)*R,1);
B = vec2fac(xx,szX);

% Choice of the step 1 for algorithm (before calling ODE solver)
% alg = 'als2';
% alg = 'hals2';
alg = algo_Sel;

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

        x1_als = max(tX1/(H1+delta*eye(R)),0);% x1als = max(0,Y1*khatr(B3,B2)*inv(H1))
        x2_als = max(tX2/(H2+delta*eye(R)),0);
        x3_als = max(tX3/(H3+delta*eye(R)),0);

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

        opts = ncp_hals;
        opts.init = B;
        opts.maxiters = 1;
        opts.tol = 1e-10;
        [Yx,~] = ncp_hals(tensor(X,[1 size(X)]),R,opts);
        %[Yx,out] = ncp_fLM(tensor(X),R,opts);

        Bnew = Yx.U;
        Bnew{1} = Bnew{1}*diag(Yx.lambda);
        dB = B;
        for k = 1:numel(Bnew)
            dB{k} = (Bnew{k}-B{k})/epsilon.eps_1;
        end
        dB = dB([2:end 1]);
        dydt = fac2vec(dB);


end


end