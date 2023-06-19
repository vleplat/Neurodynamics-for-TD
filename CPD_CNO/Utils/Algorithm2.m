function [err,cpu_time,A_1,A_2,A_3,lamTab] = Algorithm2(X,options)
%Algorithm2 solves the system of ODE's, aka the three-time scale 
% neurodynamics with
% 1. a full Forward Euler discretization with preconditioning 
% 2. a full Forward Euler discretization with Cubic Regularized Newton
% Methods
% 3. semi-implicit discretization

%-------------------------------------------------------------
% Inputs:
%       - X: the input tensor
%       - options: misc options
%-------------------------------------------------------------
% Outputs:
%       - err : the evolution of relative Frobenius error
%       - cpu_time: the cpu time.
%       - A_i : arrays of cells with factors of the CPD
%       - lamTab : final values for lambda_i
%-------------------------------------------------------------
%--------------------------------------------------------------------------
% Loading parameters
%--------------------------------------------------------------------------
AlgoSel = options.AlgoSel;
preCondiSel = options.preCondiSel;
maxIter = options.maxIter;
alpha = options.alpha;
beta = options.beta;
R = options.R;
initType = options.initType;
delta = options.delta;
tX = tensor(X);

%--------------------------------------------------------------------------
% Init variables
%--------------------------------------------------------------------------
%%% Init. for the factors
switch initType
    case 1
        disp('Random init.')
        % Random
        sizeX = size(X);
        A_10 = rand(sizeX(1),R);
        A_20 = rand(sizeX(2),R);
        A_30 = rand(K,R);
        A_1{1} = A_10;
        A_2{1} = A_20;
        A_3{1} = A_30;

    case 2
        disp('Init. using mixed strategies')
        opts = cp_init();
        opts.init = {'rand' 'fiber' 'fiber' }; %rand fiber fiber
        Ui = cp_init(tensor(X),R,opts);
        A_10 = Ui{1};
        A_20 = Ui{2};
        A_30 = Ui{3};
        A_1{1} = A_10;
        A_2{1} = A_20;
        A_3{1} = A_30;

    case 3
        disp('Init. given')
        A_1{1} = options.U0{1};
        A_2{1} = options.U0{2};
        A_3{1} = options.U0{3};
        
    otherwise
        disp('wrong selection for init method')
        return
end

%%% Errors arrays
err=[];
X_app=(double(full(ktensor(ones(R,1),A_1{1},A_2{1},A_3{1}))));
Error=X-X_app;
relaerror=norm(Error(:))/norm(X(:));
err=[err,relaerror];
tim = tic;

%--------------------------------------------------------------------------
% Main Loop
%--------------------------------------------------------------------------
for i=1:maxIter
    %%% Restart values for parameters lambda_i
    lambda_1=1;
    lambda_2=1;
    lambda_3=1;
   
    %----------------------------------------------------------------------
    %%% Update of factor one
    %----------------------------------------------------------------------
    x3tx3 = A_3{i}'*A_3{i};
    x2tx2 = A_2{i}'*A_2{i};
    H1 = ((x3tx3).*((x2tx2)));
    H=H1+delta*eye(R);
    gradient_1=A_1{i}*H1-mttkrp(tX,{A_1{i},A_2{i},A_3{i}},1);
    switch AlgoSel
        case 1
            if preCondiSel == 1
                P = H;
            else
                P = eye(R);
            end
            A_1{i+1}=A_1{i}+lambda_1*(-A_1{i}+max(A_1{i}-gradient_1/P,0));
        case 2
            [alphak,step] = step_CubReg_Newton(H1,gradient_1);
            A_1{i+1}=A_1{i}+lambda_1*(-A_1{i}+max(A_1{i}-alphak*step,0));
        case 3
            if preCondiSel == 1
                P = H;
            else
                P = eye(R);
            end
            A_1{i+1}=(A_1{i}+lambda_1*(+max(A_1{i}-gradient_1/P,0)))/(1+lambda_1);
        otherwise
            disp('wrong selection for Algorithm 2')
            return
    end
    
    % Line-search
    Error_11=X-double(full(ktensor(ones(R,1),A_1{i},A_2{i},A_3{i})));
    Error_12=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
    
    if (norm(Error_12(:))^2-norm(Error_11(:))^2)>=alpha*lambda_1*gradient_1(:)'*(A_1{i+1}(:)-A_1{i}(:))
        while (norm(Error_12(:))^2-norm(Error_11(:))^2)>alpha*lambda_1*gradient_1(:)'*(A_1{i+1}(:)-A_1{i}(:))
            lambda_1=lambda_1*beta;
            switch AlgoSel
                case 1
                    A_1{i+1}=A_1{i}+lambda_1*(-A_1{i}+max(A_1{i}-gradient_1/P,0));
                case 2
                    [alphak,step] = step_CubReg_Newton(H1,gradient_1);
                    A_1{i+1}=A_1{i}+lambda_1*(-A_1{i}+max(A_1{i}-alphak*step,0));            
                case 3
                    A_1{i+1}=(A_1{i}+lambda_1*(+max(A_1{i}-gradient_1/P,0)))/(1+lambda_1);
                otherwise
                    disp('wrong selection for Algorithm 2')
                    return
            end
            Error_11=X-double(full(ktensor(ones(R,1),A_1{i},A_2{i},A_3{i})));
            Error_12=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
            if lambda_1<1e-9
                break
            end
        end
    end
    
    %----------------------------------------------------------------------
    %%% Update of factor two
    %----------------------------------------------------------------------
    x1tx1 = A_1{i+1}'*A_1{i+1};
    H2 = ((x3tx3).*((x1tx1)));
    H  = H2+delta*eye(R);
    gradient_2 = A_2{i}*H2-mttkrp(tX,{A_1{i+1},A_2{i},A_3{i}},2);
    
    switch AlgoSel
        case 1
            if preCondiSel == 1
                P = H;
            else
                P = eye(R);
            end
            A_2{i+1}=A_2{i}+lambda_2*(-A_2{i}+max(A_2{i}-gradient_2/P,0));
        case 2
            [alphak,step] = step_CubReg_Newton(H2,gradient_2);
            A_2{i+1}=A_2{i}+lambda_2*(-A_2{i}+max(A_2{i}-alphak*step,0));         
        case 3
            if preCondiSel == 1
                P = H;
            else
                P = eye(R);
            end
            A_2{i+1}=(A_2{i}+lambda_2*(+max(A_2{i}-gradient_2/P,0)))/(1+lambda_2);
        otherwise
            disp('wrong selection for Algorithm 2')
            return
    end

    
    % Line-search
    Error_21=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
    Error_22=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i})));
    if (norm(Error_22(:))^2-norm(Error_21(:))^2)>=alpha*lambda_2*gradient_2(:)'*(A_2{i+1}(:)-A_2{i}(:))
        while (norm(Error_22(:))^2-norm(Error_21(:))^2)>alpha*lambda_2*gradient_2(:)'*(A_2{i+1}(:)-A_2{i}(:))
            lambda_2=lambda_2*beta;
            switch AlgoSel
                case 1
                    A_2{i+1}=A_2{i}+lambda_2*(-A_2{i}+max(A_2{i}-gradient_2/P,0));
                case 2
                    [alphak,step] = step_CubReg_Newton(H2,gradient_2);
                    A_2{i+1}=A_2{i}+lambda_2*(-A_2{i}+max(A_2{i}-alphak*step,0));           
                case 3
                    A_2{i+1}=(A_2{i}+lambda_2*(+max(A_2{i}-gradient_2/P,0)))/(1+lambda_2);
                otherwise
                    disp('wrong selection for Algorithm 2')
                    return
            end
            
            Error_21=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i},A_3{i})));
            Error_22=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i})));
            if lambda_2<1e-9
                break
            end
        end
    end
  
    %----------------------------------------------------------------------
    %%% Update of factor three
    %----------------------------------------------------------------------
    x2tx2 = A_2{i+1}'*A_2{i+1};
    H3 = ((x2tx2).*((x1tx1)));
    H = H3+delta*eye(R);
    gradient_3 = A_3{i}*H3-mttkrp(tX,{A_1{i+1},A_2{i+1},A_3{i}},3);

    switch AlgoSel
        case 1
            if preCondiSel == 1
                P = H;
            else
                P = eye(R);
            end
            A_3{i+1}=A_3{i}+lambda_3*(-A_3{i}+max(A_3{i}-gradient_3/P,0));
        case 2
            [alphak,step] = step_CubReg_Newton(H3,gradient_3);
            A_3{i+1}=A_3{i}+lambda_3*(-A_3{i}+max(A_3{i}-alphak*step,0));       
        case 3
            if preCondiSel == 1
                P = H;
            else
                P = eye(R);
            end
            A_3{i+1}=(A_3{i}+lambda_3*(max(A_3{i}-gradient_3/P,0)))/(1+lambda_3);
        otherwise
            disp('wrong selection for Algorithm 2')
            return
    end


    % Line-search
    Error_31=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i})));
    Error_32=X-double(full(ktensor(ones(R,1),A_1{i+1},A_2{i+1},A_3{i+1})));
    if (norm(Error_32(:))^2-norm(Error_31(:))^2)>=alpha*lambda_3*gradient_3(:)'*(A_3{i+1}(:)-A_3{i}(:))
    while (norm(Error_32(:))^2-norm(Error_31(:))^2)>alpha*lambda_3*gradient_3(:)'*(A_3{i+1}(:)-A_3{i}(:))
        lambda_3=lambda_3*beta;
        switch AlgoSel
            case 1
                A_3{i+1}=A_3{i}+lambda_3*(-A_3{i}+max(A_3{i}-gradient_3/P,0));
            case 2
                [alphak,step] = step_CubReg_Newton(H3,gradient_3);
                A_3{i+1}=A_3{i}+lambda_3*(-A_3{i}+max(A_3{i}-alphak*step,0));           
            case 3
                A_3{i+1}=(A_3{i}+lambda_3*(max(A_3{i}-gradient_3/P,0)))/(1+lambda_3);
            otherwise
                disp('wrong selection for Algorithm 2')
                return
        end
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
    if i > 1 && mod(i,50)==0 && options.verbose == 1
        fprintf('Iter: %d   | Rel.Error  %.4d  | min(lam): %.4d  |max(lam): %d \n',i,err(end), min(lamTab),max(lamTab))
    end
end

%--------------------------------------------------------------------------
% Return variables
%--------------------------------------------------------------------------
cpu_time = toc(tim);
end