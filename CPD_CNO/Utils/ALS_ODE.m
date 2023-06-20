function [err_ode,err_ode2,cpu_time,B] = ALS_ODE(X,epsilon,theta,tspan,options)
%ALS_ODE computes the nonnegative CPD of an input tensor X
%-------------------------------------------------------------
% Inputs:
%       - X: the input tensor
%       - epsilon: the time constant for n-scale neurodynamics
%       - theta: the init.
%       - tspan: time interval for ODE solver
%       - options: misc options
%-------------------------------------------------------------
% Outputs:
%       - err_ode: 
%       - err_ode2:
%       - cpu_time: the cpu time.
%       - B: cell of arrays with the factors of the CPD
%-------------------------------------------------------------

% Loading parameters
R = options.R;
maxKrun = options.maxKrun;
algo_Sel = options.algo_Sel;

% Init variables
err_ode = [];
normX = norm(X(:));
tim = tic;

% Main Loop
for krun = 1:maxKrun
    [~,yode] = ode45(@(t,xx) ODE_sol(t,xx,X,epsilon,R,algo_Sel), tspan, theta,options);

    clear err_ode2;
    for k = 1:size(yode,1)
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
     
end

% Return variables

cpu_time = toc(tim);
B = Bnew;

end