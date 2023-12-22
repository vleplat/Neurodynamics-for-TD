function [ Gc,Uc,Vc,Wc,err ] = HOSVD( Y,dim_u,dim_c,random)
%-------------------------------------------------------------------------%
% This function computes the HOSVD of tensor T with approximations R1 R2 R3
% for ranks of unfolding matrices. SVD are can be done by randomized
% approximation.
% 
%-------- Outputs -----%
% - Gc : compressed tensor of size R1xR2xR3
% - Uc,Vc,Wc : unitary transformation matrix, size KxR1, LxR2 and MxR3. Note that R1,R2,R3 are upper bounded by respectively LxM, KxM and KxL 
% - err : error of compression G-Gc where Gc is filled with zero-padding
%
%-------- Inputs ------% 
% - Y : tensor of data
% - dim_u : vector containing the dimensions of Y
% - dim_c : vector containing the dimensions of the compressed space.
% - random : Put 1 for fast but approximate SVDs, and 0 for exact svd.
%------------------------------------------------------------------------

prod_u     =     prod(dim_u)./dim_u;

%--------- Unfolding matrices --------%

Y_mode     =     cell(3);

Y_mode{1}  =     reshape(Y,dim_u(1),prod_u(1));
Y_mode{2}  =     reshape(permute(Y,[2,1,3]),dim_u(2),prod_u(2));
Y_mode{3}  =     reshape(permute(Y,[3,1,2]),dim_u(3),prod_u(3));

%--------- Singular Value Decompositions -------%

switch(random)

    case (0)
        
    [U,~,~]     =     svd(Y_mode{1},'econ');

    [V,~,~]     =     svd(Y_mode{2},'econ');

    [W,~,~]     =     svd(Y_mode{3},'econ');
    
    %--------- Low rank approximation --------%

    Uc    =     U(:,1:dim_c(1));
    Vc    =     V(:,1:dim_c(2));
    Wc    =     W(:,1:dim_c(3));


    case(1)

    [Uc,~,~]=rsvd(Y_mode{1}, dim_c(1), 1);

    [Vc,~,~]=rsvd(Y_mode{2}, dim_c(2), 1);

    [Wc,~,~]=rsvd(Y_mode{3}, dim_c(3), 1);
    
end


%--------- Core tensors computation -------%

%G=MultProd(Y,U',V',W');
Gc=MultProd(Y,Uc',Vc',Wc');

%--------- Approximation error (only for exact svd)-------%

% U1   =     [U(:,1:dim_c(1)),zeros(dim_u(1),min(dim_u(1),prod_u(1))-dim_c(1))];
% 
% V1   =     [V(:,1:dim_c(2)),zeros(dim_u(2),min(dim_u(2),prod_u(2))-dim_c(2))];
% 
% W1   =     [W(:,1:dim_c(3)),zeros(dim_u(3),min(dim_u(3),prod_u(3))-dim_c(3))];
% 
% err  =     norm(vec(MultProd(Y,U1',V1',W1')-G)); 
err=NAN;
end

